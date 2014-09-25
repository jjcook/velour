//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

// T O D O: copyright statement and GPL (or whatever) license header

//
// Velour -- memory efficient de Bruijn graph-based short read de novo
//             DNA assembly front-end tool
//
//
// Authors:  Jeffrey J. Cook      jeffrey.j.cook@gmail.com
//           Craig Zilles         zilles@illinois.edu
//

#include "types.h"

#include "minikmer.h"  // XXX remove this

//
// global variables
//

#ifdef VELOUR_TBB
static int p__TBB_THREAD_LIMIT = tbb::task_scheduler_init::automatic;

static tbb::atomic<unsigned> thread_next_index;
__thread unsigned tls_thread_index = 0;

static tbb::cache_aligned_allocator<bool> abort_allocator;
bool *thread_aborts[THREADID_MAX+2] = {NULL}; // TODO: make this an atomic vector?
__thread bool *tls_thread_abort = NULL;

__thread PrivateSeqGraphSet *tls_private_sgraph_set = NULL;

__thread TLS_KmerNodeAllocator *tls_kmernode_allocator = NULL;
__thread TLS_SeqNodeAllocator *tls_seqnode_allocator = NULL;

class thread_index_observer : public tbb::task_scheduler_observer
{
    public:
        thread_index_observer() : tbb::task_scheduler_observer() { observe(true); }
        virtual void on_scheduler_entry( bool is_worker )
        {
            assert(tls_thread_index == 0);
            tls_thread_index = thread_next_index++;
            assert( thread_next_index <= THREADID_MAX );
            printf("new thread appears, with id: %u\n", tls_thread_index);

            assert(thread_aborts[tls_thread_index] == NULL);
            thread_aborts[tls_thread_index] = abort_allocator.allocate(sizeof(bool));
            assert(tls_thread_abort == NULL);
            tls_thread_abort = thread_aborts[tls_thread_index];
            *tls_thread_abort = false;

            assert(tls_private_sgraph_set == NULL);
            tls_private_sgraph_set = new PrivateSeqGraphSet();

            tls_kmernode_allocator = g__KMERNODE_ALLOCATOR->GetTlsAllocator(tls_thread_index);
            tls_seqnode_allocator = g__SEQNODE_ALLOCATOR->GetTlsAllocator(tls_thread_index);
        }
};

#endif // VELOUR_TBB

//
// command-line parameters and derived values
//

static bool p__CLOBBER_WORK_DIRECTORY = true;

static bool        p__SINGLEPASS           = true;  // default: directly emit pregraph
static char *      p__PREGRAPH_FILENAME    = NULL;

static bool        p__PGPART               = false;
static bool        p__PGDIST               = false;
unsigned           g__PGDIST_PARTITIONS    = 0;
unsigned           g__PGDIST_FILTER        = 0;

static bool        p__LOOMING              = false;
static bool        p__QUILTING             = false;
static bool        p__FLOW                 = false;

//
// static functions
//

// Print out program usage information
static void
printUsage() {
  puts("Usage:");
  puts("./velour directory kmer_length <options> {[-file_format][-read_type] filenames}");
  puts("     directory             = directory name for output files");
  puts("     kmer_length           = odd integer (if even, it will be decremented)");
  puts("     file_format           = fasta | fastq | fasta.gz | fastq.gz | eland | gerald ");
  puts("     read_type             = short | shortPaired | short2 | shortPaired2 | long");  // FIXME: we don't use this
  puts("     filenames             = name of input files or directory");
  puts("   -drop                 = perform single-copy kmer removal");
  puts("   -noclip               = do not perform tip clipping");
  puts("   -stat                 = enables full statistics reporting (slow!)");
  puts("");
  puts("  Experimental options:");
  puts("    -part partitions minikmer_length  = partition to reduce max memory usage");
  puts("         partitions                       = power-of-4 integer (if not, will be rounded up)");
  puts("         minikmer_length                  = odd integer <= kmer_length");
  puts("       -noemit                            = don't emit loom files (just emit statistics)");
  puts("    -loom");
  puts("    -flow");
  puts("    -quilt");
  puts("    -break");
  puts("");
  puts("Outputs:");
  puts("  - directory/PreGraph");
  puts("");
}

unsigned parseUnsigned(char *h); // forward decl

void initializePartitionIndexFromInputFilename(file_object_vector *file_objects); // forward decl

// parse kmer_length from string
static inline unsigned
getKmerLength(char *h)
{
    unsigned kmer_length = parseUnsigned(h);

    if (kmer_length > MAXKMERLENGTH) {
        fprintf(stderr,"ERROR: k-mer length of %u is larger than MAXKMERLENGTH of %u\n"
                       "       please recompile and increase MAXKMERLENGTH in Makefile.\n", kmer_length, MAXKMERLENGTH);
        exit(EXIT_FAILURE);
    } else if (kmer_length % 2 == 0) {
        fprintf(stderr,"ERROR: Velour requires an odd-length k-mer size, %u is even.\n", kmer_length);
        exit(EXIT_FAILURE);
    }
    return kmer_length;
}

// parse minikmer length from string
static inline unsigned
getMiniKmerLength(char *h)
{
    unsigned minikmer_length = parseUnsigned(h);

#ifdef PART_MINIKMER_BITS
    unsigned limit = g__FULLKMER_LENGTH * 2;
    if (minikmer_length > limit) {
        fprintf(stderr,"ERROR: mini k-mer bits of %u cannot be larger than double the full k-mer length %u\n",
            minikmer_length, limit);
        exit(EXIT_FAILURE);
    }
#else
    unsigned limit = g__FULLKMER_LENGTH;
    if (minikmer_length > limit) {
        fprintf(stderr,"ERROR: mini k-mer of %u cannot be larger than the full k-mer length %u\n",
            minikmer_length, limit);
        exit(EXIT_FAILURE);
    }
#endif

#ifdef PART_MIDDLE_MINIKMER
    if (minikmer_length % 2 == 0) {
        fprintf(stderr,"ERROR: mini k-mers are required to be an odd-length, %u is even.\n", minikmer_length);
        exit(EXIT_FAILURE);
    }
#endif

    return minikmer_length;
}

// parse partition count from string and correct (round-up) to power-of-4
static inline unsigned
getNumPartitions(char *h)
{
  unsigned numPartitions = parseUnsigned(h);

  // FIXME: is this only for metis?
#ifdef PART_POWEROF4
  if( !is_power_of_4(numPartitions) )
  {
    unsigned tmp = round_up_power_of_2(numPartitions);
    unsigned pow4 = tmp & 0x55555555UL ? tmp : tmp << 1;

    fprintf(stderr,"WARNING: partitions specified %u is not a power of 4.  Using %u instead.\n",
        numPartitions, pow4);
    numPartitions = pow4;
  }
#endif

  unsigned maxPartitions = (1 << (8 * sizeof(color_t))) - 2;

#ifdef PART_POWEROF4
    if (numPartitions > maxPartitions) {
        unsigned tmp = numPartitions >> 2;
        fprintf(stderr,"WARNING: partitions %u is too large.  Using %u instead.\n", numPartitions, tmp);
        numPartitions = tmp;
    }
#else
    if (numPartitions > maxPartitions) {
        fprintf(stderr,"WARNING: partitions %u is too large.  Using %u instead.\n", numPartitions, maxPartitions);
        numPartitions = maxPartitions;
    }
#endif

  return numPartitions;
}

// parse the command line, building a vector of descriptors for the files to be processed.
static file_object_vector *
parseCommandLine(int argc, char **argv)
{
  assert( argc >= 4 );

  g__WORK_BASE_DIRECTORY = argv[1];
  g__WORK_LOOM_DIRECTORY = argv[1]; // TODO FIXME
  g__WORK_QUILT_DIRECTORY = argv[1]; // TODO FIXME
  g__WORK_INBOX_ROOT_DIRECTORY = argv[1]; // TODO FIXME

  g__FULLKMER_LENGTH = getKmerLength(argv[2]);

  p__PREGRAPH_FILENAME = static_cast<char *>( calloc((PATH_MAX+1), sizeof(char)) );
  strcpy(p__PREGRAPH_FILENAME, g__WORK_BASE_DIRECTORY);
  strcat(p__PREGRAPH_FILENAME, "/PreGraph");

  file_object_t current_entry;
  file_object_vector *work_to_do = new file_object_vector;

  current_entry.filetype = FASTA;
  current_entry.fileindex = 0;
  current_entry.cat = 0;
  current_entry.length = -1;

  for (int argIndex = 3; argIndex < argc; argIndex++) {
    if (argv[argIndex][0] == '-') {
      // FILE FORMAT QUALIFIERS
      if (strcmp(argv[argIndex], "-fastq") == 0) {
		  current_entry.filetype = FASTQ;
          fprintf(stderr,"ERROR: Velour read format FASTQ not yet supported.\n");
          exit(EXIT_FAILURE);
      } else if (strcmp(argv[argIndex], "-fasta") == 0) {
		  current_entry.filetype = FASTA;
      } else if (strcmp(argv[argIndex], "-gerald") == 0) {
		  current_entry.filetype = GERALD;
          fprintf(stderr,"ERROR: Velour read format GERALD not yet supported.\n");
          exit(EXIT_FAILURE);
      } else if (strcmp(argv[argIndex], "-eland") == 0) {
		  current_entry.filetype = ELAND;
          fprintf(stderr,"ERROR: Velour read format ELAND not yet supported.\n");
          exit(EXIT_FAILURE);
      } else if (strcmp(argv[argIndex], "-fastq.gz") == 0) {
		  current_entry.filetype = FASTQ_GZ;
          fprintf(stderr,"ERROR: Velour read format FASTQ_GZ not yet supported.\n");
          exit(EXIT_FAILURE);
      } else if (strcmp(argv[argIndex], "-fasta.gz") == 0) {
		  current_entry.filetype = FASTA_GZ;
          fprintf(stderr,"ERROR: Velour read format FASTA_GZ not yet supported.\n");
          exit(EXIT_FAILURE);
      } else if (strcmp(argv[argIndex], "-maq.gz") == 0) {
		  current_entry.filetype = MAQ_GZ;
          fprintf(stderr,"ERROR: Velour read format MAQ_GZ not yet supported.\n");
          exit(EXIT_FAILURE);
      }
      // CATEGORY QUALIFIERS
      else if (strcmp(argv[argIndex], "-short") == 0) {
		  current_entry.cat = 0;
      } else if (strcmp(argv[argIndex], "-shortPaired") == 0) {
		  current_entry.cat = 1;
      } else if (strncmp(argv[argIndex], "-shortPaired", 12) == 0) { // must come before the -shortX check
          short short_var;
          sscanf(argv[argIndex], "-shortPaired%hd", &short_var);
		  current_entry.cat = short_var;
          if (current_entry.cat < 1 || current_entry.cat > CATEGORIES) {
            printf("Invalid option integer (%hd): %s\n", current_entry.cat, argv[argIndex]);
            exit(EXIT_FAILURE);
          }
          current_entry.cat --;
          current_entry.cat *= 2;
          current_entry.cat ++;
      }	else if (strncmp(argv[argIndex], "-short", 6) == 0) {
          short short_var;
          sscanf(argv[argIndex], "-short%hd", &short_var);
		  current_entry.cat = short_var;
          if (current_entry.cat < 1 || current_entry.cat > CATEGORIES) {
            printf("Invalid option integer (%hd): %s\n", current_entry.cat, argv[argIndex]);
            exit(EXIT_FAILURE);
          }
          current_entry.cat --;
          current_entry.cat *= 2;
      } else if (strcmp(argv[argIndex], "-long") == 0) {
		  current_entry.cat = CATEGORIES * 2;
      } else if (strcmp(argv[argIndex], "-longPaired") == 0) {
		  current_entry.cat = (CATEGORIES * 2) + 1;
      }
      // GENERAL OPTIONS
      else if (strcmp(argv[argIndex], "-stat") == 0) {
            g__FULL_STATISTICS = true;
      } else if (strcmp(argv[argIndex], "-mem") == 0) {
#ifdef ARCH_64BIT
            unsigned megabytes = parseUnsigned(argv[++argIndex]);
            g__MEMORY_FOOTPRINT_LIMIT = megabytes * 1024 * 1024ULL;
#else // 32-bit
            unsigned megabytes = parseUnsigned(argv[++argIndex]);
            if (megabytes >= 4096) {
                fprintf(stderr,"WARNING: Memory footprint limit of %u MB is too big, as Velour was built in 32-bit mode...\n", megabytes);
                megabytes = 4000;
                fprintf(stderr,"WARNING:   using %u MB instead.\n", megabytes);
            }
            g__MEMORY_FOOTPRINT_LIMIT = megabytes * 1024 * 1024UL;
#endif
#ifdef VELOUR_TBB
      } else if (strcmp(argv[argIndex], "-thr") == 0) {
            p__TBB_THREAD_LIMIT = parseUnsigned(argv[++argIndex]);
#endif // VELOUR_TBB
      }
      // PARTITIONING SPECIFIC
      else if (strcmp(argv[argIndex], "-part") == 0) {
            g__PARTITIONING = true;
			p__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            g__PARTITION_COUNT = getNumPartitions(argv[++argIndex]);
            g__MINIKMER_LENGTH = getMiniKmerLength(argv[++argIndex]);
      } else if (strcmp(argv[argIndex], "-noemit") == 0) {
            g__PARTITION_NOEMIT = true;
      } else if (strcmp(argv[argIndex], "-minfootprint") == 0) {
            g__MINIMIZE_FOOTPRINT = true;
      } else if (strcmp(argv[argIndex], "-loom") == 0) {
            p__LOOMING = true;
			p__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = LOOM;
      } else if (strcmp(argv[argIndex], "-quilt") == 0) {
            p__QUILTING = true;
			p__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = QUILT;
      } else if (strcmp(argv[argIndex], "-bucket") == 0) {
            current_entry.filetype = BUCKET;
      } else if (strcmp(argv[argIndex], "-flow") == 0) {
            p__FLOW = true;
			p__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = LOOM;
            g__PARTITION_COUNT = getNumPartitions(argv[++argIndex]);
      } else if (strcmp(argv[argIndex], "-slice") == 0) {
            g__SLICING = true;
      } else if (strcmp(argv[argIndex], "-noslice") == 0) {
            g__SLICING = false;
      } else if (strcmp(argv[argIndex], "-noclip") == 0) {
            g__NO_TIP_CLIPPING = true;
      } else if (strcmp(argv[argIndex], "-pgpart") == 0) {
            p__PGPART = true;
      } else if (strcmp(argv[argIndex], "-pgdist") == 0) {
            p__PGDIST = true;
			p__SINGLEPASS = false;
            g__PGDIST_PARTITIONS = parseUnsigned(argv[++argIndex]);
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = QUILT;
      } else if (strcmp(argv[argIndex], "-pgfilter") == 0) {
            g__PGDIST_FILTER = parseUnsigned(argv[++argIndex]);
            assert( g__PGDIST_FILTER < 38 ); // FIXME
      } else if (strcmp(argv[argIndex], "-cov_cutoff") == 0) {
            g__COVCUTOFF_MIN = atof(argv[++argIndex]);
      } else if (strcmp(argv[argIndex], "-max_coverage") == 0) {
            g__COVCUTOFF_MAX = atof(argv[++argIndex]);
      } else if (strcmp(argv[argIndex], "-bubble_removal") == 0) {
            g__BUBBLE_POPPING = true;
      }
	  else {
          fprintf(stderr, "ERROR: Unknown option: '%s'\n", argv[argIndex]);
          exit(EXIT_FAILURE);
      }

      continue;
    }

    // input file(s): test if its a directory or a filename
    DIR *dir = opendir(argv[argIndex]);
    if (dir != NULL) { // TODO: support subdirectories???
      // if a directory, add all regular files to the input list
      if( p__LOOMING ) {
          fprintf(stderr, "ERROR: Can't '-loom' a directory (yet).  Exiting...\n");
          exit(1);
      }
      struct dirent * ep;
      struct stat st;
      while ((ep = readdir(dir)) != NULL) {
        char filename[PATH_MAX+1];
        strcpy(filename, argv[argIndex]);

        int pathlen = strlen(filename);
        if(filename[pathlen-1] != '/' ) { filename[pathlen] = '/'; filename[pathlen+1] = '\0'; }

        strcat(filename, ep->d_name);

        if (ep->d_type == DT_REG || ep->d_type == DT_LNK ||
            (ep->d_type == DT_UNKNOWN && !stat(filename, &st) && (S_ISREG(st.st_mode) || S_ISLNK(st.st_mode))) )
        {
          char *work_filename = (char *) malloc((PATH_MAX+1) * sizeof(char));
          strcpy(work_filename, filename);

          printf("Using directory input file: %s\n", work_filename);

          current_entry.filename = work_filename;
          work_to_do->push_back(current_entry);
          current_entry.fileindex++;
        }
      }
      closedir(dir);
    } else {
      // else it's a filename
      current_entry.filename = argv[argIndex];
      work_to_do->push_back(current_entry);
      current_entry.fileindex++;
    }
  }

  // check input option combinations
  if( g__PARTITIONING + p__LOOMING + p__FLOW + p__QUILTING > 1 ) {
      fprintf(stderr, "Invalid combination of '-part' '-loom' '-flow' and/or '-quilt'.  Exiting...\n");
      exit(EXIT_FAILURE);
  }

  if (g__COVCUTOFF_MAX != 0.0 && g__COVCUTOFF_MIN >= g__COVCUTOFF_MAX) {
      fprintf(stderr, "Invalid combination of -cov_cutoff and -max_coverage -- Exiting...\n");
      exit(EXIT_FAILURE);
  }

  return work_to_do;
}


// writes command line arguments to log file
static void
logInstructions(int argc, char **argv, char *directory) {
  int index;
  char *logFilename = (char *) malloc((PATH_MAX+1) * sizeof(char));
  FILE *logFile;
  time_t date;
  char *string;

  time(&date);
  string = ctime(&date);

  strcpy(logFilename, directory);
  strcat(logFilename, "/Log");
  logFile = fopen(logFilename, "a");

  if (logFile == NULL) {
    printf("Could not open file %s, exiting...\n", logFilename);
    exit(EXIT_FAILURE);
  }

  fprintf(logFile, "%s", string);

  for (index = 0; index < argc; index++) {
    fprintf(logFile, " %s", argv[index]);
  }

  fprintf(logFile, "\n");

  fclose(logFile);
  free(logFilename);
}

// create the directory if it doesn't exist; if it exists blow away any old files
void
makeDirectories(char *directory) {
  DIR *dir = opendir(directory);
  if (dir == NULL) {  // TODO: create loom/etc subdirectories
    mkdir(directory, 0777);
  } else if( p__CLOBBER_WORK_DIRECTORY ) { // TODO FIXME: work subdirectory clobbering
    assert( PATH_MAX > strlen(directory) + NAME_MAX );
    char *buf = (char *) malloc((PATH_MAX+1) * sizeof(char));

    struct dirent * ep;
    while((ep = readdir(dir)) != NULL)
      if( !strncmp(ep->d_name, "PreGraph", 8) ||
          !strncmp(ep->d_name, "Graph", 5)    ||
          !strncmp(ep->d_name, "Subsequences", 12)           ) // TODO
      {
        sprintf(buf, "%s/%s", directory, ep->d_name);
        remove(buf);
      }

    sprintf(buf, "%s/Log", directory);
    remove(buf);

    free(buf);
    closedir(dir); // TODO: close even if no clobber
  }
}

void printLinuxProcStatus(int signal)
{
    fflush(stderr);
    fflush(stdout);
    char vmInstrumentation[100];
    sprintf(vmInstrumentation, "cat /proc/%"PRIdMAX"/status | grep -i vm | perl -e 'while(<>){print \"(LINUX) $_\";}'", static_cast<intmax_t>(getpid()));
    system(vmInstrumentation);
    fflush(stdout);
    fflush(stderr);
    if (signal > 0) {
        exit(signal);
    }
}


void runPartitioning(file_object_vector *file_objects)
{
    initializePartitionerPreObservation();

#ifdef PART_TRAINS_ON_INPUT
    g__PHASE = PHASE_OBSERVE_SEQUENCES;
    process_files(*file_objects);
    completeObservation();
#endif // PART_TRAINS_ON_INPUT

    initializePartitionerPostObservation();

    char velvet_sequences_filename[PATH_MAX];
    sprintf(velvet_sequences_filename, "%s/Sequences", g__WORK_BASE_DIRECTORY);
    /*g__VELVET_SEQUENCES_FILE = fopen(velvet_sequences_filename, "w");
    if (g__VELVET_SEQUENCES_FILE == NULL) {
        printf("Could not create file %s, exiting...\n", velvet_sequences_filename);
        exit(EXIT_FAILURE);
    }*/

    g__PHASE = PHASE_DISTRIBUTE_SEQUENCES;
    process_files(*file_objects);

    if (g__VELVET_SEQUENCES_FILE != NULL) { fclose(g__VELVET_SEQUENCES_FILE); }
    g__VELVET_SEQUENCES_FILE = NULL;

    if( g__FULL_STATISTICS ) {
        printf("Generating detailed statistics...\n");
        printMiniKmerStats();
    }
    printPartitionerStatistics();

    destroyPartitioner();

    // communicate number of actual partitions back to environment
    char filename[PATH_MAX+1];
    sprintf(filename, "%s/common.partitions", g__WORK_BASE_DIRECTORY);
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        fprintf(stderr, "ERROR: fopen()\n");
        assert(false);
        exit(EXIT_FAILURE);
    }
    fprintf(f, "%u", g__PARTITION_COUNT);
    fclose(f);
}

void runLooming(file_object_vector *file_objects)
{
    initializePartitionIndexFromInputFilename(file_objects);

    // next, load the subsequences into kmer graph
    load_loom_files(*file_objects);

    printf("%" PRIuPTR " kmer nodes built from Loom file(s)\n", g__KG_HASHTABLE->node_count);

    if( g__FULL_STATISTICS ) {
        kg_stat_components(g__KG_HASHTABLE, stdout);
    }

    // then, optimize the kmer graph
    remove_tips(g__KG_HASHTABLE);

    // last, concatenate directly to a quilt file
    char quilt_filename[PATH_MAX+1];
    sprintf(quilt_filename, "%s/Partition-%u.quilt", g__WORK_QUILT_DIRECTORY, g__PARTITION_INDEX);
    sg_dump_quilt_from_kmergraph(g__KG_HASHTABLE, quilt_filename);
}

void runQuilt(file_object_vector *file_objects)
{
    g__PARTITION_INDEX = UINT_MAX; // TODO: should be g__PARTITION_COUNT+1
    quilt_files(g__SG_HASHTABLE, *file_objects);

    // initialize the read count for pregraph emittal
    char filename[PATH_MAX+1];
    sprintf(filename, "%s/common.sequences", g__WORK_BASE_DIRECTORY);
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "ERROR: fopen()\n");
        assert(false);
        exit(EXIT_FAILURE);
    }
    fscanf(f, "%"SCNu64"", &g__READ_COUNT);
    fclose(f);

    sg_dump_pregraph(g__SG_HASHTABLE, p__PREGRAPH_FILENAME);

    if (p__PGPART) {
        char metis_filename[PATH_MAX+1];
        sprintf(metis_filename, "%s/PreGraph.metis", g__WORK_BASE_DIRECTORY);
        pregraph_partitioning(g__SG_HASHTABLE, metis_filename);

        char quilt_filename[PATH_MAX+1];
        sprintf(quilt_filename, "%s/PreGraph.quilt", g__WORK_BASE_DIRECTORY);
        sg_dump_quilt(g__SG_HASHTABLE, quilt_filename); // TODO FIXME: this only works if 16-bit colors alias vid!
    }
}

void runDirect(file_object_vector *file_objects)
{
    char velvet_sequences_filename[PATH_MAX];
    sprintf(velvet_sequences_filename, "%s/Sequences", g__WORK_BASE_DIRECTORY);
    /*g__VELVET_SEQUENCES_FILE = fopen(velvet_sequences_filename, "w");
    if (g__VELVET_SEQUENCES_FILE == NULL) {
        printf("Could not create file %s, exiting...\n", velvet_sequences_filename);
        exit(EXIT_FAILURE);
    }*/

    process_files(*file_objects);  // FIXME: TBB version

    if (g__VELVET_SEQUENCES_FILE != NULL) { fclose(g__VELVET_SEQUENCES_FILE); }
    g__VELVET_SEQUENCES_FILE = NULL;

    printf("%" PRIuPTR " kmer nodes built from input file(s)\n", g__KG_HASHTABLE->node_count); fflush(stdout);

    //print_hashtable_histogram(g__KG_HASHTABLE);

    if (g__DOTGRAPH) {
      char dot_filename[PATH_MAX+1];
      sprintf(dot_filename, "%s/Direct-initial-kmergraph.dot", g__WORK_BASE_DIRECTORY);
      emit_graphviz(g__KG_HASHTABLE, dot_filename);
    }

    if (g__FULL_STATISTICS) {
        //print_hashtable_histogram(g__KG_HASHTABLE);
        kg_stat_components(g__KG_HASHTABLE, stdout);
    }

    remove_tips(g__KG_HASHTABLE);  // FIXME:  TBB version

    if (g__DOTGRAPH) {
      char dot_filename[PATH_MAX+1];
      sprintf(dot_filename, "%s/Direct-clipped-kmergraph.dot", g__WORK_BASE_DIRECTORY);
      emit_graphviz(g__KG_HASHTABLE, dot_filename);
    }

    if (p__PGPART) {
        char quilt_filename[PATH_MAX+1];
        sprintf(quilt_filename, "%s/Direct.quilt", g__WORK_QUILT_DIRECTORY);
        sg_dump_quilt_from_kmergraph(g__KG_HASHTABLE, quilt_filename);
    } else {
        sg_dump_pregraph_from_kmergraph(g__KG_HASHTABLE, p__PREGRAPH_FILENAME);
    }
}

void runFlow(file_object_vector *file_objects, KmerGraph *resident_kgraph, SeqGraph *resident_sgraph); // forward decl

void runPregraphDistribution(file_object_vector *file_objects, SeqGraph *sgraph); // forward decl

//********************************************************************************
//****************************    Main Program    ********************************
//********************************************************************************
int
main(int argc, char **argv)
{
    struct timeval wallclock_start_time;
    gettimeofday(&wallclock_start_time, NULL);

/*
#ifdef VELOUR_TBB
    tbb::tick_count main_start_time = tbb::tick_count::now();
#endif
*/

  // parse the command line
  //
  if (argc < 4) {
    puts("velour - sequence assembly front-end");
    printf("Version %i.%i.%2.2i\n", VERSION_NUMBER, RELEASE_NUMBER, UPDATE_NUMBER);
    printUsage();
    exit(EXIT_FAILURE);
  }
  file_object_vector *file_objects = parseCommandLine(argc, argv);

  if (g__MEMORY_FOOTPRINT_LIMIT == 0) { // not user specified
#ifdef ARCH_64BIT
      g__MEMORY_FOOTPRINT_LIMIT = 2 * 1024 * 1024 * 1024ULL; // TODO: base on machine physical memory
#else // 32-bit
      g__MEMORY_FOOTPRINT_LIMIT = 2 * 1024 * 1024 * 1024UL; // TODO: base on machine physical memory
#endif
  }

#ifdef VERIFY
  printf("VERIFY enabled!\n");
#endif

  makeDirectories(g__WORK_BASE_DIRECTORY);
  logInstructions(argc, argv, g__WORK_BASE_DIRECTORY);

  // initialization
  //

  // NOTE: node allocators must be constructed before TBB task scheduler initialization
  g__NODE_ALLOCATORS = new NodeAllocators((size_t)((g__MEMORY_FOOTPRINT_LIMIT * 0.90)) /* TODO */);
  g__KMERNODE_ALLOCATOR = &( g__NODE_ALLOCATORS->get_kmernode_allocator() );
  g__SEQNODE_ALLOCATOR = &( g__NODE_ALLOCATORS->get_seqnode_allocator() );

#ifdef VELOUR_TBB
  printf("Using TBB.\n");
#  ifdef VERIFY
    const int max_threads = 1; // XXX: for now, single thread when verifying
#  else
    const int max_threads = p__TBB_THREAD_LIMIT;
#  endif
  tbb::task_scheduler_init init( max_threads );
  tbb::tick_count time0, time1;
  thread_next_index = 1;
  thread_index_observer tbb_tid_observer;
#endif
#ifdef USE_TBB_ALLOC
  printf("Using TBB allocation.\n");
#else
  printf("Using internal allocation.\n");
  velour_alloc_init();
#endif

#ifdef RANDOMSEED
  printf("Using srandom seed = %d\n",RANDOMSEED);
  srandom(RANDOMSEED);
#else
  FILE *urandom = fopen("/dev/urandom","r");
  if( urandom == NULL )
  {
    fprintf(stderr,"Open of /dev/urandom failed... exiting.\n");
    exit(EXIT_FAILURE);
  }
  unsigned long theseed;
  fread(&theseed, 1, sizeof(unsigned long), urandom);
  fclose(urandom);
  srandom(theseed);
#endif

  initializeBaseMap();

  printf("  kmer node size %zu\n", sizeof(KmerNode));
  printf("  sequence node size %zu\n", sizeof(SeqNode) + Sequence::MIN_BYTES);

  //
  // for now, choose the hashtable size based on the footprint limit:
  //   -- using 1/32 (about 3%) of footprint limit for each hashtable
  //
  uintptr_t buckets = 0;
  if (is_power_of_2(g__MEMORY_FOOTPRINT_LIMIT)) {
    buckets = (g__MEMORY_FOOTPRINT_LIMIT >> 5) / sizeof(void *);
  } else {
    size_t temp = round_up_power_of_2(g__MEMORY_FOOTPRINT_LIMIT) >> 1; // round down power of 2
    buckets = (temp >> 5) / sizeof(void *);
  }
  assert( is_power_of_2(buckets) );
  printf("%" PRIuPTR " hashtable buckets\n", buckets);

  if (!g__PARTITIONING) {
    if( p__SINGLEPASS || p__LOOMING ) {
        g__KG_HASHTABLE = new KmerGraph(buckets);
    } else if (p__FLOW) {
        g__KG_HASHTABLE = new KmerGraph(buckets);
        g__SG_HASHTABLE = new SeqGraph(buckets);
    } else {
        g__SG_HASHTABLE = new SeqGraph(buckets);
    }
  }
  g__KMERNODE_ALLOCATOR->set_kmergraph(g__KG_HASHTABLE);
  g__SEQNODE_ALLOCATOR->set_seqgraph(g__SG_HASHTABLE);

  pid_t velour_pid = getpid();
  printf("process id %" PRIdMAX "\n", static_cast<intmax_t>(velour_pid));

    // before&after /proc/vmstat provides vm swapping data
    /*fflush(stdout);
    fflush(stderr);
    char vmstatInstrumentation[100];
    sprintf(vmstatInstrumentation, "cat /proc/vmstat | perl -e 'while(<>){print \"(PROC) $_\";}'");
    system(vmstatInstrumentation);
    fflush(stdout);
    fflush(stderr);*/

/*#ifdef __linux__
    signal(SIGINT, &printLinuxProcStatus);
#endif*/

    //
    // PRIMARY FRONT-END ASSEMBLY CONTROL SEQUENCE
    //

    if (g__PARTITIONING) {
        runPartitioning(file_objects);
    } else if (p__LOOMING) {
        runLooming(file_objects);
    } else if (p__QUILTING) {
        runQuilt(file_objects);
    } else if (p__FLOW) {
        runFlow(file_objects, g__KG_HASHTABLE, g__SG_HASHTABLE);
    } else if (p__PGDIST) {
        runPregraphDistribution(file_objects, g__SG_HASHTABLE);
    } else {
        runDirect(file_objects);
    }

#ifdef VELOUR_TBB
    printf("  %u total TBB threads\n", thread_next_index-1);
#endif

/*
#ifdef VELOUR_TBB
    tbb::tick_count main_stop_time = tbb::tick_count::now();
    tbb::tick_count::interval_t total_time = main_stop_time - main_start_time;
    printf("%lld seconds TBB TIME\n", (long long int)total_time.seconds());
#endif
*/

    struct timeval wallclock_stop_time;
    gettimeofday(&wallclock_stop_time, NULL);
    printf("%lld seconds WALL CLOCK TIME\n", (long long int)(wallclock_stop_time.tv_sec - wallclock_start_time.tv_sec));

#ifndef USE_TBB_ALLOC
    velour_alloc_done();
#endif

    delete g__NODE_ALLOCATORS;

#ifdef __linux__
    printLinuxProcStatus(0);
#endif

    // printf("Velour run completed successfully.\n");
    return 0;
}

