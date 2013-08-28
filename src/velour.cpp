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

// TODO: make these private variables
//
bool        g__SINGLEPASS           = true;  // default: directly emit pregraph
char *      g__PREGRAPH_FILENAME    = NULL;

bool        g__PGPART               = false;
bool        g__PGDIST               = false;
unsigned    g__PGDIST_PARTITIONS    = 0;
unsigned    g__PGDIST_FILTER        = 0;

bool        g__LOOMING              = false;
bool        g__SPLIT                = false;
bool        g__ISPLIT               = false;
bool        g__IESPLIT              = false;
bool        g__CIESPLIT             = false;
bool        g__QUILTING             = false;
bool        g__FLOW                 = false;

//
// static functions
//

// Print out program usage information 
static void 
printUsage() {
  puts("Usage:");
  puts("./velour directory kmer_length <options> {[-file_format][-read_type] filenames}");
  puts("     directory             = directory name for output files");
  puts("     kmer_length           = odd integer (if even, it will be decremented) <= 31 (if above, will be reduced)");
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
  puts("    -quilt");
  puts("    -break");
  puts("");
  puts("Outputs:");
  puts("  - directory/PreGraph");
  puts("");
}

unsigned parseUnsigned(char *h); // forward decl

void initializePartitionIndexFromInputFilename(file_object_vector *file_objects); // forward decl

// parse kmer_length from string and correct to make odd and less than 32 bases long
static inline unsigned
getKmerLength(char *h)
{
  unsigned kmer_length = parseUnsigned(h);

  if (kmer_length > 31) {
    fprintf(stderr,"WARNING: Velour can't handle k-mers as long as %u!  Using k-mer length of 31 instead.\n",
	   kmer_length);
    kmer_length = 31;
  } else if (kmer_length % 2 == 0) {
    fprintf(stderr,"WARNING: Velour can't work with even length k-mers, such as %u.  Using k-mer length of %u instead.\n",
	   kmer_length, kmer_length - 1);
    kmer_length--;
  }
  return kmer_length;
}

// parse minikmer length from string and correct it as necessary
static inline unsigned
getMiniKmerLength(char *h)
{
    unsigned minikmer_length = parseUnsigned(h);

#ifdef PART_MINIKMER_BITS
    unsigned limit = g__FULLKMER_LENGTH * 2;
    if (minikmer_length > limit) {
        fprintf(stderr,"WARNING: mini k-mer bits cannot be larger than double the full k-mer length %u!  Using mini k-mer bits of %u instead.\n",
            g__FULLKMER_LENGTH, limit);
        minikmer_length = limit;
    }
#else
    unsigned limit = g__FULLKMER_LENGTH;
    if (minikmer_length > limit) {
        fprintf(stderr,"WARNING: mini k-mers cannot be larger than the full k-mer length %u!  Using mini k-mer length of %u instead.\n",
            g__FULLKMER_LENGTH, limit);
        minikmer_length = limit;
    }
#endif

#ifdef PART_MIDDLE_MINIKMER
    if (minikmer_length % 2 == 0) {
        fprintf(stderr,"WARNING: mini k-mers can't be of even length, such as %u.  Using mini k-mer length of %u instead.\n",
            minikmer_length, minikmer_length-1);
        -- minikmer_length;
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

  g__PREGRAPH_FILENAME = static_cast<char *>( calloc((PATH_MAX+1), sizeof(char)) );
  strcpy(g__PREGRAPH_FILENAME, g__WORK_BASE_DIRECTORY);
  strcat(g__PREGRAPH_FILENAME, "/PreGraph");

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
      } else if (strcmp(argv[argIndex], "-fasta") == 0) {
		  current_entry.filetype = FASTA;
      } else if (strcmp(argv[argIndex], "-gerald") == 0) {
		  current_entry.filetype = GERALD;
      } else if (strcmp(argv[argIndex], "-eland") == 0) {
		  current_entry.filetype = ELAND;
      } else if (strcmp(argv[argIndex], "-fastq.gz") == 0) {
		  current_entry.filetype = FASTQ_GZ;
      } else if (strcmp(argv[argIndex], "-fasta.gz") == 0) {
		  current_entry.filetype = FASTA_GZ;
      } else if (strcmp(argv[argIndex], "-maq.gz") == 0) {
		  current_entry.filetype = MAQ_GZ;
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
			g__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            g__PARTITION_COUNT = getNumPartitions(argv[++argIndex]);
            g__MINIKMER_LENGTH = getMiniKmerLength(argv[++argIndex]);
      } else if (strcmp(argv[argIndex], "-noemit") == 0) {
            g__PARTITION_NOEMIT = true;
      } else if (strcmp(argv[argIndex], "-minfootprint") == 0) {
            g__MINIMIZE_FOOTPRINT = true;
      } else if (strcmp(argv[argIndex], "-loom") == 0) {
            g__LOOMING = true;
			g__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = LOOM;
      } else if (strcmp(argv[argIndex], "-quilt") == 0) {
            g__QUILTING = true;
			g__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = QUILT;
      } else if (strcmp(argv[argIndex], "-split") == 0) {
            g__SPLIT = true;
			g__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = QUILT;
            g__PARTITION_COUNT = getNumPartitions(argv[++argIndex]);
      } else if (strcmp(argv[argIndex], "-isplit") == 0) {
            g__ISPLIT = true;
			g__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = QUILT;
            g__PARTITION_COUNT = getNumPartitions(argv[++argIndex]);
      } else if (strcmp(argv[argIndex], "-iesplit") == 0) {
            g__ISPLIT = true;
            g__IESPLIT = true;
			g__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = QUILT;
            g__PARTITION_COUNT = getNumPartitions(argv[++argIndex]);
      } else if (strcmp(argv[argIndex], "-ciesplit") == 0) {
            g__CIESPLIT = true;
			g__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = QUILT;
            g__PARTITION_COUNT = getNumPartitions(argv[++argIndex]);
      } else if (strcmp(argv[argIndex], "-bucket") == 0) {
            current_entry.filetype = BUCKET;
      } else if (strcmp(argv[argIndex], "-combine") == 0) {
            g__COMBINING = true;
			g__SINGLEPASS = false;
            g__PSEUDO_NODES_PRESENT = true;
            p__CLOBBER_WORK_DIRECTORY = false;
            current_entry.filetype = QUILT;
      } else if (strcmp(argv[argIndex], "-flow") == 0) {
            g__FLOW = true;
			g__SINGLEPASS = false;
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
            g__PGPART = true;
      } else if (strcmp(argv[argIndex], "-pgdist") == 0) {
            g__PGDIST = true;
			g__SINGLEPASS = false;
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
          fprintf(stderr, "Unknown option: '%s'\n", argv[argIndex]);
          exit(EXIT_FAILURE);
      }

      continue;
    }

    // input file(s): test if its a directory or a filename
    DIR *dir = opendir(argv[argIndex]);
    if (dir != NULL) { // TODO: support subdirectories???
      // if a directory, add all regular files to the input list
      if( g__LOOMING ) {
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
  if( g__PARTITIONING + g__LOOMING + g__SPLIT + g__ISPLIT + g__CIESPLIT + g__COMBINING + g__QUILTING > 1 ) {
      fprintf(stderr, "Invalid combination of '-part' '-loom' '-combine' '-split,etc' and/or '-quilt'.  Exiting...\n");
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

void runSplit(file_object_vector *file_objects)
{
    initializePartitionIndexFromInputFilename(file_objects);

    // next, load the quilt file(s)
    for(file_object_vector::iterator itr = file_objects->begin() ; itr != file_objects->end() ; ++itr) {
        switch (itr->filetype) {
            case BUCKET:
                sg_load_bucket(g__SG_HASHTABLE, itr->filename);
                break;
            case QUILT:
                sg_load_quilt(g__SG_HASHTABLE, itr->filename);
                break;
            default:
                fprintf(stderr, "ERROR: Cannot use %s file %s as quilt input.  Exiting...\n",
                        FILE_TYPES[itr->filetype], itr->filename);
                exit(EXIT_FAILURE);
        }
    }

    uintptr_t peak_nodes = g__SG_HASHTABLE->node_count;
    printf("%"PRIuPTR" peak sequence nodes (bulk split).\n", peak_nodes);

    if( g__FULL_STATISTICS ) {
        sg_stat_components(g__SG_HASHTABLE, stdout);
    }

    if (g__DOTGRAPH) {
        char dot_filename[PATH_MAX+1];
        sprintf(dot_filename, "%s/Split-%u.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX);
        emit_graphviz(g__SG_HASHTABLE, dot_filename);
    }


    if (file_objects->size() > 1) {
#ifdef VERIFY
        g__SG_HASHTABLE->verify();
#endif
        sg_remove_tips(g__SG_HASHTABLE);
#ifdef VERIFY
        g__SG_HASHTABLE->verify();
#endif
        sg_concatenate(g__SG_HASHTABLE);
#ifdef VERIFY
        g__SG_HASHTABLE->verify();
#endif
    }

    if (g__DOTGRAPH) {
        char dot_filename[PATH_MAX+1];
        sprintf(dot_filename, "%s/Split-%u-simplify.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX);
        emit_graphviz(g__SG_HASHTABLE, dot_filename);
    }

    if (g__SLICING) {
        slice2_graph(g__SG_HASHTABLE, g__PARTITION_INDEX);
        printf("%"PRIuPTR" nodes sliced out to final.\n", g__SLICE2_FINAL_NODE_COUNT);
        printf("%"PRIuPTR" nodes sliced out to others.\n", g__SLICE2_NODE_COUNT);

        if (g__DOTGRAPH) {
            char dot_filename[PATH_MAX+1];
            sprintf(dot_filename, "%s/Split-%u-sliced.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX);
            emit_graphviz(g__SG_HASHTABLE, dot_filename);
        }
    }

    // last, form components and emit to buckets
    SplitBuckets *buckets = new SplitBuckets(false);
    buckets->split(g__SG_HASHTABLE);
    buckets->printStatistics();
    delete buckets;
}

void runISplit(file_object_vector *file_objects)
{
    initializePartitionIndexFromInputFilename(file_objects);

    // fully load the first file
    file_object_vector::iterator itr = file_objects->begin();
    assert(itr != file_objects->end());
    switch (itr->filetype) {
        case BUCKET:
            sg_load_bucket(g__SG_HASHTABLE, itr->filename);
            break;
        case QUILT:
            sg_load_quilt(g__SG_HASHTABLE, itr->filename);
            break;
        default:
            fprintf(stderr, "ERROR: Cannot use %s file %s as quilt input.  Exiting...\n",
                    FILE_TYPES[itr->filetype], itr->filename);
            exit(EXIT_FAILURE);
    }
    ++itr;

    // we just loaded the 'self' bucket, no simplification yet to do

    if( g__FULL_STATISTICS ) {
        sg_stat_components(g__SG_HASHTABLE, stdout);
    }

    uintptr_t self_nodes = g__SG_HASHTABLE->node_count;
    uintptr_t ceiling = static_cast<uintptr_t>(1.2 * self_nodes); // TODO FIXME: want 1.2x of initial nodes, NOT self nodes -- maybe global value set during looming?

    uintptr_t min_limit = static_cast<uintptr_t>(0.05 * self_nodes); // always load at least 5% of initial graph size
    min_limit = (min_limit == 0 ? 1 : min_limit);

    uintptr_t peak_nodes = g__SG_HASHTABLE->node_count;

    uintptr_t total_bucket_nodes = 0;

    SplitBuckets *buckets = new SplitBuckets(true);

    // incrementally load the other files
    for( ; itr != file_objects->end() ; ++itr) {
        unsigned round = 1;
        assert( itr->filetype == BUCKET );
        int filedes = open(itr->filename, O_RDONLY);
        FILE* file = fdopen(filedes, "r"); // FIXME: performance -- setbuffer() a larger buffer?
        while(!feof(file)) {
            printf("=== Incrementally loading INBOX bucket: Round %u ===\n", round);

            uintptr_t limit;
            if (g__IESPLIT && g__SG_HASHTABLE->node_count < ceiling) {
                limit = ceiling - g__SG_HASHTABLE->node_count;
            } else {
                if (g__SG_HASHTABLE->node_count < (0.95 * peak_nodes)) {
                    limit = peak_nodes - g__SG_HASHTABLE->node_count;
                } else {
                    limit = static_cast<uintptr_t>(1.1 * peak_nodes) - g__SG_HASHTABLE->node_count;
                }
            }
            limit = max(limit, min_limit);

            printf("%"PRIuPTR" nodes in working graph before loading from stream bucket.\n", g__SG_HASHTABLE->node_count);
            printf("%"PRIuPTR" load limit.\n", limit);
            total_bucket_nodes += sg_load_stream_bucket(g__SG_HASHTABLE, limit, file);
            peak_nodes = max(peak_nodes, g__SG_HASHTABLE->node_count);

            printf("%"PRIuPTR" nodes in working graph after loading from stream bucket.\n", g__SG_HASHTABLE->node_count);

            if( g__FULL_STATISTICS ) {
                sg_stat_components(g__SG_HASHTABLE, stdout);
            }

            if (g__DOTGRAPH) {
                char dot_filename[PATH_MAX+1];
                sprintf(dot_filename, "%s/SplitBucket-%u-%u.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
                emit_graphviz(g__SG_HASHTABLE, dot_filename);
            }

#ifdef VERIFY
            g__SG_HASHTABLE->verify();
#endif
            sg_remove_tips(g__SG_HASHTABLE);
#ifdef VERIFY
            g__SG_HASHTABLE->verify();
#endif
            sg_concatenate(g__SG_HASHTABLE);
#ifdef VERIFY
            g__SG_HASHTABLE->verify();
#endif

            if (g__DOTGRAPH) {
                char dot_filename[PATH_MAX+1];
                sprintf(dot_filename, "%s/SplitBucket-%u-%u-simplify.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
                emit_graphviz(g__SG_HASHTABLE, dot_filename);
            }

            // emit sub-components that are no longer relevant to  
            if (g__SLICING) {
                slice2_graph(g__SG_HASHTABLE, g__PARTITION_INDEX);

                if (g__DOTGRAPH) {
                    char dot_filename[PATH_MAX+1];
                    sprintf(dot_filename, "%s/SplitBucket-%u-%u-sliced.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
                    emit_graphviz(g__SG_HASHTABLE, dot_filename);
                }
            }

            // emit components that are no longer relevant to this partition
            if (g__IESPLIT && !feof(file)) {
                buckets->split(g__SG_HASHTABLE);
            } else {
                g__SG_HASHTABLE->resetFlags(); // XXX: redundant?
            }

            ++round;
        }
        fclose(file);
    }

    printf("%"PRIuPTR" peak sequence nodes (incremental split).\n", peak_nodes);
    printf("%"PRIuPTR" total sequence nodes loaded from stream bucket(s).\n", total_bucket_nodes);

    if (g__SLICING) {
        printf("%"PRIuPTR" nodes sliced out to final.\n", g__SLICE2_FINAL_NODE_COUNT);
        printf("%"PRIuPTR" nodes sliced out to others.\n", g__SLICE2_NODE_COUNT);
    }

    if( g__FULL_STATISTICS ) {
        sg_stat_components(g__SG_HASHTABLE, stdout);
    }

    // last, form components and emit to buckets
    buckets->split(g__SG_HASHTABLE);
    buckets->printStatistics();
    delete buckets;
}

void sg_serial_import_related_components(SeqGraph *working_graph, SeqGraph *resident_graph); // forward decl

void runCIESplit(file_object_vector *file_objects)
{
    initializePartitionIndexFromInputFilename(file_objects);

    // fully load the first file
    file_object_vector::iterator itr = file_objects->begin();
    assert(itr != file_objects->end());
    switch (itr->filetype) {
        case BUCKET:
            sg_load_bucket(g__SG_HASHTABLE, itr->filename);
            break;
        case QUILT:
            sg_load_quilt(g__SG_HASHTABLE, itr->filename);
            break;
        default:
            fprintf(stderr, "ERROR: Cannot use %s file %s as quilt input.  Exiting...\n",
                    FILE_TYPES[itr->filetype], itr->filename);
            exit(EXIT_FAILURE);
    }
    ++itr;

    // we just loaded the 'self' bucket, no simplification yet to do

    if( g__FULL_STATISTICS ) {
        sg_stat_components(g__SG_HASHTABLE, stdout);
    }

    uintptr_t peak_nodes = g__SG_HASHTABLE->node_count;

    uintptr_t total_bucket_nodes = 0;

    SeqGraph *g__WORK_GRAPH = new SeqGraph(PRIVATE_BUCKETS);

    SplitBuckets *buckets = new SplitBuckets(true);

    // incrementally load components from the other files
    for( ; itr != file_objects->end() ; ++itr) {
        unsigned round = 1;
        assert( itr->filetype == BUCKET );
        int filedes = open(itr->filename, O_RDONLY);
        FILE* file = fdopen(filedes, "r"); // FIXME: performance -- setbuffer() a larger buffer?
        while(!feof(file)) {

            total_bucket_nodes += sg_load_stream_component(g__WORK_GRAPH, file);
            peak_nodes = max(peak_nodes, g__SG_HASHTABLE->node_count + g__WORK_GRAPH->node_count);

            sg_serial_import_related_components(g__WORK_GRAPH, g__SG_HASHTABLE);

            if (g__DOTGRAPH) {
                char dot_filename[PATH_MAX+1];
                sprintf(dot_filename, "%s/SplitBucket-%u-%u-resident.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
                emit_graphviz(g__SG_HASHTABLE, dot_filename);

                char dot_filename2[PATH_MAX+1];
                sprintf(dot_filename2, "%s/SplitBucket-%u-%u-work.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
                emit_graphviz(g__WORK_GRAPH, dot_filename2);
            }

            /*if( g__FULL_STATISTICS ) {
              sg_stat_components(g__WORK_GRAPH, stdout);
              }*/

#ifdef VERIFY
            g__WORK_GRAPH->verify(true);
#endif
            sg_remove_tips(g__WORK_GRAPH, true);
#ifdef VERIFY
            g__WORK_GRAPH->verify(true);
#endif
            sg_concatenate(g__WORK_GRAPH, true);
#ifdef VERIFY
            g__WORK_GRAPH->verify(true);
#endif

            /*if (g__DOTGRAPH) {
              char dot_filename[PATH_MAX+1];
              sprintf(dot_filename, "%s/SplitBucket-%u-%u-simplify.dot", g__WORK_BASE_DIRECTORY, currentPartitionIndex, round);
              emit_graphviz(g__WORK_GRAPH, dot_filename);
              }*/

            // emit sub-components that are no longer relevant to  
            if (g__SLICING) {
                slice2_graph(g__WORK_GRAPH, g__PARTITION_INDEX);

                if (g__DOTGRAPH) {
                    char dot_filename[PATH_MAX+1];
                    sprintf(dot_filename, "%s/SplitBucket-%u-%u-worksliced.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
                    emit_graphviz(g__WORK_GRAPH, dot_filename);
                }
            }

            // emit components that are no longer relevant to the working graph
            if (!feof(file)) {
                buckets->split(g__WORK_GRAPH);
            }

            if (g__DOTGRAPH) {
                char dot_filename[PATH_MAX+1];
                sprintf(dot_filename, "%s/SplitBucket-%u-%u-resident-postemit.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
                emit_graphviz(g__SG_HASHTABLE, dot_filename);

                char dot_filename2[PATH_MAX+1];
                sprintf(dot_filename2, "%s/SplitBucket-%u-%u-work-postemit.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
                emit_graphviz(g__WORK_GRAPH, dot_filename2);
            }

            // move remaining nodes from working graph into resident graph
            g__WORK_GRAPH->bulkMoveAllNodes(g__SG_HASHTABLE);
            assert( g__WORK_GRAPH->node_count == 0 );

            //g__SG_HASHTABLE->verify(true);

            ++ round;
        }
        fclose(file);
    }

    printf("%"PRIuPTR" peak sequence nodes (incremental split).\n", peak_nodes);
    printf("%"PRIuPTR" total sequence nodes loaded from stream bucket(s).\n", total_bucket_nodes);

    if (g__SLICING) {
        printf("%"PRIuPTR" nodes sliced out to final.\n", g__SLICE2_FINAL_NODE_COUNT);
        printf("%"PRIuPTR" nodes sliced out to others.\n", g__SLICE2_NODE_COUNT);
    }

    /*if( g__FULL_STATISTICS ) {
      sg_stat_components(g__SG_HASHTABLE, stdout);
      }*/

    // last, form components and emit to buckets
    buckets->split(g__SG_HASHTABLE);
    buckets->printStatistics();
    delete buckets;
}

void runCombine(file_object_vector *file_objects)
{
    initializePartitionIndexFromInputFilename(file_objects);

    char filename[PATH_MAX+1];
    strncpy(filename, file_objects->front().filename, PATH_MAX);
    filename[PATH_MAX] = '\0';

    char *end = strstr(filename, ".quilt");
    strcpy(end, ".combined.quilt");
    quilt_files(g__SG_HASHTABLE, *file_objects);
    sg_dump_quilt(g__SG_HASHTABLE, filename);
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

    sg_dump_pregraph(g__SG_HASHTABLE, g__PREGRAPH_FILENAME);

    if (g__PGPART) {
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

    if (g__PGPART) {
        char quilt_filename[PATH_MAX+1];
        sprintf(quilt_filename, "%s/Direct.quilt", g__WORK_QUILT_DIRECTORY);
        sg_dump_quilt_from_kmergraph(g__KG_HASHTABLE, quilt_filename);
    } else {
        sg_dump_pregraph_from_kmergraph(g__KG_HASHTABLE, g__PREGRAPH_FILENAME);
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
    if( g__SINGLEPASS || g__LOOMING ) {
        g__KG_HASHTABLE = new KmerGraph(buckets);
    } else if (g__FLOW) {
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
    } else if (g__LOOMING) {
        runLooming(file_objects);
    } else if (g__SPLIT) {
        runSplit(file_objects);
    } else if (g__ISPLIT) {
        runISplit(file_objects);
    } else if (g__CIESPLIT) {
        runCIESplit(file_objects);
    } else if (g__COMBINING) {
        runCombine(file_objects);
    } else if (g__QUILTING) {
        runQuilt(file_objects);
    } else if (g__FLOW) {
        runFlow(file_objects, g__KG_HASHTABLE, g__SG_HASHTABLE);
    } else if (g__PGDIST) {
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

