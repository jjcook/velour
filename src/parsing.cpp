//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// parsing.cpp
//
// parses reads files
//

#include "types.h"

#include "minikmer.h"

static unsigned g_tolerable_unknowns = 0; // TODO: jjcook configurable parameter

// note: keep consistent with the enum in types.h
const char *FILE_TYPES[] = {"Unknown", "FastQ", "FastA", "Solexa", "Eland", "FastQ_GZ", "FastA_GZ", "MAQ_GZ", "PreGraph", "Loom", "Quilt", "Bucket", "NoFormat"};

static uint64_t readFastaFile(char *file, off_t length, Category category, KmerGraph *hashtable, int kmer_length);
static uint64_t readSolexaFile(char *file, off_t length, Category category, KmerGraph *hashtable, int kmer_length);
static uint64_t readLoomFile(char *file, off_t length, Category category, KmerGraph *hashtable, int kmer_length);

//********************************************************************************
//*******************    Base Translation / Manipulations    *********************
//********************************************************************************

// maps characters (e.g., 'A' and 'a') to a base (ADENINE)
Nucleotide BASE_MAP[256];
char CHAR_BASE_MAP[4];
char CONNECTION_COUNT_MAP[16];
char CONNECTION_MAP[16];

void initializeBaseMap(void) {
#ifdef VELVET_EMULATION
  memset(BASE_MAP, ADENINE, 256);  // set all unknown characters to 'A'
#else
  memset(BASE_MAP, -1, 256);  // make all unknown characters invalid
#endif
  BASE_MAP[int('A')] = BASE_MAP[int('a')] = ADENINE;
  BASE_MAP[int('C')] = BASE_MAP[int('c')] = CYTOSINE;
  BASE_MAP[int('G')] = BASE_MAP[int('g')] = GUANINE;
  BASE_MAP[int('T')] = BASE_MAP[int('t')] = THYMINE;

  CHAR_BASE_MAP[ADENINE] = 'A';
  CHAR_BASE_MAP[CYTOSINE] = 'C';
  CHAR_BASE_MAP[GUANINE] = 'G';
  CHAR_BASE_MAP[THYMINE] = 'T';
  
  memset(CONNECTION_MAP, -1, 16);
  CONNECTION_MAP[(1<<0)] = 0;  CONNECTION_MAP[(1<<1)] = 1;
  CONNECTION_MAP[(1<<2)] = 2;  CONNECTION_MAP[(1<<3)] = 3;
  for (int i = 0 ; i < 16 ; ++ i) {
	 CONNECTION_COUNT_MAP[i] = ((i&1)?1:0) + ((i&2)?1:0) + ((i&4)?1:0) + ((i&8)?1:0);
  }	 
}


// mmaps the file (read-only) specified by the file_object pointed to
// by the iterator.  As a side effect it sets the length of the file
// in the iterator.
static char *
//open_file(file_object_vector::iterator itr) {
open_file(file_object_t * itr) {
  	 struct stat file_stat;
	 int fd = open(itr->filename, O_RDONLY);
	 if (fd == -1) {
		printf("failed to open file: %s\n", itr->filename); 
		exit(1);
	 }
	 if (fstat(fd, &file_stat) != 0) {
		printf("failed to stat file: %s\n", itr->filename); 
		exit(1);
	 }
	 off_t length = file_stat.st_size;
     if (length > 0) {
        void* file = mmap(0, length, PROT_READ, MAP_PRIVATE, fd, 0);
        if (file == (void *)-1) {
            printf("failed to mmap file: %s\n", itr->filename);
            exit(1);
        }
        itr->length = length;
        return (char *)file;
     } else {
         return NULL;
     }
}

#ifdef VELOUR_TBB
void parallel_load_loom_files(file_object_vector &work_to_do);
#endif

// Serial version                                                                             
void 
load_loom_files(file_object_vector &work_to_do)
{
#ifdef VELOUR_TBB
    parallel_load_loom_files(work_to_do);
#else
    for(file_object_vector::iterator itr = work_to_do.begin() ; itr != work_to_do.end() ; ++itr) {
        // open the file using mmap
        char *file = open_file(&(*itr));
        uint64_t num_subsequences = 0;

        switch (itr->filetype) {
        case LOOM:
            num_subsequences =   readLoomFile(file, itr->length, itr->cat, g__KG_HASHTABLE, g__FULLKMER_LENGTH);
            break;
        default:
            fprintf(stderr, "ERROR: Cannot use %s file %s as loom input.  Exiting...\n",
                FILE_TYPES[itr->filetype], itr->filename);
            exit(EXIT_FAILURE);
        }

        printf("Reading %s file %s\t%" PRIu64 " subsequences found.\n",
            FILE_TYPES[itr->filetype], itr->filename, num_subsequences);
    
        // close the file using munmap
        if (file != NULL) {
            munmap((void *)file, itr->length);
        }
    }
#endif // VELOUR_TBB
}

#ifdef VELOUR_TBB
typedef struct mmap_chunk {
    file_object_t *file_info;
    char * mmap_chunk;
    size_t mmap_length;
} file_chunk_t;

class LoomLoadFilter : public tbb::filter
{
    public:
        static const size_t n_buffer = 64;  // XXX: controls the concurrency level

    private:
        file_object_vector& files_to_load_;
        int file_index_;
        off_t file_offset_;
        int filedes_;

        volatile uintptr_t touch_total_;  // fool the optimizer

        file_chunk_t buffer_[n_buffer]; // TODO OPT: align each to cache line so no false sharing
        size_t next_buffer_;

        bool advance_file()
        {
            if (filedes_ != -1) { // close current file
                if (close(filedes_)) {
                    fprintf(stderr, "ERROR: failed to close file: %s\n", files_to_load_[file_index_].filename);
                    perror("REASON: ");
                    exit(EXIT_FAILURE);
                }
                filedes_ = -1;
            }

            file_offset_ = 0;   // reset offset

            ++ file_index_;
            if (file_index_ >= static_cast<int>(files_to_load_.size())) { // no more files
                return false;
            }

            if (files_to_load_[file_index_].filetype != LOOM) {
                fprintf(stderr, "ERROR: Cannot use %s file %s as loom input.  Exiting...\n",
                    FILE_TYPES[files_to_load_[file_index_].filetype], files_to_load_[file_index_].filename);
                exit(EXIT_FAILURE);
            }

            filedes_ = open(files_to_load_[file_index_].filename, O_RDONLY);
            if (filedes_ == -1) {
                fprintf(stderr, "ERROR: failed to open file: %s\n", files_to_load_[file_index_].filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }

            struct stat file_stat;
            if (fstat(filedes_, &file_stat) != 0) {
                fprintf(stderr,"ERROR: failed to stat file: %s\n", files_to_load_[file_index_].filename);
                exit(EXIT_FAILURE);
            }

            files_to_load_[file_index_].length = file_stat.st_size;
            return true;
        }

    public:
        LoomLoadFilter(file_object_vector& work_to_do) : tbb::filter(serial_in_order),
            files_to_load_(work_to_do), file_index_(-1), file_offset_(0), filedes_(-1),
            touch_total_(0), next_buffer_(0)
        {
            advance_file();
        }

        void* operator()(void* __dummy)
        {
            if (file_offset_ >= files_to_load_[file_index_].length) { // advance to next file
                if (!advance_file()) {
                    return NULL;
                }
            }

            off_t length_to_grab = min( 1024 * 1024, (files_to_load_[file_index_].length - file_offset_) ); // XXX: configurable parameter -- don't be too small, else mmap subtly breaks!?
            assert( file_offset_ + length_to_grab <= files_to_load_[file_index_].length );

            file_chunk_t *chunk = &buffer_[next_buffer_];
            next_buffer_ = (next_buffer_ + 1) % n_buffer;

            chunk->file_info = &files_to_load_[file_index_];
            chunk->mmap_chunk = static_cast<char *>( mmap(0, length_to_grab, PROT_READ, MAP_PRIVATE, filedes_, file_offset_) );
            if (chunk->mmap_chunk == reinterpret_cast<char *>(-1)) {
                fprintf(stderr, "ERROR: failed to mmap file: %s\n", files_to_load_[file_index_].filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }
            chunk->mmap_length = length_to_grab;

            file_offset_ += length_to_grab;

            // touch the pages
            for (size_t touch = 0 ; touch < chunk->mmap_length ; touch += 4096) {
                touch_total_ += chunk->mmap_chunk[touch];
            }

            return chunk;
        }
};

class LoomProcessFilter : public tbb::filter
{
    public:
        LoomProcessFilter(KmerGraph *kgraph, unsigned kmer_length) : tbb::filter(parallel),
            kgraph_(kgraph), kmer_length_(kmer_length) {}

        void * operator()(void *item)
        {
            file_chunk_t& chunk = * static_cast<file_chunk_t *>(item);

            uint64_t num_subsequences = readLoomFile(chunk.mmap_chunk, chunk.mmap_length, chunk.file_info->cat, kgraph_, kmer_length_);

            int retval = munmap(chunk.mmap_chunk, chunk.mmap_length);
            if (retval == -1) {
                fprintf(stderr, "ERROR: failed to munmap file: %s\n", chunk.file_info->filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }
            return NULL;
        }
    private:
        KmerGraph *kgraph_;
        unsigned kmer_length_;
};

void parallel_load_loom_files(file_object_vector &work_to_do)
{
    tbb::tick_count time0, time1;
    tbb::pipeline pipeline;

    printf("  parallel loom: %zu tokens\n", LoomLoadFilter::n_buffer); fflush(stdout);

    time0 = tbb::tick_count::now();

    LoomLoadFilter lf(work_to_do);
    pipeline.add_filter( lf );

    LoomProcessFilter pf(g__KG_HASHTABLE, g__FULLKMER_LENGTH);
    pipeline.add_filter( pf );

    pipeline.run( LoomLoadFilter::n_buffer );

    time1 = tbb::tick_count::now();

    tbb::tick_count::interval_t loom_time = time1 - time0;

    printf("  loom time: %lfs\n", loom_time.seconds()); fflush(stdout);
}
#endif // VELOUR_TBB

// Serial version                                                                             
void 
process_files(file_object_vector &work_to_do)
{
    unsigned file_count = 0;
    uint64_t total_sequences = 0;
    for(file_object_vector::iterator itr = work_to_do.begin() ; itr != work_to_do.end() ; ++itr) {
        // open the file using mmap
        char *file = open_file(&(*itr));
        uint64_t file_sequences = 0;

        ++file_count;
        
        printf("Processing (%04d/%04d) %s file %s:\n", file_count, work_to_do.size(),
                FILE_TYPES[itr->filetype], itr->filename);
        fflush(stdout);

        switch (itr->filetype) {
        case FASTA:
            file_sequences =  readFastaFile(file, itr->length, itr->cat, g__KG_HASHTABLE, g__FULLKMER_LENGTH);
            break;
        case GERALD:
            file_sequences = readSolexaFile(file, itr->length, itr->cat, g__KG_HASHTABLE, g__FULLKMER_LENGTH);
            break;
        default:
            fprintf(stderr, "ERROR: Cannot use %s file %s as input.  Exiting...\n",
                FILE_TYPES[itr->filetype], itr->filename);
            exit(EXIT_FAILURE);
        }

        printf("Finished %s file %s\t%" PRIu64 " sequences found.\n",
            FILE_TYPES[itr->filetype], itr->filename, file_sequences);

        total_sequences += file_sequences;
    
        // close the file using munmap
        munmap((void *)file, itr->length);
    }

    printf("%" PRIu64 " total sequences found.\n", total_sequences);

    // communicate number of total sequences back to environment
    char filename[PATH_MAX+1];
    sprintf(filename, "%s/common.sequences", g__WORK_BASE_DIRECTORY);
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        fprintf(stderr, "ERROR: fopen()\n");
        assert(false);
        exit(EXIT_FAILURE);
    }
    fprintf(f, "%"PRIu64"", total_sequences);
    fclose(f);

    g__READ_COUNT = total_sequences;
}

static inline void
processSequence(char *seq, KmerGraph *hashtable, int kmer_length)
{
  if( g__PHASE == PHASE_BUILD_KMER_GRAPH )
    convertSequenceToKmersToPrenodes(seq,hashtable,kmer_length);
  else if( g__PHASE == PHASE_DISTRIBUTE_SEQUENCES )
    distributeSequence(seq,kmer_length);
#ifdef PART_TRAINS_ON_INPUT
  else if( g__PHASE == PHASE_OBSERVE_SEQUENCES )
    observeSequence(seq,kmer_length);
#endif // PART_TRAINS_ON_INPUT
  else {
    report_unrecoverable_error();
    exit(EXIT_FAILURE);
  }
}

static void
dealWithUnknownsAndConvert(char *seq, KmerGraph *hashtable, int kmer_length) {
  int i;
  char c;
  int outer_start = 0;
  int inner_start = 0;
  int inner_end = -1;
  unsigned num_unknowns = 0;

  for (i = 0 ; (c = seq[i]) != 0 ; i ++) {
	 if (c == 'N' || c == '.') { // unknown base -- FIXME recognize other characters
		if ((i - outer_start) < kmer_length) {
		  inner_start = outer_start = i+1;      
		  num_unknowns = 0;
		} else {
		  inner_start = i+1;
		  num_unknowns ++;
		  if (num_unknowns > g_tolerable_unknowns) {
			 if (outer_start != inner_end) {
				seq[inner_end] = 0;
				processSequence(&seq[outer_start], hashtable, kmer_length);
			 }
			 outer_start = i+1;
			 num_unknowns = 0;
		  }
		}
	 } else {        // a valid base
		if ((i - inner_start) >= kmer_length) {  // found a stretch at least "k" long!
		  inner_end = i+1;
		}
	 }
  }

  if (outer_start != i) {
	 processSequence(&seq[outer_start], hashtable, kmer_length);
  }
}

// Imports sequences from a solexa file directly into pre_nodes
static uint64_t
readSolexaFile(char *file, off_t length, Category category, KmerGraph *hashtable, int kmer_length) {
  uint64_t readIndex = 0;
  char seq[MAX_READ_LENGTH];
  char c;
  off_t i = 0;

  while ((i + kmer_length) < length) {  // bail if a partial line remains
	 for (int n = 0 ; n < 4 ; n ++) {
		while ((i < length) && (file[i++] != '\t')) { }// parse n'th number (are 4 before the string) 
	 }
	 int j = 0;
	 while (i < length) {                 // copy out the string, will truncate to MAX_READ_LENGTH length
		if ((c = file[i++]) == '\n')		     break;
		seq[j] = c;
		if (j < MAX_READ_LENGTH) {
		  j ++;
		}
	 }
	 if (j > kmer_length) {
		seq[j] = 0;
		dealWithUnknownsAndConvert(seq, hashtable, kmer_length);
		readIndex++;
	 }

	 // while (fgets(line, MAX_READ_LENGTH, file) != NULL) {
    //   // sscanf(line, "%*i\t%*i\t%*i\t%*i\t%*c%[^\n]", seq);
    //   sscanf(line, "%*i\t%*i\t%*i\t%*i\t%[^\n]", seq);
    //   // convertSequenceToKmersToPrenodes(seq, hashtable, kmer_length);
	 // 	dealWithUnknownsAndConvert(seq, hashtable, kmer_length);
    //   // newTightStringFromStringNoFree(seq);
    //   readIndex++;
    // }
  }
  return readIndex;
}

// Imports sequences from a FASTA file directly into pre_nodes
static uint64_t
readFastaFile(char *file, off_t length, Category category, KmerGraph *hashtable, int kmer_length)
{
    uint64_t read_count = 0;
    static uint64_t hack_read_count = 0; // TODO: make global, for use in Velvet sequences file

    char header[MAX_READ_HEADER_LENGTH+1];
    unsigned header_write_index = 0;

    char seq[MAX_READ_LENGTH+1];
    unsigned seq_write_index = 0;

    char c;
    off_t i = 0;
    while (i < length) {
        if ((c = file[i]) == '>') { // is a header line
            if ((read_count % 1000000) == 0) {
                printf("Processing sequence %9"PRIu64"...\n", read_count); fflush(stdout);
            }
            if (header_write_index != 0) { // emit to Velvet 'Sequences' file: header + sequence (60 bases per line)
                assert( g__VELVET_SEQUENCES_FILE != NULL );
                fprintf(g__VELVET_SEQUENCES_FILE, "%s\t%"PRIu64"\t%u\n", header, hack_read_count+1, category); // emit header + metadata
                // NOTE: seq_write_index points at the terminating null character
                unsigned emit_start_index = 0;
                while (emit_start_index < seq_write_index) {
                    unsigned line_bases = min(60, (seq_write_index - emit_start_index));
                    char *line_start = &seq[emit_start_index];
                    char *line_end   = line_start + line_bases;
                    char save_base = *line_end;
                    *line_end = '\0';
                    fprintf(g__VELVET_SEQUENCES_FILE, "%s\n", line_start);
                    *line_end = save_base; // restore
                    emit_start_index += line_bases;
                }
                header_write_index = 0;
            }
            if (seq_write_index != 0) {  // the header marked the end of a chunk of sequence
                dealWithUnknownsAndConvert(seq, hashtable, kmer_length);
                ++ read_count;
                ++ hack_read_count;
                seq_write_index = 0;
            }
            if (g__VELVET_SEQUENCES_FILE != NULL) { // buffer header line
                while ((i < length) && ((c = file[i++]) != '\n')) {
                    header[header_write_index ++] = c;
#ifdef VERIFY
                    assert( header_write_index < MAX_READ_HEADER_LENGTH );
#endif
                }
                header[header_write_index] = '\0';
            } else { // ignore header
                while ((i < length) && (file[i++] != '\n')) { } // scan to end of line
            }
        } else {
            while ((i < length) && (c = file[i++]) != '\n') {
                seq[seq_write_index ++] = c;
#ifdef VERIFY
                assert( seq_write_index < MAX_READ_LENGTH );
#endif
            }
            seq[seq_write_index] = '\0';
#ifdef VELVET_EMULATION
            for (unsigned j=0; j < seq_write_index ; ++j) { seq[j] = CHAR_BASE_MAP[BASE_MAP[seq[j]]]; } // XXX: velvetify sequence
#endif
        }
    }

    if (header_write_index != 0) { // emit to Velvet 'Sequences' file
        assert( g__VELVET_SEQUENCES_FILE != NULL );
        fprintf(g__VELVET_SEQUENCES_FILE, "%s\t%"PRIu64"\t%u\n", header, hack_read_count+1, category); // emit header + metadata
        // NOTE: seq_write_index points at the terminating null character
        unsigned emit_start_index = 0;
        while (emit_start_index < seq_write_index) {
            unsigned line_bases = min(60, (seq_write_index - emit_start_index));
            char *line_start = &seq[emit_start_index];
            char *line_end   = line_start + line_bases;
            char save_base = *line_end;
            *line_end = '\0';
            fprintf(g__VELVET_SEQUENCES_FILE, "%s\n", line_start);
            *line_end = save_base; // restore
            emit_start_index += line_bases;
        }
        header_write_index = 0;
    }
    if (seq_write_index != 0) { // not necessarily a comment after last bit of sequence
        dealWithUnknownsAndConvert(seq, hashtable, kmer_length);
        ++ read_count;
        ++ hack_read_count;
    }
  
    return read_count;
}

// Imports sequences from a LOOM file directly into pre_nodes
static uint64_t
readLoomFile(char *file, off_t length, Category category, KmerGraph *hashtable, int kmer_length) {
  uint64_t readIndex = 0;
  off_t i = 0;

  while (i < length) {
	char hasPrefixSuffix = file[i];

    if (hasPrefixSuffix == -1) {    // XXX: aligned write buffer hack
        i += 16384 - (i % 16384);   // TODO FIXME: constants!!!
        assert( i % 16384 == 0 );
        continue;
    }

	++i;

    Sequence_StackAllocated memory;
    Sequence *stack_seq = new (&memory) Sequence(Sequence::MAX_BASES);
    i += stack_seq->Load_Memory(&file[i]);

	Color prefixColor = 0;
	Color suffixColor = 0;
	if (hasPrefixSuffix & 0x1) {
		memcpy(&prefixColor, &file[i], sizeof(Color));
        assert( prefixColor != 0 );
        //assert( prefixColor != currentPartitionIndex );
		i += sizeof(Color);
	}
	if (hasPrefixSuffix & 0x2) {
		memcpy(&suffixColor, &file[i], sizeof(Color));
        assert( suffixColor != 0 );
        //assert( suffixColor != currentPartitionIndex );
		i += sizeof(Color);
	}

    /*// DEBUG: print loom entry to stdout
    printf("%01d %01d ", (hasPrefixSuffix & 0x1), (hasPrefixSuffix & 0x2));
    stack_seq->Print(); fputc('\n', stdout);*/

    convertSubsequenceToKmersToPrenodes(stack_seq, hashtable, kmer_length, hasPrefixSuffix, prefixColor, suffixColor);
    readIndex++;
  }
  
  //printf("Reading Loom file %s\t%i subsequences found.\n", filename, readIndex);

  return readIndex;
}

