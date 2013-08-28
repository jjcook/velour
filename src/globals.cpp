//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// globals.cpp
//
//   global variables
//

#include "types.h"

phase_e     g__PHASE                = PHASE_BUILD_KMER_GRAPH;
bool        g__PSEUDO_NODES_PRESENT = false;

uint64_t    g__READ_COUNT           = 0;

KmerGraph  *g__KG_HASHTABLE         = NULL;

SeqGraph   *g__SG_HASHTABLE         = NULL;

NodeAllocators *    g__NODE_ALLOCATORS = NULL;
KmerNodeAllocator * g__KMERNODE_ALLOCATOR = NULL;
SeqNodeAllocator *  g__SEQNODE_ALLOCATOR = NULL;

FILE       *g__VELVET_SEQUENCES_FILE = NULL;



// command line
//
char *      g__WORK_BASE_DIRECTORY  = NULL; //           base directory
char *      g__WORK_LOOM_DIRECTORY  = NULL; //  loom files subdirectory
char *      g__WORK_QUILT_DIRECTORY = NULL; // quilt files subdirectory
char *      g__WORK_INBOX_ROOT_DIRECTORY = NULL; //  inbox root subdirectory

unsigned    g__FULLKMER_LENGTH      = 0;
unsigned    g__MINIKMER_LENGTH      = 0;

bool        g__PARTITIONING         = false;

unsigned    g__PARTITION_COUNT      = 0;
unsigned    g__PARTITION_INDEX      = 0;
bool        g__PARTITION_NOEMIT     = false;

bool        g__FULL_STATISTICS      = false;

bool        g__COMBINING            = false;

bool        g__SLICING              = true;

bool        g__DOTGRAPH             = false;

size_t      g__MEMORY_FOOTPRINT_LIMIT = 0;

bool        g__NO_TIP_CLIPPING      = false;

bool        g__MINIMIZE_FOOTPRINT   = false;

double      g__COVCUTOFF_MIN = 0.0;
double      g__COVCUTOFF_MAX = 0.0;

bool        g__BUBBLE_POPPING = false;

// parse an unsigned value from string
unsigned parseUnsigned(char *h)
{
  long long_value;
  char *end_ptr;

  errno = 0;
  long_value = strtol(h, &end_ptr, 0); // automatically handle hexadecimal etc

  if( errno == ERANGE ) {
      fprintf(stderr, "Error parsing parameter '%s' -- out of range.\n", h);
      exit(EXIT_FAILURE);
  } else if( long_value > 0 && ((unsigned long) long_value) > UINT_MAX ) {
      fprintf(stderr, "Error parsing parameter '%s' -- too large.\n", h);
      exit(EXIT_FAILURE);
  } else if( long_value < 0 ) {
      fprintf(stderr, "Error parsing parameter '%s' -- negative value.\n", h);
      exit(EXIT_FAILURE);
  } else if( end_ptr == h ) {
      fprintf(stderr, "Error parsing parameter '%s' -- invalid numeric input.\n", h);
      exit(EXIT_FAILURE);
  } else if( *end_ptr != '\0' && *end_ptr != ' ' ) {
      fprintf(stderr, "Error parsing parameter '%s' -- extra characters proceeding input.\n", h);
      exit(EXIT_FAILURE);
  }

  return ((unsigned) long_value);
}

// extract partition index from input filename
void initializePartitionIndexFromInputFilename(file_object_vector *file_objects)
{
    char filename[PATH_MAX+1];
    strncpy(filename, file_objects->front().filename, PATH_MAX);
    filename[PATH_MAX] = '\0';
 
    char *index_start = strrchr(filename, '-') + 1;  // increment to start at next character
    if (index_start == NULL) {
        fprintf(stderr, "ERROR: strrchr() of input filename\n");
        assert(false);
        exit(EXIT_FAILURE);
    }

    char *index_stop = strchr(index_start, '.');
    if (index_stop == NULL) {
        fprintf(stderr, "ERROR: strchr() of input filename\n");
        assert(false);
        exit(EXIT_FAILURE);
    }
    *index_stop = '\0';

    g__PARTITION_INDEX = parseUnsigned(index_start);

    printf("Current partition index: %u\n", g__PARTITION_INDEX);

    // if the first partition, create the inbox subdirectories
    if (g__PARTITION_INDEX == 1) {
        for (unsigned i=1 ; i <= g__PARTITION_COUNT ; ++i) {
            char directory[PATH_MAX+1];
            sprintf(directory, "%s/inbox-for-%u", g__WORK_INBOX_ROOT_DIRECTORY, i);
            if (mkdir(directory, 0777) == -1) {
                if (errno != EEXIST) {
                    fprintf(stderr, "ERROR: failed to make directory: %s\n", directory);
                    perror("REASON: ");
                    exit(EXIT_FAILURE);
                }
                // TODO: delete directory contents?
            }
        }
    }
}

