//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

#ifndef _TYPES_H_
#define _TYPES_H_

#define VERSION_NUMBER 0
#define RELEASE_NUMBER 0
#define UPDATE_NUMBER 0

// *********************************
// *****  CODE PATH SELECTION  *****
// *****  (undef to disable)   *****
// *********************************
#ifndef OVERRIDE_CODEPATH

#define SMALL_NODES
//#define UNIQUE
//#define VELOUR_TBB 1
//#define USE_TBB_ALLOC 1
//#define VERIFY
#define VELVET_EMULATION

#endif // OVERRIDE_CODEPATH
// ************************************
// ***** END CODE PATH SELECTION  *****
// ************************************

#define __STDC_FORMAT_MACROS
#define __STDC_LIMIT_MACROS

#include <algorithm>
#include <assert.h>
#include <climits>
#include <fcntl.h>
#include <deque>
#include <dirent.h>
#include <errno.h>
#include <inttypes.h>
#include <libgen.h>
//#include <limits.h>
#include <limits>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <time.h>
#include <unistd.h>
#include <utility>
#include <vector>

#ifdef VELOUR_TBB
#  include <tbb/atomic.h>
#  include <tbb/blocked_range.h>
#  include <tbb/cache_aligned_allocator.h>
#  include <tbb/concurrent_vector.h>
#  include <tbb/parallel_for.h>
#  include <tbb/pipeline.h>
#  include <tbb/queuing_mutex.h>
#  include <tbb/scalable_allocator.h>
#  include <tbb/task_scheduler_init.h>
#  include <tbb/task_scheduler_observer.h>
#  include <tbb/tbb_thread.h>
#  include <tbb/tick_count.h>
#endif

// used to get sizeof member variables,
//   e.g. sizeof( memberType( &A::m ) )
template< typename C, typename T >
T memberType( T C::* ) ;

/*# ifdef _MSC_VER
# define inline __inline
# define DECLARE_ALIGNED( type, var, n ) __declspec(align(n)) type var
# else gcc
# define DECLARE_ALIGNED( type, var, n ) type var __attribute__((aligned(n)))
# endif*/

#ifdef VELOUR_TBB
#  define ATOMIZE(x) (*(reinterpret_cast< tbb::atomic<typeof((x))> * >(&(x))))
#else
#  define ATOMIZE(x) ((x))
#endif

#define ATOMIC_LOAD(x) ATOMIZE((x))
#define ATOMIC_STORE(x) ATOMIZE((x))

#ifdef VELOUR_TBB
#  define PACKED
#else
#  define PACKED __attribute__((packed))
#endif

typedef uint16_t threadid_t;
//typedef uint8_t threadid_t;
//static const threadid_t THREADID_MAX = std::numeric_limits<threadid_t>::max()-1; // TODO FIXME: move when not using #define max
static const threadid_t THREADID_MAX = 254; // TODO FIXME: constant

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

#include "histogram.h"

// TODO: replace define'd constants with static const int etc.?

#ifdef VERIFY
#define VERIFY_ASSERT(x) assert(x)
#else
#define VERIFY_ASSERT(x) ((void)0)
#endif

enum {ADENINE=0, CYTOSINE, GUANINE, THYMINE};
#define COMPLEMENT(x) ((x)^0x3)

#define CATEGORIES 2
typedef short Category;

typedef unsigned char byte_t;
typedef int64_t vid_t;

typedef unsigned char Nucleotide;
#include "kmer.h" // FIXME: temporary

#define PRIVATE_BUCKETS 128

// Hash table structures
#define HASH_BITS    24

// Input file types -- NOTE: keep the enum and strings consistent
enum {FASTQ=1, FASTA, GERALD, ELAND, FASTQ_GZ, FASTA_GZ, MAQ_GZ, PREGRAPH, LOOM, QUILT, BUCKET, NOFORMAT };

#define NO_SUCCESSORS (-1)
#define MULTIPLE_SUCCESSORS (-2)

#ifndef SMALL_NODES
#define LEFT_CONNECT_OFFSET 4   // left: top 4 bits, right: bottom 4 bits
#define LEFT(c) (c >> LEFT_CONNECT_OFFSET)
#define RIGHT(c) (c & 0xf)
#endif // SMALL_NODES

#define GO_LEFT 0
#define GO_RIGHT 1

#ifdef COLOR_8BIT
typedef uint8_t Color;
typedef uint8_t color_t;
typedef uint32_t four_color_t;
#else
typedef uint16_t Color;
typedef uint16_t color_t;
typedef uint64_t four_color_t;
#endif

typedef int8_t counter_t;
typedef uint32_t four_counter_t;
typedef uint64_t eight_counter_t;

#define CLIP_SINGLECOPY_COUNTER_VALUE 126
#define MAX_COUNTER_VALUE 125    // anything less than 127 should be fine
//#define MERGED_VALUE ((((((127 << 8) | 127) << 8) | 127) << 8) | 127)

static inline counter_t cnorm(counter_t x)
{
    counter_t a = abs(x);
    if (a == CLIP_SINGLECOPY_COUNTER_VALUE) {
        return 1;
    } else {
        return a;
    }
}

typedef enum {
  PHASE_BUILD_KMER_GRAPH=0,
  PHASE_OBSERVE_SEQUENCES,
  PHASE_DISTRIBUTE_SEQUENCES
} phase_e;

class KmerGraph;
class SeqGraph;
struct PrivateSeqGraphSet;

#ifdef VELOUR_TBB
extern __thread unsigned tls_thread_index;
extern bool *thread_aborts[THREADID_MAX+2];
extern __thread bool *tls_thread_abort;
extern __thread PrivateSeqGraphSet *tls_private_sgraph_set;
#endif // VELOUR_TBB


// TODO: move these functions to their own header file
//
// fast but obscure bit-level functions
//

// check if is power of 2
//
// source: http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
//
template<typename T>
bool is_power_of_2(T v)
{
  return ( !(v & (v - 1)) ) && v;
}

// check if is power of 4
//
// source: jjcook
//
static inline bool is_power_of_4(unsigned v)
{
  bool po2 = is_power_of_2(v);
  assert( sizeof(unsigned) == 4 );
  return (v & 0x55555555UL) && po2;
}

// find the log base 2 of a power of 2
//
// source: http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogDeBruijn
//
static inline uint32_t log2_of_power_of_2(uint32_t v)
{
  static const int MultiplyDeBruijnBitPosition[32] =
  {
    0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
  };

  return MultiplyDeBruijnBitPosition[static_cast<uint32_t>(v * 0x077CB531UL) >> 27];
}

// round up to the next power of 2
//   (note: if already a power of 2, does not change it)
//
// source: http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
//
template<typename T>
static inline T round_up_power_of_2(T v)
{
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;
  v++;
  return v;
}

template<>
inline unsigned round_up_power_of_2<unsigned>(unsigned v)
{
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

template<>
inline long unsigned int round_up_power_of_2<long unsigned int>(long unsigned int v)
{
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

template<unsigned value> class STATIC_ROUND_UP_POWER_OF_2
{
    public:
        enum {
            V0 = value - 1,
            V1 = V0 | (V0 >> 1),
            V2 = V1 | (V1 >> 2),
            V4 = V2 | (V2 >> 4),
            V8 = V4 | (V4 >> 8),
            V16 = V8 | (V8 >> 16),
            RESULT = V16 + 1
        };
};

extern phase_e     g__PHASE;
extern bool        g__PSEUDO_NODES_PRESENT;
extern uint64_t    g__READ_COUNT;
extern KmerGraph  *g__KG_HASHTABLE;
extern SeqGraph   *g__SG_HASHTABLE;
extern FILE       *g__VELVET_SEQUENCES_FILE;

extern char *      g__WORK_BASE_DIRECTORY;
extern char *      g__WORK_LOOM_DIRECTORY;
extern char *      g__WORK_QUILT_DIRECTORY;
extern char *      g__WORK_INBOX_ROOT_DIRECTORY;
extern unsigned    g__FULLKMER_LENGTH;
extern unsigned    g__MINIKMER_LENGTH;
extern bool        g__PARTITIONING;
extern bool        g__LOOMING;
extern bool        g__COMBINING;
extern bool        g__QUILTING;
extern unsigned    g__PARTITION_COUNT;
extern unsigned    g__PARTITION_INDEX;
extern bool        g__PARTITION_NOEMIT;
extern bool        g__FULL_STATISTICS;

extern bool        g__SLICING;
extern bool        g__DOTGRAPH;

extern size_t      g__MEMORY_FOOTPRINT_LIMIT;

extern unsigned    g__PGDIST_PARTITIONS;
extern unsigned    g__PGDIST_FILTER;

extern bool        g__NO_TIP_CLIPPING;
extern bool        g__MINIMIZE_FOOTPRINT;

extern double      g__COVCUTOFF_MIN;
extern double      g__COVCUTOFF_MAX;

extern bool        g__BUBBLE_POPPING;

#define KMER_APPEND(kmer, base, double_kmer_length)  \
  (((kmer) >> 2) | (((Kmer)(base)) << ((double_kmer_length) - 2)))

#define KMER_PREPEND(kmer, base, double_kmer_length, mask)  \
  ((((kmer) << 2) | ((Kmer)(base))) & (mask))

#define KMER_GET_HEAD_BASE(kmer, kmer_length) \
    getNucleotide(kmer, 0)

#define KMER_GET_TAIL_BASE(kmer, kmer_length) \
    getNucleotide(kmer, kmer_length - 1)

// -----
// -- allocators.cpp
// -----
void velour_alloc_init(void);
void velour_alloc_done(void);
void * velour_alloc(size_t sz);
void * tls_velour_alloc(size_t sz);
#define velour_free(v)  // does nothing.  Wait for the bulk free
void * velour_realloc(void * old_ptr, size_t old_sz, size_t new_sz);
// end allocators

class KmerGraph;
class KmerNode;
class SeqGraph;
class SeqNode;

template<typename NodeType>
void gdb_list_nodelist(std::deque<NodeType*> &nodelist)
{
    for (typename std::deque<NodeType*>::iterator nit = nodelist.begin(); nit != nodelist.end(); ++nit) {
        NodeType *node = *nit;
        printf("node: %p  claim_tid: %u\n", node, node->claim_tid);
    }
    fflush(stdout);
}

#include "kmerNode.h"
#include "seqNode.h"
#include "flags.h"

#ifdef VELOUR_TBB
typedef tbb::concurrent_vector<SeqNode *> flow_nodelist_t;
#else
typedef std::deque<SeqNode *> flow_nodelist_t;
#endif

extern size_t g__peakLiveMemory;

#include "component.h"


#include "node_allocators.h"
extern NodeAllocators *g__NODE_ALLOCATORS;
extern KmerNodeAllocator *g__KMERNODE_ALLOCATOR;
extern SeqNodeAllocator *g__SEQNODE_ALLOCATOR;
#ifdef VELOUR_TBB
extern __thread TLS_KmerNodeAllocator *tls_kmernode_allocator;
extern __thread TLS_SeqNodeAllocator *tls_seqnode_allocator;
#endif // VELOUR_TBB


#include "kmerGraph.h"
#include "seqGraph.h"

#include "sequence_inlined.h"

//#include "sg_utility.h"


typedef struct file_object {
  off_t length; 
  char *filename; 
  int filetype;
  int fileindex;
  Category cat;
} file_object_t;

typedef std::vector<file_object_t> file_object_vector;


static const unsigned MAX_READ_LENGTH = 65534;
static const unsigned MAX_READ_HEADER_LENGTH = 1000;


// *****************************************************
// ***  MONOLITHIC HEADER FILE FORWARD DECLARATIONS  ***
// *****************************************************

/*
// -----
// -- graphConstruction.cpp
// -----
  //extern node_t **NODES;
extern unsigned NODE_ALLOCATED;

void make_graph(kg_node_t **hashtable, int kmer_length);
  //node_t * newNode(void);
  //node_t * get_node(unsigned i);
void print_node(node_t *node);
void dump_graph(const char *filename, unsigned kmer_length);
*/

// kg_stats.cpp
//void kgstat_clusters(void);

// minikmer.cpp
void printMiniKmerStats(void);

// -----
// -- parsing.cpp
// -----
extern Nucleotide BASE_MAP[256]; // fixme: statically initialize this table
extern char CHAR_BASE_MAP[4];
extern char CONNECTION_COUNT_MAP[16];
extern char CONNECTION_MAP[16];
extern const char *FILE_TYPES[];

void initializeBaseMap(void);
void load_loom_files(file_object_vector &);
void process_files(file_object_vector &);

// -----
// -- partition.cpp
// -----
void initializePartitionerPreObservation(void);
void initializePartitionerPostObservation(void);
void destroyPartitioner(void);
void distributeSequence(char *seq, int kmer_length);
void printPartitionerStatistics(void);


// -----
// -- preGraphConstruction.cpp
// -----
//extern unsigned PRENODE_ALLOCATED;
//extern unsigned PRENODE_UNIQUE;

void convertSequenceToKmersToPrenodes(char *seq, KmerGraph* hashtable, int kmer_length);
void convertSubsequenceToKmersToPrenodes(Sequence *seq, KmerGraph* hashtable, int kmer_length, unsigned hasPrefixSuffix, Color prefixColor, Color suffixColor);
void print_hashtable_histogram(KmerGraph* hashtable);


// -----
// -- pg_singleCopyRemoval.cpp
// -----
//void remove_single_copy_kmers(kg_node_t **hashtable, int kmer_length);

// -----
// -- pg_tipClipping.cpp
// -----
void remove_tips(KmerGraph* hashtable);

// -----
// -- sg_tipClipping.cpp
// -----
void sg_remove_tips(SeqGraph *, bool silent=false);
void sg_nodelist_remove_tips(SeqGraph *graph, bool silent, flow_nodelist_t *nodelist);
void sg_parallel_nodelist_remove_tips(SeqGraph *graph, bool silent, flow_nodelist_t *nodelist);

// sg_bubblePopping.cpp
void sg_pop_bubbles(SeqGraph *graph, bool silent=false);
void sg_nodelist_pop_bubbles(SeqGraph *graph, bool silent, flow_nodelist_t *nodelist);

// sg_covcutoff.cpp
void sg_covcutoff(SeqGraph *graph, bool silent=false);

// -----
// -- quilt.cpp
// -----
void quilt_files(SeqGraph *, file_object_vector&);

// -----
// -- utility.cpp
// -----
//extern KmerGraph& HASHTABLE;

template<typename T> bool is_power_of_2(T v);
bool is_power_of_4(unsigned v);
unsigned log2_of_power_of_2(unsigned v);
//unsigned round_up_power_of_2(unsigned v);

void report_unrecoverable_error(void);

#ifdef SMALL_NODES
//void get_neighbor_counts(KmerNode *node, unsigned *l_count, unsigned *r_count);
void verify_node(KmerNode *node, KmerGraph* hashtable, unsigned kmer_length);
#else
void verify_node(kg_node_t * node, unsigned kmer_length);
#endif
int valid_single_successor(counter_t *counters);
//int get_node_unique(sg_node_t *node, int count);
counter_t get_counter_sum(counter_t *counters);
counter_t get_counter_max(counter_t *counters);
#ifndef SMALL_NODES
int find_connection_index(kg_node_t *one, kg_node_t *two, bool left_not_right);
void add_neighbor(kg_node_t *node, unsigned index, kg_node_t *rnode, counter_t count, bool left_not_right);
void add_right_neighbor(kg_node_t *node, unsigned index, kg_node_t *rnode, counter_t count);
void remove_right_neighbor(kg_node_t *node, unsigned index);
void remove_left_neighbor(kg_node_t *node, unsigned index);
#endif

// -----
// -- velour.cpp
// -----


// split.cpp
#include "split.h"

// graphviz.cpp
void emit_graphviz(KmerGraph *graph, char *filename);
void emit_graphviz(SeqGraph *graph, char *filename);
void emit_scoop_graphviz(SeqGraph *graph, SeqNode *root, intptr_t nodes, char *filename);

void kg_stat_components(KmerGraph *graph, FILE *output);
void sg_stat_components(SeqGraph *graph, FILE *output);

// slicing.cpp
extern uintptr_t g__SLICE_NODE_COUNT;
//void slice_component(SeqGraph *graph, Component<SeqGraph, SeqNode>& component);
void slice_graph(SeqGraph *graph);

// slicing2.cpp
extern uintptr_t g__SLICE2_FINAL_NODE_COUNT;
extern uintptr_t g__SLICE2_NODE_COUNT;
void slice2_graph(SeqGraph *graph, unsigned currentPartitionIndex);
void slice2_nodelist(SeqGraph *graph, unsigned currentPartitionIndex, flow_nodelist_t *nodelist);

// pregraph_partitining.cpp
void pregraph_partitioning(SeqGraph *sgraph, char *metis_filename);

#endif /* _TYPES_H_ */

