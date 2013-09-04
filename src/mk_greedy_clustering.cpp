//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// mk_greedy_clustering.cpp
//
//   partitioning method: greedy incremental cluster growth
//     - defers reads that don't intersect, then rescans them

#include "types.h"
#include "minikmer.h"

//#undef min
//#undef max

#include <algorithm>
#include <deque>
#include <math.h>
#include <set>

//#define MINIKMER_MODE 1

//static unsigned colorLimit;
static uint64_t *colorCount;
static uint64_t *baseCount;

static unsigned usedColors;
static Color nextColor;

static uint64_t minikmers_colored = 0;

#define NUM_SEEDS 100
#define SLOP_FACTOR 1

static uint64_t MAX_UNATTACHED; // dynamically initialized based on the free memory available
static uint64_t STRING_ALLOC_SIZE; // dynamically initialized based on the free memory available
static char **unattached;
static uint64_t num_unattached;
static char *bump_begin;
static char *bump_head;
static char *bump_end;

bool cleanup = false;
bool first = true;
void small_cluster_stats();
void merge_clusters();

unsigned *adjacencyMatrix;
unsigned g__TEMP_PARTITION_COUNT = 0;

#define ADJACENCY(i,j) (adjacencyMatrix[(i*COLORMAX)+j])
void
increment_adjacency(int i, int j) {
  if (i > j) {
	 int swap = i;	 i = j; j = swap;
  }
  ADJACENCY(i,j) ++;
}

unsigned
get_adjacency(int i, int j) {
  if (i > j) {
	 int swap = i;	 i = j; j = swap;
  }
  return ADJACENCY(i,j);
}

#ifdef MINIKMER_MODE
uint64_t greedyHash(Kmer fullKmer) { return computeMiniKmer(fullKmer, g__FULLKMER_LENGTH, g__MINIKMER_LENGTH >> 1); }
#else
#error // not fixed for LARGE_KMERS
uint64_t
greedyHash(Kmer fullKmer) {
  uint64_t hash = fullKmer ^ (fullKmer >> 21) ^ (fullKmer << 7);
  hash = hash ^ (fullKmer >> 13);
  return hash & (MINIKMER_GRAPH_SIZE - 1);
}
#endif

//
// version-dependent jazz
//
Color 
computeKmerColor(Kmer fullKmer) {
  uint64_t hash = greedyHash(fullKmer);
  Color c = getMiniKmerColor(hash); 
  assert((c != 0) && (c <= g__TEMP_PARTITION_COUNT));
  return c;
}

void 
initializeMiniKmer(void) {
  // global variables
#ifdef ARCH_64BIT
  MINIKMER_GRAPH_SIZE = 1ULL << g__MINIKMER_LENGTH;
#else
  MINIKMER_GRAPH_SIZE = 1UL << g__MINIKMER_LENGTH;
#endif

  if ((MINIKMER_GRAPH_SIZE * sizeof(struct minikmer_node_s)) >= g__MEMORY_FOOTPRINT_LIMIT) { // XXX: for STRING_ALLOC_SIZE below
      printf("GREEDYCLUSTERING -- ERROR: mini-kmer graph is larger or equal to specified memory footprint limit.\n");
      exit(1);
  }

  miniKmerGraph = (struct minikmer_node_s *) calloc(MINIKMER_GRAPH_SIZE, sizeof(struct minikmer_node_s));

#ifdef COLOR_8BIT
  if (g__PARTITION_COUNT > 254) { printf("ERROR: greedy clustering: too many partitions...\n"); exit(1); }
  g__TEMP_PARTITION_COUNT = 254;
#else
  if (g__PARTITION_COUNT > 16384) { printf("ERROR: greedy clustering: too many partitions... bad.\n"); exit(1); }
  g__TEMP_PARTITION_COUNT = max(254, 4 * g__PARTITION_COUNT - 2);
#endif

  printf("partitioning into %d initial partitions\n", g__TEMP_PARTITION_COUNT);

#ifdef MINIKMER_MODE
  printf("  greedy clustering: mini-kmer mode\n");
#else
  printf("  greedy clustering: random hash mode\n");
#endif

  COLORMIN = 1;
  COLORMAX = COLORMIN + g__TEMP_PARTITION_COUNT - 1;
  nextColor = 1;
  usedColors = 0;

  adjacencyMatrix = (unsigned *) calloc((COLORMAX+1)*(COLORMAX+1), sizeof(unsigned));

  // static variables
  //colorLimit = MINIKMER_GRAPH_SIZE / g__TEMP_PARTITION_COUNT;
  colorCount = (uint64_t *) calloc(COLORMIN + g__TEMP_PARTITION_COUNT, sizeof(uint64_t));
  baseCount  = (uint64_t *) calloc(COLORMIN + g__TEMP_PARTITION_COUNT, sizeof(uint64_t));

  printf("MKG Initialization:  mini-k = %d  full-k = %d  colors = %u  totalNodes = %"PRIuPTR"\n",
      g__MINIKMER_LENGTH, g__FULLKMER_LENGTH, g__TEMP_PARTITION_COUNT, MINIKMER_GRAPH_SIZE );

  // use all the space we can to buffer read sequences
  STRING_ALLOC_SIZE = g__MEMORY_FOOTPRINT_LIMIT - (MINIKMER_GRAPH_SIZE * sizeof(struct minikmer_node_s));
  STRING_ALLOC_SIZE -= (256 * 1024 * 1024); // misc room for the system and misc velour memory
  STRING_ALLOC_SIZE -= (STRING_ALLOC_SIZE >> 3); // make room for 'unattached' pointers

  MAX_UNATTACHED = STRING_ALLOC_SIZE / 105;

  unattached = (char **) malloc(MAX_UNATTACHED * sizeof(char *));
  if (unattached == NULL) {
      fprintf(stderr, "Allocation error in mk_greedy_clustering: unattached wanted size %"PRIu64".\n", MAX_UNATTACHED); fflush(stdout); fflush(stderr);
      assert(false && "Allocation error.");
      exit(EXIT_FAILURE);
  }
  num_unattached = 0;
  bump_begin = bump_head =  (char *)  malloc(STRING_ALLOC_SIZE *sizeof(char));
  bump_end = bump_head + STRING_ALLOC_SIZE;
}

void 
destroyMiniKmer(void) {
  free(miniKmerGraph);
  free(colorCount);
  free(baseCount);
  if (unattached != NULL) { free(unattached); }
  free(bump_begin);
}

//
// simple non-helper functions
//
uint64_t read_count = 0;

void 
printMiniKmerStats_internal(bool force) {
  if (!force) {
	 read_count ++;
	 if (read_count % 3000000) {
		return;
	 }
  }
  
  if (first) {
	 first = false;
	 small_cluster_stats();
  }

  //printf("****************************** %d ******************************\n", read_count);
  // printMiniKmerGraphAdjacencyStats();
  // for (int i = 0 ; i < COLORMAX ; i ++) {
  // 	 for (int j = 0 ; j < COLORMAX ; j ++) {
  // 		printf("%d, ", ADJACENCY(i,j));
  // 	 }
  // 	 printf("0\n");
  // }
  // printf("############################## %d ##############################", read_count);
  uint64_t totalColored = 0;
  for (int i = 1 ; i <= COLORMAX ; i ++) {
	 //printf("%"PRIu64" ", colorCount[i]);
	 totalColored += colorCount[i];
  }
  //printf("\n");
  printf("  GREEDY: %"PRIu64"/%"PRIu64" (%.2f percent) unattached\n", num_unattached, totalColored+num_unattached, (100.0*num_unattached)/(totalColored+num_unattached));
  printf("  GREEDY: %"PRIu64"/%"PRIu64" (%.2f percent) unattached array full\n", num_unattached, MAX_UNATTACHED, (100.0*num_unattached)/(MAX_UNATTACHED+num_unattached));
#ifdef MINIKMER_MODE
  printf("  GREEDY: %"PRIu64"/%"PRIuPTR" (%.2f percent) minikmer graph full\n", minikmers_colored, MINIKMER_GRAPH_SIZE/2, (100.0*minikmers_colored)/(MINIKMER_GRAPH_SIZE/2));
#else
  printf("  GREEDY: %"PRIu64"/%"PRIuPTR" (%.2f percent) table full\n", minikmers_colored, MINIKMER_GRAPH_SIZE, (100.0*minikmers_colored)/MINIKMER_GRAPH_SIZE);
#endif
}

void 
printMiniKmerStats() {
  printMiniKmerStats_internal(true);
}

bool first_go_round = true;

// look for a color that is free 
static Color
getNextFreeColor(void) {
  assert(usedColors < (NUM_SEEDS * g__TEMP_PARTITION_COUNT));
  
  // Color initialNextColor = nextColor;
  // while (colorCount[nextColor] != 0) {
  //   if (nextColor < COLORMAX)
  //     ++ nextColor;
  //   else
  //     nextColor = COLORMIN;
  // 
  //   assert(nextColor != initialNextColor);
  // }

  Color c = nextColor;
  assert(c >= COLORMIN);
  assert(c <= COLORMAX);

  if (nextColor == COLORMAX) {
    nextColor = COLORMIN;
  } else {
    nextColor ++;
  }

  usedColors ++;
  return c;
}

unsigned leastUsedCount = 0;

// look for a color that is not heavily used
static Color
getLeastUsedColor(void) {
  leastUsedCount ++;
  assert(usedColors == (NUM_SEEDS * g__TEMP_PARTITION_COUNT));
  
  unsigned best = 0;
  uint64_t best_count = UINT64_MAX;
  for (int i = COLORMIN ; i < COLORMAX ; i ++) {
	 if (colorCount[i] < best_count) {
		best = i;
		best_count = colorCount[i];
	 }
  }

  assert(best >= COLORMIN);
  assert(best <= COLORMAX);

  return best;
}

// // increase fairness of round-robin?
// static void
// returnColor(Color oldColor)
// {
//   if( nextColor == COLORMIN && oldColor == COLORMAX )
//     nextColor = COLORMAX;
//   else if( (nextColor-1) == oldColor )
//     -- nextColor;
// }

namespace {
    struct lambda_gc_process
    {
        lambda_gc_process(uint64_t &h, uint64_t *hash, Color &last_color, Color &color, Color *colors, unsigned &colored) : h(h), hash(hash), last_color(last_color), color(color), colors(colors), colored(colored) {}
        void operator()(Kmer k, Kmer canon_k, bool new_sense_reversed, Nucleotide base, Nucleotide last_base, unsigned kmer_length, bool isFirstKmer, bool isLastKmer)
        {
            if (isFirstKmer) {
                // handle first kmer
                hash[0] = h = greedyHash(canon_k);
                colors[0] = color = getMiniKmerColor(h);
                if (color != 0) {
                    colored = 1;  last_color = color;
                }

                i = kmer_length;
                return;
            }

            ++i;
            hash[i - kmer_length] = h = greedyHash(canon_k);
            colors[i - kmer_length] = color = getMiniKmerColor(h);
            if (color != 0) {
                colored ++;   last_color = color;
            }
            return;
        }
        private:
            uint64_t &h;
            uint64_t *hash;
            Color &last_color;
            Color &color;
            Color *colors;
            unsigned &colored;

            int i;
    };
}

// foreach kmer in this read: observe to create a partition map
bool
processSequence(char *seq, int kmer_length) {
  uint64_t h, hash[MAX_READ_LENGTH];
  Color last_color, color, colors[MAX_READ_LENGTH];
  unsigned colored = 0;
  bzero(colors, MAX_READ_LENGTH * sizeof(Color));

  {
    Sequence_StackAllocated memory;
    Sequence *stack_seq = new (&memory) Sequence(Sequence::MAX_BASES);
    stack_seq->InitializeWithString(seq);

    lambda_gc_process functor(h, hash, last_color, color, colors, colored);
    sequence_process_kmers(stack_seq, kmer_length, functor);
  }
/*
  Kmer kmer = 0, anti_kmer = 0;
  int i;

  // read first Kmer
  for (i = 0 ; i < kmer_length ; ++ i) {
    char c = seq[i];
    if (c == 0) { return true; }  // end of string without a full Kmer

    Nucleotide base = BASE_MAP[(int) c];
    anti_kmer >>= 2;
    kmer <<= 2;
    assert(base <= 0x3);
    anti_kmer |= ((Kmer)base ^ 0x3) << 62;
    kmer |= base;
  }

  int double_kmer_length = kmer_length << 1;
  Kmer mask = (((Kmer)1) << double_kmer_length) - 1;  // a mask with 2*kmer_length bits

  Kmer rc_kmer = (anti_kmer >> (64 - double_kmer_length));
  bool sense_reversed = rc_kmer < kmer;  
  Kmer canonical_kmer = sense_reversed ? rc_kmer : kmer;

  // handle first kmer
  hash[0] = h = greedyHash(canonical_kmer);
  colors[0] = color = getMiniKmerColor(h);
  if (color != 0) {
	 colored = 1;  last_color = color;
  }

  // make each succeeding Kmer
  char c;
  while ((c = seq[i]) != 0) {
    ++ i;

    // read the next base and extend both the kmer and the anti_kmer
    Nucleotide base = BASE_MAP[(int)c];
    anti_kmer >>= 2;
    kmer <<= 2;
    if (base > 0x3) break; // stop at first 'unknown' base
    //assert(base <= 0x3);
    anti_kmer |= ((Kmer)base ^ 0x3) << 62;
    kmer |= base;
    kmer &= mask;

    Kmer rc_kmer = (anti_kmer >> (64 - double_kmer_length));
    bool sense_reversed = rc_kmer < kmer;  
    Kmer canonical_kmer = sense_reversed ? rc_kmer : kmer;

	 hash[i - kmer_length] = h = greedyHash(canonical_kmer);
	 colors[i - kmer_length] = color = getMiniKmerColor(h);
	 if (color != 0) { 
		colored ++;   last_color = color;
	 }
  }
*/
  int i = strlen(seq)+1; // jjcook

  int end = (i - kmer_length) + 1;

  // Special cases
  if (colored == 0) {
	 if ((usedColors < (NUM_SEEDS * g__TEMP_PARTITION_COUNT)) || cleanup) { // There is a free color, use it
		Color newColor = cleanup ? getLeastUsedColor() : getNextFreeColor();
		for (int j = 0 ; j < end ; j ++) {
		  // assert(getMiniKmerColor(hash[j]) == 0);
          if (getMiniKmerColor(hash[j]) == 0) { ++ minikmers_colored; }
		  setMiniKmerColor(hash[j], newColor); 
		  colorCount[newColor] ++;
		  baseCount[newColor] ++;
		}
		return true;
	 } 
	 return false;
  }

  // walk the color array associating each kmer with the closest color, first the
  // stuff up to the first color
  int j = 0;
  for ( ; j < end ; j ++) {
	 if ((color = colors[j]) != 0) {
		baseCount[color] ++;
		break;
	 }
  }
  for (int k = 0 ; k < j ; k ++) {
     if (getMiniKmerColor(hash[k]) == 0) { ++ minikmers_colored; }
	 setMiniKmerColor(hash[k], color);
	 colorCount[color] ++;
	 baseCount[color] ++;
  }
  // go between each pair of colored nodes
  while ((i = j + 1) < end) {
	 for (; i < end ; i ++) {
		if (colors[i] != 0) {  // found the other end
		  baseCount[colors[i]] ++;

		  // color the first nodes color of the node on that end.
		  // (divide the nodes by the inverse of the size of each clusters)
		  // float divisor = 1.0 - ((float)colorCount[color]) / (colorCount[color] + colorCount[colors[i]]);
		  // int half = j + (i - j) * divisor;  // not really half anymore.
		  int half = (colorCount[color] > colorCount[colors[i]]) ? j : i;  // give them all to the smaller cluster
		  for (j = j + 1 ; j < half ; j ++) {
             if (getMiniKmerColor(hash[j]) == 0) { ++ minikmers_colored; }
			 setMiniKmerColor(hash[j], color);
			 colorCount[color] ++;
			 baseCount[color] ++;
		  }
		  if (color != colors[i]) {
			 increment_adjacency(color, colors[i]);
		  }
		  color = colors[i];
		  for ( ; j < i ; j ++) {
             if (getMiniKmerColor(hash[j]) == 0) { ++ minikmers_colored; }
			 setMiniKmerColor(hash[j], color);
			 colorCount[color] ++;
			 baseCount[color] ++;
		  }
		}
	 }
	 if (i == end) {
		for (j = j + 1 ; j < i ; j ++) {
          if (getMiniKmerColor(hash[j]) == 0) { ++ minikmers_colored; }
		  setMiniKmerColor(hash[j], color);
		  colorCount[color] ++;
		  baseCount[color] ++;
		}
	 }
  }

/* TODO TODO TODO
#ifdef MINIKMER_MODE
  if (minikmers_colored == MINIKMER_GRAPH_SIZE) { g__abortObservation = true; }
#else
  if (minikmers_colored == MINIKMER_GRAPH_SIZE) { g__abortObservation = true; }
#endif
*/

  return true;
}

void
small_cluster_stats() {
  unsigned num_small = 0;
  for (int i = 0 ; i < COLORMAX ; i ++) {
	 if (colorCount[i] < 50) {
		num_small ++;
		printf("%"PRIu64" (%d): ", colorCount[i], i);

		unsigned total = 0;
		for (int j = 0 ; j < COLORMAX ; j ++) {
		  unsigned adj = get_adjacency(i,j);
		  if (adj != 0) printf("%d, ", adj);
		  total += adj;
		}
		printf("(%d)\n", total);
	 }
  }
  printf("num_small: %d\n", num_small);
}

int fullKmer_length;

void
observeSequence(char *seq, int kmer_length) {
  fullKmer_length = kmer_length;
  bool used = processSequence(seq, kmer_length);

  if (!used) {
	 unsigned len = strlen(seq);

     if( bump_head + len + 1 >= bump_end || num_unattached == MAX_UNATTACHED ) {
		 printf("  GREEDY: unattached buffer full, forcing future sequences into table.\n");
         cleanup = true;
         bool used = processSequence(seq, kmer_length);
         assert( used == true );
         return;
     }

	 char *copy = bump_head;
	 bump_head += len + 1;
	 assert(bump_head < bump_end);
	 assert(num_unattached < MAX_UNATTACHED);
	 unattached[num_unattached ++] = copy;
	 strcpy(copy, seq);
  } else {
	 printMiniKmerStats_internal(false);
  }
}

void 
completeObservation() {
  // printf("##############################  ##############################");
  // unsigned total = 0;
  // for (int i = 0 ; i < COLORMAX ; i ++) {
  // 	 for (int j = 0 ; j < COLORMAX ; j ++) {
  // 		unsigned adj = ADJACENCY(i,j);
  // 		printf("%d, ", adj);
  // 		total += adj;
  // 	 }
  // 	 printf("0\n");
  // }
  // printf("##############################  ##############################");
  // printf("ave adjacency: %d\n", total*2/COLORMAX);
  printf("##############################  ##############################");

  small_cluster_stats();

  for (int j = 0 ; j < 2 ; j ++) {
	 uint64_t wavefront = 0;
	 for (uint64_t i = 0 ; i < num_unattached ; i ++) {
		bool used = processSequence(unattached[i], fullKmer_length);
		if (!used) {
		  char *temp = unattached[wavefront];
		  unattached[wavefront] = unattached[i];
		  unattached[i] = temp;
		  wavefront ++;
		} else {
		  printMiniKmerStats_internal(false);
		}
	 }
	 num_unattached = wavefront;
	 printf("new unattached: %"PRIu64"\n", num_unattached);
	 small_cluster_stats();
	 cleanup = true;
	 printMiniKmerStats_internal(true);
  }

  printf("least used counter: %d\n", leastUsedCount);
  merge_clusters();

  // free this memory so that it doesn't get in the way of partition file buffers later
  free(unattached);
  unattached = NULL;
}



unsigned
cost_function(unsigned color_count, unsigned base_count, unsigned connections) {
//   unsigned cost = color_count
// 
// / (connections + 100);
  return 0;
}

int mapping_compare(const void *c1, const void *c2)
{
    const unsigned color1 = *(static_cast<const unsigned*>(c1));
    const unsigned color2 = *(static_cast<const unsigned*>(c2));

    int64_t diff = static_cast<int64_t>(baseCount[color1]) - static_cast<int64_t>(baseCount[color2]);
    // since types aren't compatible to just return diff
    if (diff > 0)
        return -1;
    else if (diff < 0)
        return 1;
    else
        return 0;
}

// this is a greedy merger which tries to minimize connections between clusters.
void
merge_clusters() {
  unsigned num_clusters = g__TEMP_PARTITION_COUNT;
  unsigned target_num_clusters = g__PARTITION_COUNT;
  
  bool *alive     = (bool *) alloca((COLORMAX+1) * sizeof(bool));
  bool *mergeable = (bool *) alloca((COLORMAX+1) * sizeof(bool));
  unsigned *mapping1 = (unsigned *) alloca((COLORMAX+1) * sizeof(unsigned));

  uint64_t total_count = 0, total_base = 0;
  uint64_t smallest_count = UINT64_MAX, smallest_base = UINT64_MAX;
  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
	 alive[i] = true;  mergeable[i] = true; mapping1[i] = i;
	 total_count += colorCount[i]; 
	 total_base  +=  baseCount[i]; 
	 if (smallest_count > colorCount[i])	 smallest_count = colorCount[i];
	 if (smallest_base  > baseCount[i])     smallest_base  = baseCount[i];
  }
  uint64_t count_threshold = (uint64_t) SLOP_FACTOR * total_count / g__PARTITION_COUNT;
  uint64_t base_threshold  = (uint64_t) SLOP_FACTOR * total_base  / g__PARTITION_COUNT;
  printf("thresholds are (%"PRIu64" / %"PRIu64") and smallest are (%"PRIu64" / %"PRIu64")\n", 
			count_threshold, base_threshold, smallest_count, smallest_base);

  // if the cluster is too large to merge even with the smallest task, mark it as such
  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
	 if ((colorCount[i] + smallest_count > count_threshold) || 
		  (baseCount[i]  + smallest_base  >  base_threshold)) {
		mergeable[i] = false;
		printf("cluster %d is too large to merge (%"PRIu64" / %"PRIu64")\n", i, colorCount[i], baseCount[i]);
	 }
  }

  printf("Clustering...\n");

  while (num_clusters > target_num_clusters) {
	 // find the "best" candidate for merging
	 uint64_t new_smallest_count = UINT64_MAX, new_smallest_base = UINT64_MAX;
	 unsigned best_i = 0, best_j = 0;
	 long long best_value = -10000000LL;
	 for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
		if (!alive[i] || !mergeable[i]) {
		  continue;
		}
		if ((colorCount[i] + smallest_count > count_threshold) || 
			 (baseCount[i]  + smallest_base  >  base_threshold)) {
		  mergeable[i] = false;
		  continue;
		}
		if (new_smallest_count > colorCount[i])	 new_smallest_count = colorCount[i];
		if (new_smallest_base  > baseCount[i])     new_smallest_base  = baseCount[i];

		for (int j = i + 1 ;  j <= COLORMAX ; j ++) {
		  if (!alive[j] || !mergeable[j] || 
				(colorCount[i] + colorCount[j] > count_threshold) || 
				(baseCount[i]  +  baseCount[j] >  base_threshold)) {
			 continue;
		  }
		  long long adj = (get_adjacency(i,j) * 100000LL);
		    // - (count_threshold - (colorCount[i] + colorCount[j]))
			 // - (base_threshold - (baseCount[i]  +  baseCount[j]));

		  if (adj >= best_value) {
			 best_value = adj;
			 best_i = i; best_j = j;
		  }
		}
	 }
		
	 if ((best_j == 0) && (best_i == 0)) {
		printf("WARNING: could only reduce down to %d clusters\n", num_clusters);
		g__PARTITION_COUNT = num_clusters;
		break;
	 }
	 // do the actual merge
	 assert(best_j > best_i);
	 mapping1[best_j] = best_i;
	 alive[best_j] = false;
	 colorCount[best_i] += colorCount[best_j];
	 baseCount [best_i] += baseCount [best_j];
     colorCount[best_j] = 0; // for mapping3 below
     baseCount[best_j] = 0;
	 num_clusters --;
	 
	 // update the known smallest
	 smallest_count = new_smallest_count;		smallest_base = new_smallest_base;
  }

  printf("Remapping...\n");

  // we want to map down to the first "target_num_clusters" colors
  unsigned next_color_to_use = COLORMIN;
  unsigned *mapping2 = (unsigned *) alloca((COLORMAX + 1) * sizeof(unsigned));
  mapping2[0] = 0;
  
  for (unsigned i = COLORMIN ; i <= COLORMAX ; i ++) {
	 unsigned x = i, y;
	 while ((y = mapping1[x]) != x) {
		x = y;
	 }
	 assert(x <= i);
	 mapping2[i] = (i == x) ? next_color_to_use ++ : mapping2[x];
	 assert(mapping2[i] < (COLORMIN + num_clusters));
	 assert(next_color_to_use <= (COLORMIN + num_clusters));
  }

/*  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
      printf("mapping2[%3d]=%3d", i, mapping2[i]);
      if (alive[i]) { printf(" **"); }
      puts("");
  }*/

  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
	 if (alive[i]) {
		printf("cluster %3d: %10"PRIu64" %10"PRIu64"\n", mapping2[i], colorCount[i], baseCount[i]);
        //if (alive[i]) { printf("**"); } puts("");
	 }
  }

  printf("*** sorting partitions: largest first ***\n");

  // reorder numbering to "shape" the partition size distribution: largest partitions first
  unsigned *sorted = (unsigned *) alloca((COLORMAX + 1) * sizeof(unsigned));
  sorted[0]=0;
  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
      if (alive[i]) {
          sorted[i] = i;
      } else {
          sorted[i] = 0;
      }
  }
  qsort(&sorted[COLORMIN], COLORMAX-COLORMIN+1, sizeof(unsigned), mapping_compare);
/*  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
      //if (alive[i]) {
        printf("sorted[%3d]=%3d", i, sorted[i]);
        if (alive[i]) { printf(" **"); }
        puts("");
      //}
  }*/

  // invert the mapping
  unsigned *inv_sorted = (unsigned *) alloca((COLORMAX + 1) * sizeof(unsigned));
  inv_sorted[0]=0;
  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
      inv_sorted[sorted[i]] = i;
  }

  // apply mapping2 to sorted list
  unsigned *sorted2 = (unsigned *) alloca((COLORMAX + 1) * sizeof(unsigned));
  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
      sorted2[i] = mapping2[sorted[i]];
  }

  // invert the mapping
  unsigned *inv_sorted2 = (unsigned *) alloca((COLORMAX + 1) * sizeof(unsigned));
  inv_sorted2[0]=0;
  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
      inv_sorted2[sorted2[i]] = i;
  }

/*  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
      //if (alive[i]) {
        printf("inv_sorted[%3d]=%3d", i, inv_sorted[i]);
        if (alive[i]) { printf(" **"); }
        puts("");
      //}
  }*/

  // translate mapping2 with new numbering
  unsigned *mapping3 = (unsigned *) alloca((COLORMAX + 1) * sizeof(unsigned));
  mapping3[0]=0;
  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
      mapping3[i] = inv_sorted2[mapping2[i]];
  }

  for (int i = COLORMIN ; i <= COLORMAX ; i ++) {
      if (sorted[i]) {
		printf("cluster %3d: %10"PRIu64" %10"PRIu64"\n", i, colorCount[sorted[i]], baseCount[sorted[i]]);
	 }
  }

  // actually do the remapping
  for (uintptr_t i = 0 ; i < MINIKMER_GRAPH_SIZE ; i ++) {
	 unsigned c = mapping3[miniKmerGraph[i].color];
	 //unsigned c = mapping2[miniKmerGraph[i].color];
	 assert(c < (COLORMIN + num_clusters));
     assert( miniKmerGraph[i].color == 0 || c > 0 );
     miniKmerGraph[i].color = c;
  }
  
	 // unsigned total = 0;
	 // unsigned max = 0, max_index = 0;
  	 // for (int j = 0 ; j < COLORMAX ; j ++) {
  	 // 	unsigned adj = get_adjacency(i,j);
	 // 	total += adj;
	 // 	if (adj > max) { max = adj; max_index = j; }
	 // }
	 // printf("%4d: %7d %8d %6d (%4d: %5d %5.2f)\n", i, colorCount[i], baseCount[i], total, max_index, max, (100.0*max)/total);
}
