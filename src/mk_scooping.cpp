//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// mk_scooping.cpp
//
// partitioning method: scoop the mini-kmer graph
//

#include "types.h"
#include "minikmer.h"

#include "distance.h"

#include <deque>
#include <math.h>
#include <set>

static uintptr_t numUncoloredKmers;

//
// version-dependent jazz
//

Color computeKmerColor(Kmer fullKmer)
{
  return getMiniKmerColor( computeMiniKmer(fullKmer, g__FULLKMER_LENGTH, g__MINIKMER_LENGTH) );
}

static void colorMiniKmerGraph(void);

void initializeMiniKmer(void)
{
  MINIKMER_GRAPH_SIZE = 1UL << (2 * g__MINIKMER_LENGTH);
  miniKmerGraph = (struct minikmer_node_s *) calloc(MINIKMER_GRAPH_SIZE, sizeof(struct minikmer_node_s));

  COLORMIN = 1;
  COLORMAX = COLORMIN + g__PARTITION_COUNT - 1;

  numUncoloredKmers = MINIKMER_GRAPH_SIZE / 2;

  printf("MKG Coloring mini-kmer graph: "); fflush(stdout);
  colorMiniKmerGraph();

  printf("MKG Initialization:  mini-k = %d  full-k = %d  colors = %u  totalNodes = %u\n",
      g__MINIKMER_LENGTH, g__FULLKMER_LENGTH, g__PARTITION_COUNT, MINIKMER_GRAPH_SIZE );
}

void destroyMiniKmer(void)
{
  free(miniKmerGraph);
}


typedef std::deque<Kmer> kmer_worklist_t;

static inline bool
colorAndQueueNeighbors(Color color, Kmer seed, kmer_worklist_t *worklist)
{
  bool foundUncolored = false;

  // some direction
  for(int i=0; i < 4; ++i)
  {
    Kmer neighbor = maskKmer((seed << 2) | i, g__MINIKMER_LENGTH);
    Kmer canon = canonicalKmer(neighbor, g__MINIKMER_LENGTH);

    if( isMiniKmerUncolored(canon) )
    {
      foundUncolored = true;

      setMiniKmerColor(canon, color);
      -- numUncoloredKmers;

      worklist->push_front(canon);
    }
  }

  // other direction
  for(int i=0; i < 4; ++i)
  {
    Kmer neighbor = maskKmer((seed >> 2) | (((Kmer)i) << (2*(g__MINIKMER_LENGTH-1))), g__MINIKMER_LENGTH);
    Kmer canon = canonicalKmer(neighbor, g__MINIKMER_LENGTH);

    if( isMiniKmerUncolored(canon) )
    {
      foundUncolored = true;

      setMiniKmerColor(canon, color);
      -- numUncoloredKmers;

      worklist->push_front(canon);
    }
  }

  return foundUncolored;
}

static void colorMiniKmerGraph(void)
{
  Kmer colorSeeds[COLORMAX+1];
  kmer_worklist_t colorWorklist[COLORMAX+1];

  printf("selecting seeds... "); fflush(stdout);

  // select random seeds
  for(Color c=COLORMIN; c <= COLORMAX; ++c)
  {
    Kmer seed;
    unsigned seedDistance = 0;
    unsigned fudge = (g__PARTITION_COUNT <= 64 ? 1 : 2);
    do
    {
      // select an uncolored node
      do {
        seed = random() % MINIKMER_GRAPH_SIZE;
        seed = canonicalKmer(seed, g__MINIKMER_LENGTH);
      } while( isMiniKmerColored(seed) );

      // check against existing seeds for minimum distance requirement
      unsigned minDistance = g__MINIKMER_LENGTH+1;
      for(unsigned j=COLORMIN; j < c; ++j)
      {
        unsigned thisDistance = getDistance(seed, colorSeeds[j], g__MINIKMER_LENGTH);
        minDistance = min(minDistance, thisDistance);
      }
      seedDistance = minDistance;
#ifdef VERBOSE
      printf("MKG Candidate seed 0x%08llx of min-distance %u\n", seed, seedDistance);
#endif
    } while( seedDistance < (g__MINIKMER_LENGTH-fudge) ); // FIXME: what value to use?

#ifdef VERBOSE
    printf("MKG Selected seed 0x%08llx\n", seed);
#endif
    colorSeeds[c] = seed;

    setMiniKmerColor(seed, c);
    -- numUncoloredKmers;

    colorWorklist[c].push_front(seed);
  }

  printf("coloring... "); fflush(stdout);

  // iterate until all canonical nodes colored
  unsigned iteration = 0;
  while( numUncoloredKmers > 0 )
  {
    bool worklistsEmpty = true;
    for(Color c=COLORMIN; c <= COLORMAX; ++c)
    {
      bool foundUncolored = false;
      while (!colorWorklist[c].empty() && !foundUncolored)
      {
        worklistsEmpty = false;

        Kmer nextKmer = colorWorklist[c].back();
        colorWorklist[c].pop_back();

        foundUncolored = colorAndQueueNeighbors(c, nextKmer, &(colorWorklist[c]));
      }
    }
    assert( worklistsEmpty == false && "Infinite loop!?" );
    ++ iteration;
  }

  printf("done.\n"); fflush(stdout);
}

void printMiniKmerStats(void)
{
  printMiniKmerGraphAdjacencyStats();
}

