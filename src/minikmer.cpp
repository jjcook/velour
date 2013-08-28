//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// minikmer.cpp
//
// common minikmer code

#include "types.h"
#include "minikmer.h"

#undef min
#undef max

#include <algorithm>
#include <deque>
#include <math.h>
#include <set>

//
// global variables
// 

struct minikmer_node_s *miniKmerGraph = NULL;
uintptr_t MINIKMER_GRAPH_SIZE = 0;

Color COLORMAX;
Color COLORMIN;


//
// statistics functions
//

// Adjacency statistics for mini-kmer graph
// - utilization of finite field
//
// histograms: (in percent)
//   - number of adjacencies
//   - number of adjacencies of any different color
//   - number of adjacent different colors
void printMiniKmerGraphAdjacencyStats(void)
{
  unsigned PARTITIONS = COLORMAX - COLORMIN + 1;
  assert( PARTITIONS == g__PARTITION_COUNT );

  uintptr_t numNodes = 0;
  uintptr_t numCanonNodes = 0;

  unsigned adjHistogram[9];
  bzero(adjHistogram, 9 * sizeof(unsigned));

  unsigned adjDiffColorHistogram[9];
  bzero(adjDiffColorHistogram, 9 * sizeof(unsigned));

  unsigned adjNumDiffColorHistogram[9];
  bzero(adjNumDiffColorHistogram, 9 * sizeof(unsigned));

  unsigned balance[COLORMAX+1];
  bzero(balance, (COLORMAX+1) * sizeof(unsigned));

  for(uintptr_t k=0; k < MINIKMER_GRAPH_SIZE; ++k)
  {
    if( k == canonicalKmer(k, g__MINIKMER_LENGTH) )
      ++ numCanonNodes;
	else
		continue;

    Color thisColor = getMiniKmerColor(k);

    ++ numNodes;
    ++ balance[thisColor];

    std::set<Color> colorSet;

    uintptr_t numAdjacencies = 0;
    uintptr_t numDifferentColor = 0;

    // one direction
    for(int i=0; i < 4; ++i)
    {
      Kmer neighbor = maskKmer((k << 2) | i, g__MINIKMER_LENGTH);
      Kmer canon = canonicalKmer(neighbor, g__MINIKMER_LENGTH);
      Color thatColor = getMiniKmerColor(canon);

	  ++ numAdjacencies;

      if( thisColor != thatColor )
      {
        ++ numDifferentColor;
        colorSet.insert(thatColor);
      }
    }

    // other direction
    for(int i=0; i < 4; ++i)
    {
      Kmer neighbor = maskKmer((k >> 2) | (i << ((g__MINIKMER_LENGTH-1) << 1)), g__MINIKMER_LENGTH);
      Kmer canon = canonicalKmer(neighbor, g__MINIKMER_LENGTH);
      Color thatColor = getMiniKmerColor(canon);

	  ++ numAdjacencies;

      if( thisColor != thatColor )
      {
        ++ numDifferentColor;
        colorSet.insert(thatColor);
      }
    }

    ++ adjHistogram[numAdjacencies];

    ++ adjDiffColorHistogram[numDifferentColor];

    assert( colorSet.size() < 9 );
    ++ adjNumDiffColorHistogram[colorSet.size()];
  }

  printf("MKG Utilization:                      %0.2f%% -- %" PRIuPTR " / %" PRIuPTR " of %" PRIuPTR "\n",
      100.0 * ((float)numNodes / (float)numCanonNodes), numNodes, numCanonNodes, MINIKMER_GRAPH_SIZE);

  // mini-kmer color balance
  {
    uintptr_t idealBalance = numNodes / g__PARTITION_COUNT;
    double maxdiff = 0.0;
    double sum = 0.0;
    for(unsigned i=COLORMIN; i <= COLORMAX; ++i)
    {
      double diff;
      if( balance[i] > idealBalance )
        diff = balance[i] - idealBalance;
      else
        diff = idealBalance - balance[i];

      sum += diff * diff;

      if( diff > maxdiff)
        maxdiff = diff;
    }
    sum /= g__PARTITION_COUNT;
    double stddev = sqrt(sum);
    printf("MKG Color Balance - Std Deviation:");
    printf("  %0.2lf%%\n", 100.0 * stddev / idealBalance);

    printf("MKG Color Balance - Max Deviation:");
    printf("  %0.2lf%%\n", 100.0 * maxdiff / idealBalance);
  }

  printf("MKG Adjacency:                      ");
  {
    float avg = 0;
    for(int i=0; i < 9; ++i)
      avg += i * ((float)adjHistogram[i] / (float)numNodes);
    printf("  %0.2f --", avg);

    for(int i=0; i < 9; ++i)
      printf("  [%d] %5.2f", i, 100.0 * ((float)adjHistogram[i] / (float)numNodes) );
    printf("\n");
  }

  printf("MKG Adjacency - Any Different Color:");
  {
    float avg = 0;
    for(int i=0; i < 9; ++i)
      avg += i * ((float)adjDiffColorHistogram[i] / (float)numNodes);
    printf("  %0.2f --", avg);

    for(int i=0; i < 9; ++i)
      printf("  [%d] %5.2f", i, 100.0 * ((float)adjDiffColorHistogram[i] / (float)numNodes) );
    printf("\n");
  }

  printf("MKG Adjacency - Num Different Color:");
  {
    float avg = 0;
    for(int i=0; i < 9; ++i)
      avg += i * ((float)adjNumDiffColorHistogram[i] / (float)numNodes);
    printf("  %0.2f --", avg);

    for(int i=0; i < 9; ++i)
      printf("  [%d] %5.2f", i, 100.0 * ((float)adjNumDiffColorHistogram[i] / (float)numNodes) );
    printf("\n");
  }
}

