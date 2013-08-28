//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// mk_random.cpp
//
// partitioning method: random mini-kmer assignment
//

#include "types.h"
#include "minikmer.h"

#include <math.h>
#include <set>

//
// version-dependent jazz
//

Color computeKmerColor(Kmer fullKmer)
{
  return getMiniKmerColor( computeMiniKmer(fullKmer, g__FULLKMER_LENGTH, g__MINIKMER_LENGTH) );
}

void initializeMiniKmer(void)
{
  MINIKMER_GRAPH_SIZE = 1UL << (2 * g__MINIKMER_LENGTH);
  miniKmerGraph = (struct minikmer_node_s *) calloc(MINIKMER_GRAPH_SIZE, sizeof(struct minikmer_node_s));

  COLORMIN = 1;
  COLORMAX = COLORMIN + g__PARTITION_COUNT - 1;

  for(unsigned k=0; k < MINIKMER_GRAPH_SIZE; ++k)
  {
    if( getMiniKmerColor(k) == 0 && k < reverseComplement(k, g__MINIKMER_LENGTH) )
      setMiniKmerColor(k, COLORMIN + (random() % g__PARTITION_COUNT));
  }

  printf("MKG Initialization:  mini-k = %d  full-k = %d  colors = %u  totalNodes = %u\n",
      g__MINIKMER_LENGTH, g__FULLKMER_LENGTH, g__PARTITION_COUNT, MINIKMER_GRAPH_SIZE );
}

void destroyMiniKmer(void)
{
  free(miniKmerGraph);
}


void printMiniKmerStats(void)
{
  printMiniKmerGraphAdjacencyStats();
}

