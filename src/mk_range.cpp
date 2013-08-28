//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// mk_range.cpp
//
// partitioning method: range
//

#include "types.h"
#include "minikmer.h"

//#include <math.h>
//#include <set>

//
// version-dependent jazz
//

Color computeKmerColor(Kmer fullKmer)
{
	Kmer miniKmer = computeMiniKmer(fullKmer, g__FULLKMER_LENGTH, g__MINIKMER_LENGTH);

    // TODO: optimize using log2 etc s.t. don't need division
	uintptr_t rangeSize = 1UL << (2 * g__MINIKMER_LENGTH); // mini-kmer space size
	uintptr_t subrangeSize = rangeSize / g__PARTITION_COUNT;
	return 1 + ( miniKmer / subrangeSize );
}

void initializeMiniKmer(void)
{
}

void destroyMiniKmer(void)
{
}

void printMiniKmerStats(void)
{
  //printMiniKmerGraphAdjacencyStats(); // no graph
}

