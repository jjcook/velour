//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// distance.cpp
//
// k-mer space distance calculation using dynamic programming
//

// FIXME: getNucleotide() changed which end of kmer is zero'th base, maybe broke something here

#include "types.h"
#include "distance.h"

#include "kmer.h"

/*
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>


static void print_row(int* row, int row_length, unsigned is, Kmer source, Kmer target)
{
  printf("S(%c) ", CHAR_BASE_MAP[getNucleotide(source, is)]);
  for(int i=0; i < row_length; ++i)
    printf("[%d:%c]=%d ", i, CHAR_BASE_MAP[getNucleotide(target, i)], row[i]);
  printf("\n");
}
*/

unsigned getDistance(Kmer source, Kmer target, unsigned kmer_length)
{
  unsigned minDistance = kmer_length;
  unsigned dpArray[2][kmer_length];

  bzero(dpArray, kmer_length * 2 * sizeof(unsigned));

/*
#ifdef VERBOSE
  printf("(dist) source k-mer: "); print_kmer(source, kmer_length); printf("\n");
  printf("(dist) target k-mer: "); print_kmer(target, kmer_length); printf("\n");
#endif // VERBOSE
*/

  for(unsigned is=0; is < kmer_length; ++is)
  {
    unsigned *previousRow = &dpArray[((is+1) % 2)][0];
    unsigned *currentRow = &dpArray[is % 2][0];
    bzero(currentRow, kmer_length * sizeof(int));

    for(unsigned it=0; it < kmer_length; ++it)
    {
      if( getNucleotide(target, it) == getNucleotide(source, is) )
      {
/*
#ifdef VERBOSE
        printf("(dist) is=%d isN=%c  it=%d itN=%c\n", is, CHAR_BASE_MAP[getNucleotide(source,is)], it, CHAR_BASE_MAP[getNucleotide(target,it)]);
#endif // VERBOSE
*/
        unsigned previousCount = (it == 0 ? 0 : previousRow[it-1]);
        unsigned currentCount = previousCount + 1;
        currentRow[it] = currentCount;

        // calculate distance
        unsigned ms = min(is, kmer_length - is - currentCount);
        unsigned mt = (ms < is ? it : kmer_length - it - currentCount);
        unsigned distance = ms + (kmer_length - currentCount) + mt;

        if( distance < minDistance )
        {
/*
#ifdef VERBOSE
          printf("previous: "); print_row(previousRow, 21, is-1, source, target);
          printf(" current: "); print_row(currentRow, 21, is, source, target);
          printf("(dist) is=%d  it=%d  ms=%d  mt=%d  len=%d  distance=%d\n", is, it, ms, mt, currentCount, distance);
#endif // VERBOSE
*/
          minDistance = distance;
        }
      }
    }
  }
/*
#ifdef VERBOSE
  printf("(dist) minDistance = %d\n", minDistance);
#endif // VERBOSE
*/
  return minDistance;
}

