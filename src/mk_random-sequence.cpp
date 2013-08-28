//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// mk_random-sequence.cpp
//
// partitioning method: random mini-kmer sequence

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

  printf("MKG Initialization:  mini-k = %d  full-k = %d  colors = %u  totalNodes = %u\n",
      g__MINIKMER_LENGTH, g__FULLKMER_LENGTH, g__PARTITION_COUNT, MINIKMER_GRAPH_SIZE );
}

void destroyMiniKmer(void)
{
  free(miniKmerGraph);
}

// foreach kmer in this read: observe to create a partition map
void
observeSequence(char *seq, int kmer_length) {
    assert( false && "TODO FIX! REVERSED KMER ENDIANNESS" );
  Kmer kmer = 0, anti_kmer = 0;
  int i;

  Color thisColor = (COLORMIN + (random() % (COLORMAX-COLORMIN+1)));

  // read first Kmer
  for (i = 0 ; i < kmer_length ; ++ i) {
    char c = seq[i];
    if (c == 0) { return; }  // end of string without a full Kmer
    Nucleotide base = BASE_MAP[(int) c];
    anti_kmer >>= 2;
    kmer <<= 2;
    assert(base <= 0x3);
    anti_kmer |= ((Kmer)base ^ 0x3) << 62;
    kmer |= base;
  }

  int double_kmer_length = kmer_length << 1;
  Kmer mask = (((Kmer)1) << double_kmer_length) - 1;  // a mask with 2*kmer_length bits

  // lookup first Kmer to get first node
  Kmer rc_kmer = (anti_kmer >> (64 - double_kmer_length));
  bool sense_reversed = rc_kmer < kmer;
  Kmer canonical_kmer = sense_reversed ? rc_kmer : kmer;

  // observe
  Kmer miniKmer = computeMiniKmer(canonical_kmer, g__FULLKMER_LENGTH, g__MINIKMER_LENGTH);
  setMiniKmerColor(miniKmer, thisColor);

  // make each succeeding Kmer
  char c;
  while ((c = seq[i]) != 0) {
    ++ i;

    // read the next base and extend both the kmer and the anti_kmer
    Nucleotide base = BASE_MAP[(int)c];
    anti_kmer >>= 2;
    kmer <<= 2;
    if( base > 0x3 ) break; // stop at first 'unknown' base
    //assert(base <= 0x3);
    anti_kmer |= ((Kmer)base ^ 0x3) << 62;
    kmer |= base;

    // which base just fell off the kmer?
    //Nucleotide last_base = (kmer >> double_kmer_length) & 0x3;
    kmer &= mask;
    //bool new_sense_reversed = false;

    Kmer new_rc_kmer = (anti_kmer >> (64 - double_kmer_length));
    bool new_sense_reversed_check = new_rc_kmer < (kmer & mask);
    canonical_kmer = new_sense_reversed_check ? new_rc_kmer : (kmer & mask);

    // observe
    Kmer miniKmer = computeMiniKmer(canonical_kmer, g__FULLKMER_LENGTH, g__MINIKMER_LENGTH);
    setMiniKmerColor(miniKmer, thisColor);
  }
}

void printMiniKmerStats(void)
{
  printMiniKmerGraphAdjacencyStats();
}

