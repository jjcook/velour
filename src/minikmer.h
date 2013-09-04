//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

#ifndef _MINIKMER_H_
#define _MINIKMER_H_

//typedef uint16_t Color;

extern Color COLORMAX;
extern Color COLORMIN;

//
// mini-kmer finite field
//

struct minikmer_node_s {
  Color color;
} PACKED;


/*
struct minikmer_node_s {
  Color marked : 1;
  Color color  : 15;
} __attribute__((packed));
*/

/*
struct minikmer_node_s {
  Color marked : 1;
  Color color  : 15;
  uint8_t leftA : 1;
  uint8_t leftC : 1;
  uint8_t leftG : 1;
  uint8_t leftT : 1;
  uint8_t rightA : 1;
  uint8_t rightC : 1;
  uint8_t rightG : 1;
  uint8_t rightT : 1;
} __attribute__((packed));
*/

extern struct minikmer_node_s *miniKmerGraph;
extern uintptr_t MINIKMER_GRAPH_SIZE;

//
// inline functions
//

// extract center mini k-mer from full kmer
static inline Kmer computeMiniKmer(Kmer fullKmer, unsigned fullKmer_length, unsigned miniKmer_length)
{
  assert( fullKmer_length >= miniKmer_length );
  int shift_amount = fullKmer_length - miniKmer_length; // * 2 / 2 -- half difference in bases
  assert( shift_amount % 2 == 0 );
  Kmer miniKmer = maskKmer(fullKmer >> shift_amount, miniKmer_length);
  return canonicalKmer( miniKmer, miniKmer_length );
}

static inline Color getMiniKmerColor(Kmer miniKmer)
{
  uint64_t miniKmerInt = convertKmerToUint64(miniKmer);
  assert( miniKmerGraph != NULL );
  assert( miniKmerInt < MINIKMER_GRAPH_SIZE );
  return miniKmerGraph[miniKmerInt].color + (COLORMIN == 0 ? 1 : 0); // jjcook FIXME: hack!
}

static inline void setMiniKmerColor(Kmer miniKmer, Color newColor)
{
  uint64_t miniKmerInt = convertKmerToUint64(miniKmer);
  assert( miniKmerGraph != NULL );
  assert( miniKmerInt < MINIKMER_GRAPH_SIZE );
  //assert( getMiniKmerColor(miniKmer) == 0 );
  assert( newColor >= COLORMIN );
  assert( newColor <= COLORMAX );
  miniKmerGraph[miniKmerInt].color = newColor;
}

static inline bool isMiniKmerColored(Kmer miniKmer) { return (getMiniKmerColor(miniKmer) != 0); }
static inline bool isMiniKmerUncolored(Kmer miniKmer) { return !isMiniKmerColored(miniKmer); }

/*
// marks

static inline void clearAllMarks(void)
{
  assert( miniKmerGraph != NULL );
  for(unsigned k=0; k < MINIKMER_GRAPH_SIZE; ++k)
    miniKmerGraph[k].marked = 0;
}

static inline bool isKmerMarked(Kmer miniKmer)
{
  assert( miniKmerGraph != NULL );
  assert( miniKmer < MINIKMER_GRAPH_SIZE );
  return (miniKmerGraph[miniKmer].marked == 1);
}

static inline bool isKmerUnmarked(Kmer miniKmer) { return !isKmerMarked(miniKmer); }

static inline void clearKmerMark(Kmer miniKmer)
{
  assert( miniKmerGraph != NULL );
  assert( miniKmer < MINIKMER_GRAPH_SIZE );
  miniKmerGraph[miniKmer].marked = 0;
}

static inline void setKmerMark(Kmer miniKmer)
{
  assert( miniKmerGraph != NULL );
  assert( miniKmer < MINIKMER_GRAPH_SIZE );
  assert( miniKmerGraph[miniKmer].marked == 0 );
  miniKmerGraph[miniKmer].marked = 1;
}
*/

//
// non-inline functions
//

void initializeMiniKmer(void);
void destroyMiniKmer(void);

Color computeKmerColor(Kmer fullKmer);

#ifdef PART_TRAINS_ON_INPUT
void observeSequence(char *seq, int fullKmer_length);
void completeObservation(void);
#endif // PART_TRAINS_ON_INPUT

void printMiniKmerStats(void);
void printMiniKmerGraphAdjacencyStats(void);

#endif /* _MINIKMER_H_ */
