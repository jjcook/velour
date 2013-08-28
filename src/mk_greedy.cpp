//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// mk_greedy.cpp
//
// partitioning method: greedy incremental

#include "types.h"
#include "minikmer.h"

#undef min
#undef max

#include <algorithm>
#include <deque>
#include <math.h>
#include <set>

static Color nextColor;

static uintptr_t colorLimit;
static uintptr_t *colorCount;

//
// version-dependent jazz
//

Color computeKmerColor(Kmer fullKmer)
{
  return getMiniKmerColor( computeMiniKmer(fullKmer, g__FULLKMER_LENGTH, g__MINIKMER_LENGTH) );
}

void initializeMiniKmer(void)
{
  // global variables
  MINIKMER_GRAPH_SIZE = 1UL << (2 * g__MINIKMER_LENGTH);
  miniKmerGraph = (struct minikmer_node_s *) calloc(MINIKMER_GRAPH_SIZE, sizeof(struct minikmer_node_s));

  COLORMIN = 1;
  COLORMAX = COLORMIN + g__PARTITION_COUNT - 1;

  // static variables
  nextColor = COLORMIN;
  assert( (MINIKMER_GRAPH_SIZE / 2) % g__PARTITION_COUNT == 0 ); // ensure colorLimit is good
  colorLimit = (MINIKMER_GRAPH_SIZE / 2) / g__PARTITION_COUNT;
  colorCount = (uintptr_t *) calloc(COLORMIN+g__PARTITION_COUNT, sizeof(unsigned));

  printf("MKG Initialization:  mini-k = %d  full-k = %d  colors = %u  totalNodes = %u\n",
      g__MINIKMER_LENGTH, g__FULLKMER_LENGTH, g__PARTITION_COUNT, MINIKMER_GRAPH_SIZE );
}

void destroyMiniKmer(void)
{
  free(miniKmerGraph);
}


// round-robin
static Color
getNextFreeColor(void)
{
  Color initialNextColor = nextColor;
  while( colorCount[nextColor] >= colorLimit )
  {
    if( nextColor < COLORMAX )
      ++ nextColor;
    else
      nextColor = COLORMIN;

    assert( nextColor != initialNextColor );
  }

  Color c = nextColor++;
  assert( c >= COLORMIN );
  assert( c <= COLORMAX );

  if( nextColor > COLORMAX )
    nextColor = COLORMIN;

  return c;
}

// increase fairness of round-robin?
static void
returnColor(Color oldColor)
{
  if( nextColor == COLORMIN && oldColor == COLORMAX )
    nextColor = COLORMAX;
  else if( (nextColor-1) == oldColor )
    -- nextColor;
}


// foreach kmer in this read: observe to create a partition map
void
observeSequence(char *seq, int kmer_length)
{
    assert( false && "TODO FIX! REVERSED KMER ENDIANNESS" );
  std::deque<Kmer> tipQueue;
  std::deque<std::pair<Color, Color> > adjacencyQueue;

  Kmer kmer = 0, anti_kmer = 0;
  int i;

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

  Kmer rc_kmer = (anti_kmer >> (64 - double_kmer_length));
  bool sense_reversed = rc_kmer < kmer;  
  Kmer canonical_kmer = sense_reversed ? rc_kmer : kmer;

  Kmer miniKmer = computeMiniKmer(canonical_kmer, g__FULLKMER_LENGTH, g__MINIKMER_LENGTH);

  Color currentColor = 0;
  bool isTip = false;

  // handle first mini-kmer
  if( getMiniKmerColor(miniKmer) == 0 )
  {
    currentColor = getNextFreeColor();

    setMiniKmerColor(miniKmer, currentColor);
    ++ colorCount[currentColor];

    isTip = true;
    tipQueue.push_front(miniKmer);
  }
  else
    currentColor = getMiniKmerColor(miniKmer);
    
  assert( getMiniKmerColor(miniKmer) != 0 );
  assert( currentColor != 0 );


  // make each succeeding Kmer
  char c;
  while((c = seq[i]) != 0)
  {
    ++ i;

    // read the next base and extend both the kmer and the anti_kmer
    Nucleotide base = BASE_MAP[(int)c];
    anti_kmer >>= 2;
    kmer <<= 2;
    if( base > 0x3 ) break; // stop at first 'unknown' base
    //assert(base <= 0x3);
    anti_kmer |= ((Kmer)base ^ 0x3) << 62;
    kmer |= base;

    kmer &= mask;

    Kmer rc_kmer = (anti_kmer >> (64 - double_kmer_length));
    bool sense_reversed = rc_kmer < kmer;  
    Kmer canonical_kmer = sense_reversed ? rc_kmer : kmer;

    Kmer miniKmer = computeMiniKmer(canonical_kmer, g__FULLKMER_LENGTH, g__MINIKMER_LENGTH);

    if( getMiniKmerColor(miniKmer) == 0 )
    {
      // limit the size of a color's membership
      if( colorCount[currentColor] >= colorLimit )
      {
        //Color oldColor = currentColor;
        currentColor = getNextFreeColor();

        // forced color transition induces tip and adjacency
        isTip = true;
        //colorAdjacency.set(adjacencyIndex(oldColor,currentColor));
        //adjacencyQueue.push_front(std::make_pair(oldColor,currentColor));
      }

      setMiniKmerColor(miniKmer, currentColor);
      ++ colorCount[currentColor];

      if( isTip )
        tipQueue.push_front(miniKmer);
    }
    else
    {
      Color oldColor = currentColor;
      currentColor = getMiniKmerColor(miniKmer);

      if( currentColor != oldColor )
      {
        returnColor(oldColor); // try to recycle color
        isTip = false; // tip ends at color change

        /*
        // set adjacency, even if it is a tip we might recolor
        colorAdjacency.set(adjacencyIndex(oldColor,currentColor));
        adjacencyQueue.push_front(std::make_pair(oldColor,currentColor));
        */

        // recolor trailing mini-kmers of tip -- starting with the end
        while( !tipQueue.empty() && colorCount[currentColor] < colorLimit )
        {
          Kmer someKmer = tipQueue.back();
          tipQueue.pop_back();

          // de-color the mini-kmer
          //assert( colorCount[getMiniKmerColor(someKmer)] > 0 );
          --colorCount[getMiniKmerColor(someKmer)];
          miniKmerGraph[someKmer].color = 0; // explicitly uncolor mini-kmer

          // re-color
          setMiniKmerColor(someKmer, currentColor);
          ++ colorCount[currentColor];
        }
      }
      else
        tipQueue.clear();  // tip is same color, so just drop entries
    }

    assert( getMiniKmerColor(miniKmer) != 0 );
    assert( currentColor != 0 );

    // TODO ?? arc connection stuff
  }

  /*
  // delete added adjacencies if either color is no longer present
  while( !adjacencyQueue.empty() )
  {
    std::pair<Color, Color> cp = adjacencyQueue.back();
    adjacencyQueue.pop_back();

    Color a = cp.first;
    Color b = cp.second;

    if( colorCount[a] == 0 || colorCount[b] == 0 )
      colorAdjacency.reset(adjacencyIndex(a,b));
  }
  */
}


//
// simple non-helper functions
//

void printMiniKmerStats(void)
{
  printMiniKmerGraphAdjacencyStats();
}

