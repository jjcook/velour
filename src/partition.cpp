//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// partition.cpp
//
//   common framework for all partitioning methods
//

#include "types.h"

#include "minikmer.h"

#include "alignedWriteBuffer.h"

// file descriptors for subsequence buckets
static AlignedWriteBuffer **partitionFiles = NULL;

// partitioner statistics
//
#define HISTBUCKETS 21
static unsigned combinedHistogram[HISTBUCKETS];
static uint64_t numSubsequenceKmers = 0;
static uint64_t numSubsequences = 0;


// forward declarations
static void createSubsequenceFiles(void);
static void closeSubsequenceFiles(void);
// end forward declarations

void initializePartitionerPreObservation(void)
{
    initializeMiniKmer();
}

void initializePartitionerPostObservation(void)
{
    if( !g__PARTITION_NOEMIT )
        createSubsequenceFiles();

    bzero(combinedHistogram, HISTBUCKETS * sizeof(unsigned));
}

void destroyPartitioner(void)
{
  if( !g__PARTITION_NOEMIT )
    closeSubsequenceFiles();
  
  // TODO
  destroyMiniKmer();
}


static void
createSubsequenceFiles(void)
{
    char buf[PATH_MAX+1];

    partitionFiles = new AlignedWriteBuffer* [g__PARTITION_COUNT+1];

    for(unsigned fileNo = 1; fileNo <= g__PARTITION_COUNT; ++ fileNo )
    {
        sprintf(buf, "%s/Subsequences-%d.loom", g__WORK_LOOM_DIRECTORY, fileNo);
        partitionFiles[fileNo] = new AlignedWriteBuffer(buf);
    }
}

static void
closeSubsequenceFiles(void)
{
    for(unsigned fileNo = 1; fileNo <= g__PARTITION_COUNT; ++ fileNo ) {
        delete partitionFiles[fileNo];
    }
    delete [] partitionFiles;
}

static void
emitSubsequence(Color color, Sequence *subsequence, bool hasPrefix, Color prefixColor, bool hasSuffix, Color suffixColor, int kmer_length)
{
  assert( color != 0 );
  assert( color <= g__PARTITION_COUNT );

  // emit subsequence to file
  if( !g__PARTITION_NOEMIT )
  {
    assert( color >= 1 );
    assert( color <= g__PARTITION_COUNT );
    assert( !hasPrefix || prefixColor != color );
    assert( !hasSuffix || suffixColor != color );

    size_t write_size = sizeof(char) + subsequence->GetSerializedBytes() + (hasPrefix ? sizeof(Color) : 0) + (hasSuffix ? sizeof(Color) : 0);
    char * memory = static_cast<char *>( partitionFiles[color]->RequestBuffer(write_size) );
    char * memory_end = memory + write_size;

	unsigned char mode = (hasPrefix ? 1 : 0) + (hasSuffix ? 2 : 0);
    *(memory++) = mode;

    memory += subsequence->Save_Memory(memory, write_size); // FIXME: don't use write_size here

	if (hasPrefix) {
        Color * memory_prefix = reinterpret_cast<Color *>(memory); // FIXME: use memcpy instead
        *memory_prefix = prefixColor;
        memory += sizeof(Color);
	}
	if (hasSuffix) {
        Color * memory_suffix = reinterpret_cast<Color *>(memory); // FIXME: use memcpy instead
        *memory_suffix = suffixColor;
        memory += sizeof(Color);
    }

    /*// DEBUG: print subsequence to stdout
    printf("%01d %01d ", hasPrefix, hasSuffix);
    subsequence->Print(); fputc('\n', stdout);*/

    assert( memory == memory_end );
  }

  // update summary statistics
  unsigned numKmers = subsequence->get_length() - kmer_length - (hasPrefix?1:0) - (hasSuffix?1:0) + 1;
  assert( numKmers > 0 );

  ++ numSubsequences;
  numSubsequenceKmers += numKmers;

  if( numKmers > HISTBUCKETS-1 )
    numKmers = HISTBUCKETS-1;
  ++ combinedHistogram[numKmers];
}


void printPartitionerStatistics(void)
{
  uint64_t totalKmers = numSubsequenceKmers;

  float avgKmers = (float)numSubsequenceKmers / (float)numSubsequences;
  printf("Average full-kmers per subsequence:  %0.2f --", avgKmers);

  unsigned runningSumKmers = 0;
  for(int i=1; i < HISTBUCKETS-1; ++i)
  {
    unsigned bucketKmers = i * combinedHistogram[i];
    runningSumKmers += bucketKmers;
    printf("  [%d] %.1f%%", i, 100.0 * (((float)bucketKmers) / (float)totalKmers) );
  }
  printf("  [%d] %.1f%%", HISTBUCKETS-1, 100.0 * ((float)(numSubsequenceKmers-runningSumKmers) / (float)totalKmers) );
  printf("\n");
}

namespace {
    struct lambda_distribute
    {
        lambda_distribute(Sequence *subseq) : subsequence_(subseq) {}
        void operator()(Kmer k, Kmer canon_k, bool new_sense_reversed, Nucleotide base, Nucleotide last_base, unsigned kmer_length, bool isFirstKmer, bool isLastKmer)
        {
            if (isFirstKmer) {
                subsequence_->InitializeWithKmer_Unsafe(k, kmer_length);
                hasPrefix_ = false;
                Color kmerColor = computeKmerColor(canon_k);
                currentColor_ = kmerColor;
                return;
            }

            subsequence_->AppendBase_Unsafe(base);

            // calculate color for this kmer
            Color kmerColor = computeKmerColor(canon_k);

            // detect and handle color transition
            if( currentColor_ != kmerColor ) {
                bool hasSuffix = true;
                Color suffixColor = kmerColor;
                emitSubsequence(currentColor_, subsequence_, hasPrefix_, prefixColor_, hasSuffix, suffixColor, kmer_length);

                // keep the last two kmers, 'previous' and 'current'
                subsequence_->PopFront(subsequence_->get_length() - (kmer_length + 1));
                assert( subsequence_->get_length() == (kmer_length + 1) );

                hasPrefix_ = true;
                prefixColor_ = currentColor_;

                currentColor_ = kmerColor;
            }

            // handle end-of-read
            if (isLastKmer) {
                bool hasSuffix = false;
                Color suffixColor = 0;
                emitSubsequence(currentColor_, subsequence_, hasPrefix_, prefixColor_, hasSuffix, suffixColor, kmer_length);
            }
            return;
        }
        private:
            Sequence *subsequence_;
            bool hasPrefix_;
            Color prefixColor_;
            Color currentColor_;
    };
}

// foreach kmer in this read: look-up in partition map
//   and emit longest subsequence of same partition (+ prefix/suffix)
void distributeSequence(char *seq, int kmer_length)
{
    Sequence_StackAllocated memory;
    Sequence *stack_seq = new (&memory) Sequence(Sequence::MAX_BASES);
    stack_seq->InitializeWithString(seq);
    
    Sequence_StackAllocated submemory;
    Sequence *stack_sub_seq = new (&submemory) Sequence(Sequence::MAX_BASES);

    lambda_distribute functor(stack_sub_seq);
    sequence_process_kmers(stack_seq, kmer_length, functor);
}

