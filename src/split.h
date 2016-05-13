//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// split.h
//
//   distribute (sub-)components to respective buckets
//

#ifndef __SPLIT_H
#define __SPLIT_H

struct SplitBuckets
{
    SplitBuckets(bool skipSelfBucket);
    ~SplitBuckets();

    void split(SeqGraph *);
    void split_nodelist(SeqGraph *graph, flow_nodelist_t *nodelist);

    void printStatistics();
    void resetStatistics();

  //private:
    void openInboxBucketFile(unsigned target);
    void closeInboxBucketFile(unsigned target);

    FILE *selfBucket;
    FILE *finalBucket;
    FILE **inboxBucket;

#ifdef VELOUR_MPI
    off_t *inboxBucketProgress;
#endif // VELOUR_MPI

#ifdef VELOUR_TBB
    tbb::queuing_mutex selfMutex;
    tbb::queuing_mutex finalMutex;
    tbb::queuing_mutex *inboxMutex;
#endif

    uintptr_t stat_selfNodes;
    uintptr_t stat_finalNodes;
    uintptr_t stat_inboxNodes;

    uintptr_t stat_maxRedistFinalComponentSize;
    uintptr_t stat_maxRedistInboxComponentSize;

    size_t FILE_BUFFER_SIZE;
};


#endif // __SPLIT_H
