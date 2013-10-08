//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// split.cpp
//

#include "types.h"
#include "split.h"

SplitBuckets::SplitBuckets(bool skipSelfBucket) :
    selfBucket(NULL), finalBucket(NULL), inboxBucket(new FILE *[g__PARTITION_COUNT+1]),
#ifdef VELOUR_TBB
    inboxMutex(new tbb::queuing_mutex [g__PARTITION_COUNT+1]),
#endif // VELOUR_TBB
    stat_selfNodes(0), stat_finalNodes(0), stat_inboxNodes(0),
    stat_maxRedistFinalComponentSize(0), stat_maxRedistInboxComponentSize(0)
{
    // TODO: choose size more safely (allocator doesn't know about this) and more optimally
    size_t FILE_BUFFER_SIZE = (g__MEMORY_FOOTPRINT_LIMIT >> 4) / g__PARTITION_COUNT;
    if (FILE_BUFFER_SIZE > 2*1024*1024UL) {
        FILE_BUFFER_SIZE = 2*1024*1024UL;
    }

    // null the inbox file pointers
    for (unsigned i=0; i <= g__PARTITION_COUNT ; ++i) {
        inboxBucket[i] = NULL;
    }

    // self bucket: truncate
    if (!skipSelfBucket)
    {
        char selfbucket_filename[PATH_MAX+1];
        sprintf(selfbucket_filename, "%s/inbox-for-%u/SelfBucket-%u.bucket", g__WORK_INBOX_ROOT_DIRECTORY, g__PARTITION_INDEX, g__PARTITION_INDEX);
        int filedes = open(selfbucket_filename, O_CREAT | O_TRUNC | O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
        if (filedes == -1) {
            fprintf(stderr, "FAILED: open() of %s\n", selfbucket_filename);
            perror(NULL);
            exit(1);
        }
        FILE *fp = fdopen(filedes, "w");
        if (fp == NULL) {
            fprintf(stderr, "FAILED: fdopen() of %s\n", selfbucket_filename);
            perror(NULL);
            exit(1);
        }
        selfBucket = fp;

        // configure buffer for writes
        char * buffer = static_cast<char*>( malloc( FILE_BUFFER_SIZE ) );
        if (buffer != NULL) {
            setbuffer(fp, buffer, FILE_BUFFER_SIZE);
        }
    }

    // final bucket: truncate
    {
        char finalbucket_filename[PATH_MAX+1];
        sprintf(finalbucket_filename, "%s/FinalBucket-from-%u.bucket", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX);
        int filedes = open(finalbucket_filename, O_CREAT | O_TRUNC | O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
        if (filedes == -1) {
            fprintf(stderr, "FAILED: open() of %s\n", finalbucket_filename);
            perror(NULL);
            exit(1);
        }
        FILE *fp = fdopen(filedes, "w");
        if (fp == NULL) {
            fprintf(stderr, "FAILED: fdopen() of %s\n", finalbucket_filename);
            perror(NULL);
            exit(1);
        }
        finalBucket = fp;

        // configure buffer for writes
        char * buffer = static_cast<char*>( malloc( FILE_BUFFER_SIZE ) );
        if (buffer != NULL) {
            setbuffer(fp, buffer, FILE_BUFFER_SIZE);
        }
    }

    // inbox buckets: truncate
    {
        for (unsigned i=g__PARTITION_INDEX+1; i <= g__PARTITION_COUNT ; ++i) {
            char inboxbucket_filename[PATH_MAX+1];
            sprintf(inboxbucket_filename, "%s/inbox-for-%u/InboxBucket-from-%u.bucket", g__WORK_INBOX_ROOT_DIRECTORY, i, g__PARTITION_INDEX);
            int filedes = open(inboxbucket_filename, O_CREAT | O_TRUNC | O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
            if (filedes == -1) {
                fprintf(stderr, "FAILED: open() of %s\n", inboxbucket_filename);
                perror(NULL);
                exit(EXIT_FAILURE);
            }
            FILE *fp = fdopen(filedes, "w");
            if (fp == NULL) {
                fprintf(stderr, "FAILED: fdopen() of %s\n", inboxbucket_filename);
                perror(NULL);
                exit(EXIT_FAILURE);
            }
            inboxBucket[i] = fp;

            // configure buffer for writes
            char * buffer = static_cast<char*>( malloc( FILE_BUFFER_SIZE ) );
            if (buffer != NULL) {
                setbuffer(fp, buffer, FILE_BUFFER_SIZE);
            }
        }
    }
}

SplitBuckets::~SplitBuckets()
{
    // close files
    if (selfBucket != NULL) { fclose(selfBucket); }
    if (finalBucket != NULL) { fclose(finalBucket); }
    for (unsigned i=g__PARTITION_INDEX+1; i <= g__PARTITION_COUNT ; ++i) {
        if (inboxBucket[i] != NULL) {
            fclose(inboxBucket[i]);
        }
    }

    // deallocations
    delete [] inboxBucket;
#ifdef VELOUR_TBB
    delete [] inboxMutex;
#endif // VELOUR_TBB
}

namespace {

struct lambda_split {
    lambda_split(SeqGraph *graph, SplitBuckets *sb) : graph(graph), sb(sb) {}

    void operator()(SeqNode *root) {
        SerialComponent<SeqGraph, SeqNode> component(graph, root); // forms the component

        FILE *targetBucket = NULL;
        uintptr_t *targetBucketCounter = NULL;
        uintptr_t *targetMaxRedistComponentSize = NULL;

#ifdef VELOUR_TBB
        tbb::queuing_mutex *targetMutex;
#endif

        assert( component.get_min_frontier_index() != g__PARTITION_INDEX );
        //assert( component.maxFrontierIndex != FINALBUCKET || component.get_min_frontier_index() == FINALBUCKET );
        if (component.get_min_frontier_index() < g__PARTITION_INDEX) {
            if (sb->selfBucket == NULL) return;
            //sb->stat_maxRedistInboxComponentSize = max(component.get_node_count(), sb->stat_maxRedistInboxComponentSize);
            targetBucket = sb->selfBucket;
            targetBucketCounter = &sb->stat_selfNodes;
#ifdef VELOUR_TBB
            targetMutex = &sb->selfMutex;
#endif
        } else if (component.get_min_frontier_index() == (g__PARTITION_COUNT+1)) {
            targetMaxRedistComponentSize = &sb->stat_maxRedistFinalComponentSize;
            targetBucket = sb->finalBucket;
            targetBucketCounter = &sb->stat_finalNodes;
#ifdef VELOUR_TBB
            targetMutex = &sb->finalMutex;
#endif
        } else {
            targetMaxRedistComponentSize = &sb->stat_maxRedistInboxComponentSize;
            assert( component.get_min_frontier_index() > g__PARTITION_INDEX );
            assert( component.get_min_frontier_index() < (g__PARTITION_COUNT+1));
            targetBucket = sb->inboxBucket[component.get_min_frontier_index()];
            targetBucketCounter = &sb->stat_inboxNodes;
#ifdef VELOUR_TBB
            targetMutex = &sb->inboxMutex[component.get_min_frontier_index()];
#endif
        }

        uintptr_t numNodes = component.nodes_.size();

        // compute the serialized size of the component entry
        size_t serialized_bytes = sizeof(serialized_bytes) + sizeof(numNodes);
        for (std::deque<SeqNode*>::iterator it = component.nodes_.begin(); it != component.nodes_.end(); ++it) {
            SeqNode *node = *it;
            serialized_bytes += node->GetNodeSerializedBytes();
        }

#ifdef VELOUR_TBB
        tbb::queuing_mutex::scoped_lock lock( *targetMutex );
#endif

        // stat: max component size for final or inbox -- TBB safe since after lock
        if (targetMaxRedistComponentSize != NULL && component.get_node_count() > *targetMaxRedistComponentSize) {
            *targetMaxRedistComponentSize = component.get_node_count();
        }

        if (fwrite(&serialized_bytes, sizeof(serialized_bytes), 1, targetBucket) != 1) {
            fprintf(stderr, "ERROR: failed to fwrite()\n");
            exit(EXIT_FAILURE);
        }

        if (fwrite(&numNodes, sizeof(numNodes), 1, targetBucket) != 1) {
            fprintf(stderr, "ERROR: fwrite() to bucket in lambda_split.\n");
            exit(EXIT_FAILURE);
        }

        size_t component_bytes = sizeof(serialized_bytes) + sizeof(numNodes);

        while(!component.nodes_.empty()) {
            SeqNode *node = component.nodes_.back();
            component.nodes_.pop_back();

            component_bytes += node->emitToFile(targetBucket, BUCKET);
            ATOMIC_ADD(*targetBucketCounter,1);

            //if (selfBucket == NULL) { // delete emitted nodes when incremental emittal
            setNodeDead<SeqNode>(node); // causes node to be removed from graph and deallocated by the graph iterator
        }
        assert( serialized_bytes == component_bytes );
    }
  private:
    SeqGraph *graph;
    SplitBuckets *sb;
};

} // namespace: anonymous

void SplitBuckets::split(SeqGraph *graph)
{
    graph->resetFlags();

    // form and emit components
    lambda_split s(graph, this);
    sg_for_each_mutable(graph, s);

    // FIXME: fold this into above
    struct lambda {
        static void id(SeqNode *) {}
    };
    sg_for_each_mutable(graph, lambda::id); // make sure things are deleted
}

void SplitBuckets::split_nodelist(SeqGraph *graph, flow_nodelist_t *nodelist)
{
    for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
        resetNodeMergingFlags<SeqNode>(*it);
    }

    // form and emit components
    lambda_split s(graph, this);
    for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
        SeqNode *trav = *it;
        if (!isNodeDead<SeqNode>(trav)) {
            s(trav);
            bool remove_now = isNodeDead<SeqNode>(trav);
            if (remove_now) {
                graph->removeNode(trav); // remove node
    #ifdef VELOUR_TBB
                tls_seqnode_allocator->DeallocateNodeMemory(trav);
    #else
                g__SEQNODE_ALLOCATOR->DeallocateNodeMemory(trav);
    #endif // VELOUR_TBB
            }
        }
    }
}

void SplitBuckets::printStatistics(void)
{
    printf("%"PRIuPTR" nodes split into self.\n", stat_selfNodes);
    printf("%"PRIuPTR" nodes split into final.\n", stat_finalNodes);
    printf("%"PRIuPTR" nodes split into inbox.\n", stat_inboxNodes);

    printf("%"PRIuPTR" max redistributed into final component size.\n", stat_maxRedistFinalComponentSize);
    printf("%"PRIuPTR" max redistributed into inbox component size.\n", stat_maxRedistInboxComponentSize);
    printf("%"PRIuPTR" max redistributed component size.\n", max(stat_maxRedistFinalComponentSize, stat_maxRedistInboxComponentSize));
}

void SplitBuckets::resetStatistics(void)
{
    stat_selfNodes = 0;
    stat_finalNodes = 0;
    stat_inboxNodes = 0;

    stat_maxRedistFinalComponentSize = 0;
    stat_maxRedistInboxComponentSize = 0;
}

