//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// seqNode.h
//
//   sequence graph node: graph-specific node for variable sized sequences
//

// TODO: alignment/packing macros -- NOTE: must respect alignment requirements for tbb::atomic

#ifndef __SEQUENCE_NODE_H
#define __SEQUENCE_NODE_H

#include "sequence.h"

extern SeqNode * GrowSeqNode(SeqGraph * sgraph, SeqNode *node, unsigned bases, bool skip_head_update, bool skip_tail_update);

class SeqNode {  // IMPORTANT: this must be a POD type due to its variable length!
  public:
    union {
        struct {
            counter_t left_count[4];
            counter_t right_count[4];
        };
        struct {
            four_counter_t left_side;
            four_counter_t right_side;
        };
        eight_counter_t connections;
    };

    union {
        struct {
            color_t left_color[4];
            color_t right_color[4];
        };
        struct {
            four_color_t left_side_colors;
            four_color_t right_side_colors;
        };
#ifdef COLOR_8BIT
        #error fixme
#else
        // FIXME: forces colors to be 16-bits
        struct {
            vid_t left_vid;
            vid_t right_vid;
        };
#endif
    };

    SeqNode *head_next;
    SeqNode *tail_next;

    uint8_t flags;
    uint8_t slice_flags;
    threadid_t claim_tid;
    color_t slice_color;

    uint64_t kmer_occurrences;

    Sequence sequence;  // *** IMPORTANT: this must be the _LAST_ member as it is variable length! ***

  public:
    SeqNode(bool) : sequence(false) {}  // used for stack allocation

    SeqNode(const uint16_t alloc_bases) :
		connections(0), left_side_colors(0), right_side_colors(0),
        head_next(NULL), tail_next(NULL), 
        flags(0), slice_flags(0), claim_tid(0), slice_color(0), kmer_occurrences(0), sequence(alloc_bases) {}

    SeqNode(const SeqNode* node, const uint16_t alloc_bases) :
		connections(node->connections), left_side_colors(node->left_side_colors),
        right_side_colors(node->right_side_colors),
        head_next(NULL), tail_next(NULL),
        flags(node->flags), // FIXME: don't copy all flag bits?
        slice_flags(0), claim_tid(0), slice_color(0), kmer_occurrences(node->kmer_occurrences), sequence(node->sequence, alloc_bases) {}

    private:
        // private declarations with no definition to disable compiler generated
        SeqNode();

    public:
        ~SeqNode() { }

    static size_t ComputeNodeAllocatedBytes(uint16_t total_bases)
        { return sizeof(SeqNode) + Sequence::ComputeAllocatedBytes(total_bases); }
    size_t GetNodeAllocatedBytes() const
        { return sizeof(SeqNode) + sequence.GetAllocatedBytes(); }

    size_t GetNodeSerializedBytes() const;

    void loadFromFile(FILE *, int fileFormat);
    size_t emitToFile(FILE *, int fileFormat, uintptr_t nodeIndex=17);

    size_t loadFromMemory(char *memory);
    size_t saveToMemory(char *memory, size_t space_remaining);

    void getNeighborCounts(unsigned *l_count, unsigned *r_count); // const;

    double getNodeKmerCoverage(unsigned kmer_length) {
        uint64_t node_kmers = sequence.get_length() - kmer_length + 1;
        double node_kmer_coverage = ((double)kmer_occurrences) / node_kmers;
        return node_kmer_coverage;
    }

} PACKED;

//
// use this class to stack allocate the variable length SeqNode of maximal length
//
class SeqNode_StackAllocated {  // use instead of 'alloca()' since inlining alloca is dangerous
    public:
        SeqNode_StackAllocated() : __dummy_seqnode(false), __dummy_stacksequence() {}

    private:
        SeqNode __dummy_seqnode;
        Sequence_StackAllocated __dummy_stacksequence;
} PACKED;

//
// helper functions to explicitly touch the slice_flags field
//
static inline void resetNodeSliceFlags(SeqNode *node) {
    // assert( !isNodeDead<SeqNode>(node) );
	node->slice_flags = 0;
}

static inline uint8_t getNodeLeftSliceFlags(SeqNode *node) {
    return static_cast<uint8_t>( (node->slice_flags >> 4) & 0xf );
}

static inline uint8_t getNodeRightSliceFlags(SeqNode *node) {
    return static_cast<uint8_t>( (node->slice_flags & 0xf) );
}

static inline void setNodeLeftSliceFlag(SeqNode *node, unsigned i) {
    assert( i < 4 );
    node->slice_flags |= 0x1 << (4+i);
}

static inline void setNodeRightSliceFlag(SeqNode *node, unsigned i) {
    assert( i < 4 );
    node->slice_flags |= 0x1 << i;
}

static inline bool getNodeLeftSliceFlag(SeqNode *node, unsigned i) {
    assert( i < 4 );
    return static_cast<bool>( (node->slice_flags >> (4+i)) & 0x1 );
}

static inline bool getNodeRightSliceFlag(SeqNode *node, unsigned i) {
    assert( i < 4 );
    return static_cast<bool>( (node->slice_flags >> (i)) & 0x1 );
}

#endif // __SEQUENCE_NODE_H
