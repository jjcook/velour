//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// flags.h
//
// kg_node and sg_node flags
//

#ifndef __FLAGS_H
#define __FLAGS_H

#define NODE_MERGE_LEFT   (1 << 0)
#define NODE_MERGE_RIGHT  (1 << 1)
#define NODE_MERGED       (1 << 2)
#define NODE_DEAD         (1 << 3)

// NOTE: below values chosen carefully, do not modify carelessly
#define NODE_WEIGHT_NONE          (0x0 << 4)    // no internal kmers
#define NODE_WEIGHT_SINGLECOPY    (0x1 << 4)    // internal kmers are all single copy
#define NODE_WEIGHT_NONSINGLECOPY (0x2 << 4)    // internal kmers are all non-single copy
#define NODE_WEIGHT_MIXED         (0x3 << 4)    // internal kmers are mixed single + non-single copy
#define NODE_WEIGHT_MASK          (0x3 << 4)    // mask used to access flag
//#define DONOTUSE        (1 << 5)          // used by NODE_WEIGHT flags

template<typename T>
static inline bool
isNodeMerged(T *node) {
    return node->flags & NODE_MERGED;
}

template<typename T>
static inline bool
isNodeMergingLeft(T *node) {
    return node->flags & NODE_MERGE_LEFT;
}

template<typename T>
static inline bool
isNodeMergingRight(T *node) {
    return node->flags & NODE_MERGE_RIGHT;
}

template<typename T>
static inline void
setNodeMerged(T *node) {
    assert( !isNodeMerged<T>(node) );
    node->flags |= NODE_MERGED;
}

template<typename T>
static inline void
setNodeMergingLeft(T *node) {
    node->flags |= NODE_MERGE_LEFT;
}

template<typename T>
static inline void
setNodeMergingRight(T *node) {
    node->flags |= NODE_MERGE_RIGHT;
}

template<typename T>
static inline void
resetNodeMergingFlags(T *node) {
    node->flags &= ~(NODE_MERGE_LEFT | NODE_MERGE_RIGHT | NODE_MERGED);
}

template<typename T>
static inline bool
isNodeDead(T *node) {
    return node->flags & NODE_DEAD;
}

template<typename T>
static inline bool
atomic_isNodeDead(T *node) {
    return ATOMIC_LOAD(node->flags) & NODE_DEAD;
}

template<typename T>
static inline void
setNodeDead(T *node) {
    //assert( !isNodeDead<T>(node) ); // don't check this, as can re-enter in structural loop
    node->flags |= NODE_DEAD;
}

template<typename T>
static inline void
atomic_setNodeDead(T *node) {
    //assert( !isNodeDead<T>(node) ); // don't check this, as can re-enter in structural loop
    while(1) {
        uint8_t old_flags = ATOMIC_LOAD(node->flags);
        if (ATOMIZE(node->flags).compare_and_swap((old_flags | NODE_DEAD), old_flags) == old_flags) break;
    }
}

/*template<typename T>
static inline void
resetNodeFlags(T *node) {
    // assert( !isNodeDead<T>(node) );
	node->flags = 0;
}*/

static inline bool
isNodeWeight_None(SeqNode *node) {
    return ( (node->flags & NODE_WEIGHT_MASK) == NODE_WEIGHT_NONE );
}

static inline bool
isNodeWeight_AllSingleCopy(SeqNode *node) {
    return ( (node->flags & NODE_WEIGHT_MASK) == NODE_WEIGHT_SINGLECOPY );
}

static inline bool
isNodeWeight_HasSingleCopy(SeqNode *node) {
    return ( ((node->flags & NODE_WEIGHT_MASK) == NODE_WEIGHT_SINGLECOPY) ||
             ((node->flags & NODE_WEIGHT_MASK) == NODE_WEIGHT_MIXED)         );
}

static inline bool
isNodeWeight_AllNonSingleCopy(SeqNode *node) {
    return ( (node->flags & NODE_WEIGHT_MASK) == NODE_WEIGHT_NONSINGLECOPY );
}

static inline bool
isNodeWeight_Mixed(SeqNode *node) {
    return ( (node->flags & NODE_WEIGHT_MASK) == NODE_WEIGHT_MIXED );
}

// use when concatenating nodes
static inline void
updateNodeWeight(SeqNode *node, counter_t new_edge_weight) {
    assert( new_edge_weight > 0 );
    if (new_edge_weight > 1) {
        node->flags |= NODE_WEIGHT_NONSINGLECOPY; // can result in either NONSINGLECOPY or MIXED
    } else { // new_edge_weight == 1
        node->flags |= NODE_WEIGHT_SINGLECOPY; // can result in either SINGLECOPY or MIXED
    }
}

#endif // __FLAGS_H
