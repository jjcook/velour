//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// kmerNode.h
//
//   kmer graph node: graph-specific implementation
//

// TODO: alignment/packing macros -- NOTE: must respect alignment requirements for tbb::atomic

#ifndef __KMER_NODE_H
#define __KMER_NODE_H

class KmerNode {
  public:
    Kmer kmer;

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
    };

    union {
        threadid_t claim_tid;
        color_t slice_color; // XXX unused: present for templatized serial component
        uint16_t __initializer; // should be as large as the largest union member
    };

    uint8_t flags;

    uint64_t kmer_occurrences;

    KmerNode *next;

  public:
    KmerNode(const Kmer k) :
		kmer(k), connections(0),
        left_side_colors(0),
        right_side_colors(0),
        __initializer(0), flags(0), kmer_occurrences(1), next(NULL) {}

    void getNeighborCounts(unsigned *l_count, unsigned *r_count); // const;
} PACKED;

#endif // __KMER_NODE_H
