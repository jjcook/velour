//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

#include "types.h"
#include "kmerNode.h"

void KmerNode::getNeighborCounts(unsigned *l_count, unsigned *r_count) // const
{
    unsigned left_count = 0, right_count = 0;
#ifdef VELOUR_TBB // XXX: only needed for directly constructed version that parallel tip-clips
    union {
        counter_t count[4];
        four_counter_t side;
    } local_left, local_right;
    local_left.side = ATOMIC_LOAD(this->left_side);
    local_right.side = ATOMIC_LOAD(this->right_side);
    for (int i = 0 ; i < 4 ; i ++) {
        left_count += (local_left.count[i] != 0);
        right_count += (local_right.count[i] != 0);
    }
#else
    for (int i = 0 ; i < 4 ; i ++) {
        left_count += (this->left_count[i] != 0);
        right_count += (this->right_count[i] != 0);
    }
#endif
    *l_count = left_count;
    *r_count = right_count;
}


