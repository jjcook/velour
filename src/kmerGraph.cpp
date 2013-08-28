//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// kmerGraph.cpp
//

#include "types.h"
#include "kmerGraph.h"

KmerGraph::KmerGraph(uintptr_t _buckets) :
    buckets(round_up_power_of_2(_buckets)), hashmask(buckets-1),
    table(static_cast<KmerNode**>(calloc(buckets, sizeof(KmerNode *)))), node_count(0) {}

KmerGraph::~KmerGraph() { free(table); }


