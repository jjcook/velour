//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

// 
// pg_tipClipping.cpp
//

#include "types.h"

// TODO: configurable parameter?
#define TIP_REMOVAL_LENGTH (2 * g__FULLKMER_LENGTH)


static counter_t g_minimum_edge_weight = 2; // TODO: configurable parameter?


//********************************************************************************
//******************************  Tip Clipping  **********************************
//********************************************************************************


//static unsigned EXTENSIONS = 0;
static uintptr_t TIP_DETACHED = 0;
//static unsigned TIP_DEALLOCATED = 0;
static uintptr_t NODES_DEALLOCATED = 0;
static uintptr_t DISCONNECTED_CHUNKS = 0;
static uintptr_t BOUNDARY_CUTS = 0;

static int REMOVE_NOW_VALUE = 2;
static int REMOVE_NEXT_VALUE = 2; // 4;

static void 
flip_to_remove_value() {
  int temp = REMOVE_NOW_VALUE;
  REMOVE_NOW_VALUE = REMOVE_NEXT_VALUE;
  REMOVE_NEXT_VALUE = temp;
}

// we know that this "node" is a node that has either no neighbors on
// one side -or- connects to a chain of nodes that leads to a single node.
//#ifdef SMALL_NODES
static bool
check_for_tip_removal(KmerNode *node, bool moving_right, bool singlecopy_tip, unsigned length) {
  length ++;
  if (length >= TIP_REMOVAL_LENGTH && !singlecopy_tip) {
    return false;  // can't remove a tip that is too long
  }

  // find out if this tip can be extended.  This has two conditions:
  // CONDITION 1: there is a single neighbor in the extension direction
/*#ifdef VELOUR_TBB
  // Normally we'd leave the data in-place and do this in valid_single_successor
  // Here, we need to load the data exactly once, do the math on local data
  union {
    struct {
      counter_t left_count[4];
      counter_t right_count[4];
    };
    unsigned long long connections;
  } local_counters, local_counters_back;
  local_counters.connections = axtomics::fetch_and_nop<unsigned long long>(&node->connections);
  //  local_counters.connections = node->connections; // single load instruction
  counter_t *counters = moving_right ? local_counters.right_count : local_counters.left_count; 
#else*/
  counter_t   *counters = moving_right ? node->right_count : node->left_count;
//#endif
  int valid_dir = valid_single_successor(counters);

  if (valid_dir == MULTIPLE_SUCCESSORS) {  
    // printf("multiple outgrowths: length %d, weight %d\n", length, 
    //        get_counter_sum(!moving_right ? node->right_count : node->left_count));
    return false; // can't remove a tip with two outgrowths (the base of a Y)
  }

  if (valid_dir == NO_SUCCESSORS) {  
    node->connections = 0;
    node->left_side_colors = 0;
    node->right_side_colors = 0;
/*#ifdef VELOUR_TBB
    axtomics::fetch_and_add<uintptr_t>(&DISCONNECTED_CHUNKS, 1);
#else*/
    ++ DISCONNECTED_CHUNKS;
//#endif
    // printf("a disconnected chunk of length %d\n", length);
    return true;
  }

  // a single successor
  // CONDITION 2: the node in the extension direction has a single
  // neighbor in the direction facing this node
  assert((valid_dir >= 0) && (valid_dir < 4));
  bool sense_changed = counters[valid_dir] < 0;

  color_t *colors = moving_right ? node->right_color : node->left_color;
  KmerNode *next = NULL;
  if (colors[valid_dir] == 0) {
      next = g__KG_HASHTABLE->findNextNode(node, valid_dir, moving_right);
  } else {
      return false;
  }
  assert( next != NULL );

  bool next_moving_right = moving_right ^ sense_changed;
/*#ifdef VELOUR_TBB
  //  local_counters.connections = next->connections;
  local_counters_back.connections = axtomics::fetch_and_nop<unsigned long long>(&(next->connections));
  counter_t *next_back_counters = next_moving_right ? 
    local_counters_back.left_count : local_counters_back.right_count;
#else*/
  counter_t *next_back_counters = next_moving_right ? next->left_count : next->right_count;
//#endif
  int next_back_valid = valid_single_successor(next_back_counters);

  if (next_back_valid != MULTIPLE_SUCCESSORS) {
    bool forceclip = (abs(next_back_counters[next_back_valid]) == CLIP_SINGLECOPY_COUNTER_VALUE);
    singlecopy_tip &= (abs(counters[valid_dir]) == 1) || (abs(counters[valid_dir]) == CLIP_SINGLECOPY_COUNTER_VALUE);
    bool remove = check_for_tip_removal(next, next_moving_right, singlecopy_tip, length);
    if (remove) {
      node->connections = 0;
      node->left_side_colors = 0;
      node->right_side_colors = 0;
    } else if (singlecopy_tip || forceclip) {
      node->connections = 0;
      node->left_side_colors = 0;
      node->right_side_colors = 0;
      next_back_counters[next_back_valid] = 0;
      (next_moving_right?next->left_color:next->right_color)[next_back_valid] = 0;
    }
    return remove || singlecopy_tip || forceclip;
  }

  // can't extend.  Need to determine if we should trim.
  counter_t max_value = get_counter_max(next_back_counters);
  counter_t value = abs(counters[valid_dir]);
  if (value == CLIP_SINGLECOPY_COUNTER_VALUE) value = 1;
  assert(value <= MAX_COUNTER_VALUE);
  if ((value >= g_minimum_edge_weight) && (value == max_value)) {
	 // printf("dominant connection: length %d, weight %d\n", length, abs(counters[valid_dir]));
	 return false;   // can't remove the dominant connection
  }

  // check for multiple single-copy predecessors in next node; mark those arcs as potential clips
  if (value == 1 && value == max_value) {
    for (int i = 0 ; i < 4 ; ++ i) {
        if (next_back_counters[i] == 1) {
            next_back_counters[i] = CLIP_SINGLECOPY_COUNTER_VALUE;
        } else if (next_back_counters[i] == -1) {
            next_back_counters[i] = - CLIP_SINGLECOPY_COUNTER_VALUE;
        }
    }
  }

  // if we are clipping the tip, then mark the current "node" for
  // removal and disconnect the "next" node from the tip.  

#ifdef VERIFY
  verify_node(node, g__KG_HASHTABLE, g__FULLKMER_LENGTH);     // FINE, BEFORE WE LET IT GET INCONSISTENT.
  verify_node(next, g__KG_HASHTABLE, g__FULLKMER_LENGTH);
#endif 

  // Note: we remove the link from "next", but not from "node".  The
  // reason why is that we don't want code to start on the other side
  // thinking that it is a tip that might need removal.  This means 
  // that from "node's" perspective, the graph looks inconsistent, but
  // that is okay, since we're going to remove "node" shortly.
  node->connections = 0;
  // printf("removal candidate: length %d, weight %d\n", length, abs(counters[valid_dir]));
  next_back_valid = (moving_right ? KMER_GET_HEAD_BASE(node->kmer, g__FULLKMER_LENGTH) : KMER_GET_TAIL_BASE(node->kmer, g__FULLKMER_LENGTH));
  if (sense_changed) { next_back_valid ^= 0x3; } // FIXME: i.e. complement??
      
  node->left_side_colors = 0;
  node->right_side_colors = 0;
  (next_moving_right?next->left_color:next->right_color)[next_back_valid] = 0;
  
  // This write is local if using TBB...
/*#ifdef VELOUR_TBB
  // ...write back to the original
  next_back_counters = next_moving_right ? next->left_count : next->right_count;
#endif*/
  next_back_counters[next_back_valid] = 0;

/*#ifdef VELOUR_TBB
  axtomics::fetch_and_add<uintptr_t>(&TIP_DETACHED, 1);
#else*/
  ++ TIP_DETACHED;
//#endif

#ifdef VERIFY
  // verify_node(node, g__FULLKMER_LENGTH);  // THIS VALIDATION WOULD FAIL, SINCE WE DON'T MODIFY NODE
  verify_node(next, g__KG_HASHTABLE, g__FULLKMER_LENGTH);
#endif 
  return true;
}

/*#else // LARGE NODES

// Note: parallelization not performed for large prenodes - this function is NOT threadsafe
bool
check_for_tip_removal(kg_node_t *node, bool moving_right, unsigned length) {
  length ++;
  if (length >= TIP_REMOVAL_LENGTH) {
    // printf("tip too long: length %d, weight %d\n", length, 
    // 		  get_counter_sum(!moving_right ? node->right_count : node->left_count));
    return false;  // can't remove a tip that is too long
  }

  // find out if this tip can be extended.  This has two conditions:
  // CONDITION 1: there is a single neighbor in the extension direction
  kg_node_t **pointers = moving_right ? node->right : node->left;
  counter_t   *counters = moving_right ? node->right_count : node->left_count;
  int valid_dir = valid_single_successor(counters);
  if (valid_dir == MULTIPLE_SUCCESSORS) {  
	 // printf("multiple outgrowths: length %d, weight %d\n", length, 
	 //        get_counter_sum(!moving_right ? node->right_count : node->left_count));
	 return false; // can't remove a tip with two outgrowths (the base of a Y)
  }
  kg_node_t *next = NULL;
  counter_t *next_back_counters = NULL;
  int next_back_valid = -4; // initialized to an invalid value
  bool next_moving_right = false;

  if (valid_dir == NO_SUCCESSORS) {  
	 node->flags |= REMOVE_NEXT_VALUE;
	 DISCONNECTED_CHUNKS ++;
	 // printf("a disconnected chunk of length %d\n", length);
	 return true;
  }

  // a single successor
  // CONDITION 2: the node in the extension direction has a single
  // neighbor in the direction facing this node
  assert((valid_dir >= 0) && (valid_dir < 4));
  bool sense_changed = counters[valid_dir] < 0;
  next = pointers[valid_dir];
  next_moving_right = moving_right ^ sense_changed;
  // next_back_pointers = next_moving_right ? next->left : next->right;
  next_back_counters = next_moving_right ? next->left_count : next->right_count;
  next_back_valid = valid_single_successor(next_back_counters);

  if (next_back_valid != MULTIPLE_SUCCESSORS) {
	 bool remove = check_for_tip_removal(next, next_moving_right, length);
	 if (remove) {
		node->flags |= REMOVE_NEXT_VALUE;
	 }
	 return remove;
  }
  
  // can't extend.  Need to determine if we should trim.
  counter_t max_value = get_counter_max(next_back_counters);
  counter_t value = abs(counters[valid_dir]);
  if ((value >= g_minimum_edge_weight) && (value == max_value)) {
	 // printf("dominant connection: length %d, weight %d\n", length, abs(counters[valid_dir]));
	 return false;   // can't remove the dominant connection
  }

  // if we are clipping the tip, then mark the current "node" for
  // removal and disconnect the "next" node from the tip.  

  // Note: we remove the link from "next", but not from "node".  The
  // reason why is that we don't want code to start on the other side
  // thinking that it is a tip that might need removal.  This means 
  // that from "node's" perspective, the graph looks inconsistent, but
  // that is okay, since we're going to remove "node" shortly.
  node->flags |= REMOVE_NEXT_VALUE;
  // printf("removal candidate: length %d, weight %d\n", length, abs(counters[valid_dir]));
  next_back_valid = (moving_right ? (int)(node->kmer >> ((g__FULLKMER_LENGTH<<1)-2)) : (int)(node->kmer & 0x3));
  if (sense_changed) { next_back_valid ^= 0x3; }
  assert(next_back_valid == find_connection_index(next, node, next_moving_right));
  ++ TIP_DETACHED;

#ifdef VERIFY
  verify_node(node, g__FULLKMER_LENGTH);     // FINE, BEFORE WE LET IT GET INCONSISTENT.
  verify_node(next, g__FULLKMER_LENGTH);
#endif 
  kg_node_t **next_back_pointers = next_moving_right ? next->left : next->right;
  next_back_pointers[next_back_valid] = NULL;
  next_back_counters[next_back_valid] = 0;
  next->connections &= ~(1 << (next_back_valid + (next_moving_right ? 4 : 0)));
#ifdef VERIFY
  // verify_node(node, g__FULLKMER_LENGTH);  // THIS VALIDATION WOULD FAIL, SINCE WE DON'T MODIFY NODE
  verify_node(next, g__FULLKMER_LENGTH);
#endif 
  return true;
}
#endif
*/
namespace {

struct serial_remove_tips_functor {
    serial_remove_tips_functor(KmerGraph *hashtable, bool *modified) : hashtable(hashtable), modified(modified) {}
    void operator()(KmerNode *node) {
#ifdef SMALL_NODES
        unsigned left_count = 0, right_count = 0;
        node->getNeighborCounts(&left_count, &right_count);
#else
        // count the number of left and right neighbors
        unsigned connections = node->connections;
        unsigned left_count = CONNECTION_COUNT_MAP[LEFT(connections)];
        unsigned right_count = CONNECTION_COUNT_MAP[RIGHT(connections)];
#endif
        unsigned total_count = left_count + right_count;
        if (total_count <= 1) {  // is a tip; check if it should be removed
            *modified |= check_for_tip_removal(node, (left_count == 0), true, 0);
        } 
    }
  private:
    KmerGraph *hashtable;
    bool *modified;
};

} // namespace: anonymous

// serial remove_tips
void
remove_tips(KmerGraph* hashtable) {
  bool modified = true;
  unsigned pass = 1;

  if (g__NO_TIP_CLIPPING) return;

/*#ifdef VELOUR_TBB
  tbb::tick_count time0, time1;
#endif*/

  while (modified) {
	 modified = false;
	 printf("==================== TIP REMOVAL PASS %d ====================\n", pass++);
/*#ifdef VELOUR_TBB
	 time0 = tbb::tick_count::now();
#endif*/

     serial_remove_tips_functor f(hashtable, &modified);
     kg_for_each_mutable(hashtable, f);

/*#ifdef VELOUR_TBB
	 time1 = tbb::tick_count::now();
#endif*/
	 printf("tip_removal: %"PRIuPTR" tips detached, %"PRIuPTR" nodes deallocated, %"PRIuPTR" floating chunks, %"PRIuPTR" boundary arcs cut\n", 
			  TIP_DETACHED, NODES_DEALLOCATED, DISCONNECTED_CHUNKS, BOUNDARY_CUTS);
/*#ifdef VELOUR_TBB
	 printf("... pass took %lfs.\n", (time1-time0).seconds());
#endif*/
	 flip_to_remove_value();
  }
}
/*
// TBB parallel remove_tips
#ifdef VELOUR_TBB
static byte_t remove_tips_modified = true;

class RemoveTipsBody {
 public:
  kg_node_t** hashtable;
  
 RemoveTipsBody(kg_node_t** _hashtable) : hashtable(_hashtable) {};
  
  void operator()(const tbb::blocked_range<int>& range) const {
    for(int i = range.begin(); i < range.end(); ++i) {
      kg_node_t** prev = &hashtable[i];
      //_mm_prefetch((hashtable[i+1]), _MM_HINT_NTA);
      for(kg_node_t *trav = *prev; trav != NULL; trav = *prev) {
	union {
	  struct {
	    counter_t left_count[4];
	    counter_t right_count[4];
	  };
	  unsigned long long connections;
	} local_counters;

	//_mm_prefetch(trav->next, _MM_HINT_NTA);
	local_counters.connections = axtomics::fetch_and_nop<unsigned long long>(&trav->connections);
	bool remove_now = (local_counters.connections == 0);
	if (remove_now) {
	  axtomics::fetch_and_add<uintptr_t>(&NODES_DEALLOCATED, 1);
	  *prev = trav->next; 
	  velour_free(trav);
	  continue;
	}

	unsigned left_count = 0, right_count = 0;
	trav->getNeighborCounts(&left_count, &right_count);
	unsigned total_count = left_count + right_count;
	if (total_count == 0) {  // dangling
	  // printf("Dangling k-mer!!!\n");
	  axtomics::fetch_and_add<uintptr_t>(&NODES_DEALLOCATED, 1);
	  *prev = trav->next; 
	  velour_free(trav);
	  continue;
	} 

	if (total_count == 1) {  // is a tip; check if it should be removed
	  //modified |= check_for_tip_removal(trav, (left_count == 0), 0);
	  axtomics::fetch_and_or<byte_t>(&remove_tips_modified, check_for_tip_removal(trav, (left_count == 0), 0));
	} 
	// keep track of pointer to next "trav" node, so that we can do removals
	// as necessary.
	prev = &(trav->next); 
      } // for this bucket
    } // for my buckets
  } // operator()
};

// Marking pass.  Buckets are broken up among threadsby TBB, and each thread is 
//  responsible for marking tips that end in their bucket.  
class MarkTipsBody {
 public:
  kg_node_t** hashtable;
  
 MarkTipsBody(kg_node_t** _hashtable) : hashtable(_hashtable) {};

  void operator()(const tbb::blocked_range<int>& range) const {
    tbb::tick_count time0, time1;
    // traverse hash buckets, marking nodes to free
    for(int i = range.begin(); i < range.end(); ++i) {
      kg_node_t **prev = &hashtable[i];
      //_mm_prefetch((hashtable[i+1]), _MM_HINT_NTA);
      for(kg_node_t *trav = *prev; trav != NULL; trav = *prev) {
	//_mm_prefetch(trav->next, _MM_HINT_NTA);

	unsigned left_count = 0, right_count = 0;
	trav->getNeighborCounts(&left_count, &right_count);
	unsigned total_count = left_count + right_count;

	if (total_count == 1) {  // is a tip; check if it should be removed
	  //remove_tips_modified |= check_for_tip_removal(trav, (left_count == 0), 0);
	  if(check_for_tip_removal(trav, (left_count == 0), 0)) {
	    remove_tips_modified = true;
	  }
	} 
	// keep track of pointer to next "trav" node, so that we can do removals
	// as necessary.
	prev = &(trav->next); 
      }
    }
    //    flip_to_remove_value();
  }
};

class FreeTipsBody {
 public:
  kg_node_t** hashtable;

 FreeTipsBody(kg_node_t** _hashtable) : hashtable(_hashtable) {};

  void operator()(const tbb::blocked_range<int>& range) const {
    //traverse hash buckets for free-able nodes
    for(int i = range.begin(); i < range.end(); ++i) {
      kg_node_t **prev = &hashtable[i];
      for(kg_node_t *trav = *prev; trav != NULL; trav = *prev) {
	
	bool remove_now = (trav->connections == 0);
	if (remove_now) {
	  //NODES_DEALLOCATED ++;
	  axtomics::fetch_and_add<unsigned>(&NODES_DEALLOCATED, 1);
	  *prev = trav->next; 
	  velour_free(trav);
	  remove_tips_modified = true;
	  continue;
	}

	unsigned left_count = 0, right_count = 0;
	trav->getNeighborCounts(&left_count, &right_count);
	unsigned total_count = left_count + right_count;

	if (total_count == 0) {  // dangling
	  // printf("Dangling k-mer!!!\n");
	  //NODES_DEALLOCATED ++;
	  axtomics::fetch_and_add<unsigned>(&NODES_DEALLOCATED, 1);
	  *prev = trav->next; 
	  velour_free(trav);
	  remove_tips_modified = true;
	  continue;
	}
	prev = &(trav->next);
      } // for this bucket
    } // for my buckets
  } // operator()
};

void
parallel_for_remove_tips(kg_node_t **hashtable) {
  unsigned pass = 1;
  tbb::tick_count time0, time1;
  remove_tips_modified = true;
  HASHTABLE = hashtable;  // fixme: yucky global hack (only needed for small prenodes)  

  while(remove_tips_modified) {
    remove_tips_modified = false;
    printf("==================== TIP REMOVAL PASS %d ====================\n", pass++);
    time0 = tbb::tick_count::now();
    //    tbb::parallel_for(tbb::blocked_range<int>(0,HASH_BUCKETS,1),
    //		      MarkTipsBody(hashtable));
    //    tbb::parallel_for(tbb::blocked_range<int>(0,HASH_BUCKETS,1),
    //		      FreeTipsBody(hashtable));
    tbb::parallel_for(tbb::blocked_range<int>(0,HASH_BUCKETS,1),
		      RemoveTipsBody(hashtable), tbb::auto_partitioner());
    time1 = tbb::tick_count::now();
    printf("tip_removal: %d tips detached, %d nodes deallocated, %d floating chunks\n", 
	   TIP_DETACHED, NODES_DEALLOCATED, DISCONNECTED_CHUNKS);
    printf("... pass took %lfs.\n", (time1-time0).seconds());
    flip_to_remove_value();
  }
}
#endif // VELOUR_TBB
*/

