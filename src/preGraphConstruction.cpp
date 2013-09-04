//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// preGraphConstruction.cpp
//

#include "types.h"


//**************************    ALLOCATE PRE-NODES    ****************************
//unsigned PRENODE_ALLOCATED = 0;
//unsigned PRENODE_UNIQUE = 0;


// Find all of the Kmer's for each read; for each Kmer make sure that
// there is a node that has been allocated, attach it to the previous Kmer,
// and increment the counters for the edges.
#ifdef SMALL_NODES
static void 
saturating_update(counter_t *counter, counter_t increment) {
#ifdef VELOUR_TBB
    while(1) {
        counter_t temp = ATOMIC_LOAD(*counter);
        if (abs(temp) < MAX_COUNTER_VALUE) {
            counter_t old_val = ATOMIZE(*counter).compare_and_swap(temp + increment, temp);
            if(old_val == temp) { // no intervening updates
                break;
            }
        } else { // don't update
            break;
        }
    }
#else // VELOUR_TBB
    if (abs(*counter) < MAX_COUNTER_VALUE) {
        *counter += increment;
    }
#endif
}

namespace {
    struct lambda_kmers_to_prenodes
    {
        lambda_kmers_to_prenodes(KmerGraph *kgraph) : kgraph_(kgraph), node_(NULL), sense_reversed_(false) {}
        void operator()(Kmer k, Kmer canon_k, bool new_sense_reversed, Nucleotide base, Nucleotide last_base, unsigned kmer_length, bool isFirstKmer, bool isLastKmer)
        {
            if (isFirstKmer) {
                node_ = kgraph_->findNodeOrAllocate(canon_k);
                sense_reversed_ = new_sense_reversed;
                return;
            }

            KmerNode *new_node = kgraph_->findNodeOrAllocate(canon_k);

            // Positive counts between nodes of the same sense, negative if sense is reversed.
            int increment = (sense_reversed_ == new_sense_reversed) ? 1 : -1;

            // We interpret the reads as left-to-right sequences.
            // Thus, if the canonical form of the previous node (called
            // "node") is reversed from how it shows up in the read, then
            // extend "node" to the left.
            if (sense_reversed_) {
                saturating_update(&node_->left_count[COMPLEMENT(base)], increment);
            } else {
                saturating_update(&node_->right_count[base], increment);
            }

            if (new_sense_reversed) {
                saturating_update(&new_node->right_count[COMPLEMENT(last_base)], increment);
            } else {
                saturating_update(&new_node->left_count[last_base], increment);
            }

#ifdef VERIFY
            verify_node(node_, kgraph_, kmer_length);
            verify_node(new_node, kgraph_, kmer_length);
#endif

            // make the "new_node" the previous "node"
            node_ = new_node;
            sense_reversed_ = new_sense_reversed;
            return;
        }
        private:
            KmerGraph *kgraph_;
            KmerNode *node_;
            bool sense_reversed_;
    };
}

void convertSequenceToKmersToPrenodes(char *seq, KmerGraph *hashtable, int kmer_length)
{
    Sequence_StackAllocated memory;
    Sequence *stack_seq = new (&memory) Sequence(Sequence::MAX_BASES);
    stack_seq->InitializeWithString(seq);

    lambda_kmers_to_prenodes functor(hashtable);
    sequence_process_kmers(stack_seq, kmer_length, functor);
}
#else
// Find all of the Kmer's for each read; for each Kmer make sure that
// there is a node that has been allocated, attach it to the previous Kmer,
// and increment the counters for the edges.
void
convertSequenceToKmersToPrenodes(char *seq, kg_node_t **hashtable, int kmer_length) {
    assert( false && "TODO FIX! REVERSED KMER ENDIANNESS" );
  Kmer kmer(0), anti_kmer(0);
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
#error  Kmer mask = (((Kmer)1) << double_kmer_length) - 1;  // a mask with 2*kmer_length bits

  // lookup first Kmer to get first node
  Kmer rc_kmer = (anti_kmer >> (64 - double_kmer_length));
  bool sense_reversed = rc_kmer < kmer;  
  Kmer canonical_kmer = sense_reversed ? rc_kmer : kmer;
  unsigned hash_index = HASH_KMER(canonical_kmer);
  kg_node_t *node = findNodeOrAllocate(&hashtable[hash_index], canonical_kmer);

  // make each succeeding Kmer
  char c;
  while ((c = seq[i]) != 0) {
    ++ i;

    // read the next base and extend both the kmer and the anti_kmer
    Nucleotide base = BASE_MAP[(int)c];
    anti_kmer >>= 2;
    kmer <<= 2;
    assert(base <= 0x3);
    anti_kmer |= ((Kmer)base ^ 0x3) << 62;
    kmer |= base;

    // which base just fell off the kmer?
    Nucleotide last_base = (kmer >> double_kmer_length) & 0x3;
    kmer &= mask;
    kg_node_t *new_node = NULL; 
    bool new_sense_reversed = false;

#ifdef VERIFY
    Kmer new_rc_kmer = (anti_kmer >> (64 - double_kmer_length));
    bool new_sense_reversed_check = new_rc_kmer < (kmer & mask);  
    canonical_kmer = new_sense_reversed_check ? new_rc_kmer : (kmer & mask);
    hash_index = HASH_KMER(canonical_kmer);
    kg_node_t *new_node_check = findNodeOrAllocate(&hashtable[hash_index], canonical_kmer);
#endif

    // We interpret the reads as left-to-right sequences.
    // Thus, if the canonical form of the previous node (called
    // "node") is reversed from how it shows up in the read, then
    // extend "node" to the left.
    if (sense_reversed) {
      Nucleotide new_rc_base = base ^ 0x3;
      new_node = node->left[(int) new_rc_base];
      if (new_node != NULL) { // "new node" already connected to previous "node"
		  VERIFY_ASSERT(new_node == new_node_check);
		  assert(new_node->kmer == canonical_kmer);
		  int count = 0;
		  while((count = node->left_count[(int) new_rc_base]) == 0); // spin until resovled
		  bool change_in_sense = count < 0;
		  new_sense_reversed = !change_in_sense;
		  VERIFY_ASSERT(new_sense_reversed == new_sense_reversed_check);
		  if (change_in_sense) { // change in sense between this node and the next, decrement both on same side
			 VERIFY_ASSERT(!new_sense_reversed_check);
			 //node->left_count[(int)new_rc_base] = count - 1;
			 axtomics::fetch_and_add<counter_t>(&(node->left_count[(int)new_rc_base]), -1);
			 //-- new_node->left_count[(int)last_base];
			 axtomics::fetch_and_add<counter_t>(&(new_node->left_count[(int)last_base]), -1);
		  } else { // no change in sense between these nodes, increment both on opposite sides
			 VERIFY_ASSERT(new_sense_reversed_check);
			 //node->left_count[(int)new_rc_base] = count + 1;
			 axtomics::fetch_and_add<counter_t>(&(node->left_count[(int)new_rc_base]), 1);
			 //++ new_node->right_count[(int)last_base ^ 0x3];
			 axtomics::fetch_and_add<counter_t>(&(new_node->right_count[(int)last_base ^ 0x3]), 1);
		  }
      } else {
		  Kmer new_rc_kmer = (anti_kmer >> (64 - double_kmer_length));
		  new_sense_reversed = new_rc_kmer < (kmer & mask);  
		  canonical_kmer = new_sense_reversed ? new_rc_kmer : (kmer & mask);
		  hash_index = HASH_KMER(canonical_kmer);
		  new_node = findNodeOrAllocate(&hashtable[hash_index], canonical_kmer);

		  node->left[(int)new_rc_base] = new_node;
		  // node->connections |= 1 << (LEFT_CONNECT_OFFSET + (int)new_rc_base);
		  axtomics::fetch_and_or<byte_t>(&(node->connections), 1 << (LEFT_CONNECT_OFFSET + (int)new_rc_base));
		  bool change_in_sense = !new_sense_reversed; // sense_reversed ^ new_sense_reversed_check;
		  if (change_in_sense) { 
		    //node->left_count[(int)new_rc_base] = -1;
		    axtomics::fetch_and_add<counter_t>(&(node->left_count[(int)new_rc_base]), -1);
		    new_node->left[(int)last_base] = node;
		    //new_node->connections |= 1 << (LEFT_CONNECT_OFFSET + last_base);
		    axtomics::fetch_and_or<byte_t>(&(new_node->connections), 1 << (LEFT_CONNECT_OFFSET + last_base));
		    //new_node->left_count[(int)last_base] = -1;
		    axtomics::fetch_and_add<counter_t>(&(new_node->left_count[(int)last_base]), -1);
		  } else { 
		    //node->left_count[(int)new_rc_base] = 1;
		    axtomics::fetch_and_add<counter_t>(&(node->left_count[(int)new_rc_base]), 1);
		    new_node->right[(int)last_base ^ 0x3] = node;
		    //new_node->connections |= 1 << (last_base ^ 0x3);
		    axtomics::fetch_and_or<byte_t>(&(new_node->connections), 1 << (last_base ^ 0x3));
		    //new_node->right_count[(int)last_base ^ 0x3] = 1;
		    axtomics::fetch_and_add<counter_t>(&(new_node->right_count[(int)last_base ^ 0x3]), 1);
		  }
      }

    } else {  // otherwise, extend "node" to the right.
      new_node = node->right[(int)base];
      if (new_node != NULL) { // "new node" already connected to previous "node"
		  VERIFY_ASSERT(new_node == new_node_check);
		  assert(new_node->kmer == canonical_kmer);
		  int count = 0;
		  while((count = node->right_count[base]) == 0); // spin until resolved
		  bool change_in_sense = count < 0;
		  new_sense_reversed = change_in_sense;
		  VERIFY_ASSERT(new_sense_reversed == new_sense_reversed_check);
		  if (change_in_sense) { // change in sense between this node and the next, decrement both on same side
			 VERIFY_ASSERT(new_sense_reversed_check);
			 //node->right_count[base] = count - 1;
			 axtomics::fetch_and_add<counter_t>(&(node->right_count[base]), -1);
			 //-- new_node->right_count[last_base ^ 0x3];
			 axtomics::fetch_and_add<counter_t>(&(new_node->right_count[last_base ^ 0x3]), -1);
		  } else { // no change in sense between these nodes, increment both on opposite sides
			 VERIFY_ASSERT(!new_sense_reversed_check);
			 //node->right_count[base] = count + 1;
			 axtomics::fetch_and_add<counter_t>(&(node->right_count[base]), 1);
			 //++ new_node->left_count[last_base];
			 axtomics::fetch_and_add<counter_t>(&(new_node->left_count[last_base]), 1);
		  }
      } else { // new_node is NULL
		  Kmer new_rc_kmer = (anti_kmer >> (64 - double_kmer_length));
		  new_sense_reversed = new_rc_kmer < (kmer & mask);  
		  canonical_kmer = new_sense_reversed ? new_rc_kmer : (kmer & mask);
		  hash_index = HASH_KMER(canonical_kmer);
		  new_node = findNodeOrAllocate(&hashtable[hash_index], canonical_kmer);

		  node->right[base] = new_node;
		  //node->connections |= (1 << base);
		  axtomics::fetch_and_or<byte_t>(&(node->connections), (1 << base));
		  bool change_in_sense = new_sense_reversed; // sense_reversed ^ new_sense_reversed;
		  if (change_in_sense) { 
		    //node->right_count[base] = -1;
		    axtomics::fetch_and_add<counter_t>(&(node->right_count[base]), -1);
		    new_node->right[last_base ^ 0x3] = node;
		    //new_node->connections |= (1 << (last_base ^ 0x3));
		    axtomics::fetch_and_or<byte_t>(&(new_node->connections), (1 << (last_base ^ 0x3)));
		    //new_node->right_count[last_base ^ 0x3] = -1;
		    axtomics::fetch_and_add<counter_t>(&(new_node->right_count[last_base ^ 0x3]), -1);
		  } else { 
		    //node->right_count[base] = 1;
		    axtomics::fetch_and_add<counter_t>(&(node->right_count[base]), 1);
		    new_node->left[last_base] = node;
		    //new_node->connections |= 1 << (LEFT_CONNECT_OFFSET + last_base);
		    axtomics::fetch_and_or<byte_t>(&(new_node->connections), 1 << (LEFT_CONNECT_OFFSET + last_base));
		    //new_node->left_count[last_base] = 1;
		    axtomics::fetch_and_add<counter_t>(&(new_node->left_count[last_base]), 1);
		  }
      }
    }

#ifdef VERIFY
    verify_node(node, kmer_length);
    verify_node(new_node, kmer_length);
#endif

    // make the "new_node" the previous "node"
    node = new_node;
    sense_reversed = new_sense_reversed;
  }
}
#endif

void print_hashtable_histogram(KmerGraph *kgraph)
{
    Histogram hist(65,1,0);

    for(uintptr_t i = 0 ; i < kgraph->buckets ; ++ i) {
        unsigned count = 0;
        for (KmerNode *trav = kgraph->table[i] ; trav != NULL ; trav = trav->next) {
            ++ count;
        }
        hist.insert(count);
    }

    hist.print("KG ", 1, 1);
}

#ifdef SMALL_NODES
namespace {
    struct lambda_subsequences
    {
        lambda_subsequences(KmerGraph *kgraph, unsigned hps, Color pc, Color sc) :
            kgraph_(kgraph), hasPrefix_(hps & 0x1), hasSuffix_(hps & 0x2), prefixColor_(pc), suffixColor_(sc), node_(NULL) {}
        void operator()(Kmer k, Kmer canon_k, bool new_sense_reversed, Nucleotide base, Nucleotide last_base, unsigned kmer_length, bool isFirstKmer, bool isLastKmer)
        {
            if (isFirstKmer) {
                node_ = NULL;
                if( !hasPrefix_ ) {
                    node_ = kgraph_->findNodeOrAllocate(canon_k);
                }
                sense_reversed_ = new_sense_reversed;
                return;
            }

            KmerNode *new_node = NULL;
            if( !isLastKmer || !hasSuffix_ )
            {
                new_node = kgraph_->findNodeOrAllocate(canon_k);
            }

            // Positive counts between nodes of the same sense, negative if sense is reversed.
            int increment = (sense_reversed_ == new_sense_reversed) ? 1 : -1;

            // We interpret the reads as left-to-right sequences.
            // Thus, if the canonical form of the previous node (called
            // "node") is reversed from how it shows up in the read, then
            // extend "node" to the left.
            if( node_ != NULL )
            {
                if (sense_reversed_) {
                    saturating_update(&node_->left_count[COMPLEMENT(base)], increment);
                    if (new_node == NULL) {
                        assert( hasSuffix_ );
                        ATOMIC_STORE(node_->left_color[COMPLEMENT(base)]) = suffixColor_;
                    }
                } else {
                    saturating_update(&node_->right_count[base], increment);
                    if (new_node == NULL) {
                        assert( hasSuffix_ );
                        ATOMIC_STORE(node_->right_color[base]) = suffixColor_;
                    }
                }
            }

            if( new_node != NULL )
            {
                if (new_sense_reversed) {
                    saturating_update(&new_node->right_count[COMPLEMENT(last_base)], increment);
                    if (node_ == NULL) {
                        assert( hasPrefix_ );
                        ATOMIC_STORE(new_node->right_color[COMPLEMENT(last_base)]) = prefixColor_;
                    }
                } else {
                    saturating_update(&new_node->left_count[last_base], increment);
                    if (node_ == NULL) {
                        assert( hasPrefix_ );
                        ATOMIC_STORE(new_node->left_color[last_base]) = prefixColor_;
                    }
                }
            }

#ifdef VERIFY
            if( node_ != NULL ) { verify_node(node_, kgraph_, kmer_length); }
            if( new_node != NULL ) { verify_node(new_node, kgraph_, kmer_length); }
#endif

            // make the "new_node" the previous "node"
            node_ = new_node;
            sense_reversed_ = new_sense_reversed;
            return;
        }
        private:
            KmerGraph *kgraph_;
            bool hasPrefix_;
            bool hasSuffix_;
            Color prefixColor_;
            Color suffixColor_;

            KmerNode *node_;
            bool sense_reversed_;
    };
}

void convertSubsequenceToKmersToPrenodes(Sequence *seq, KmerGraph *hashtable, int kmer_length, unsigned hasPrefixSuffix, Color prefixColor, Color suffixColor)
{
    lambda_subsequences functor(hashtable, hasPrefixSuffix, prefixColor, suffixColor);
    sequence_process_kmers(seq, kmer_length, functor);
}
#endif

