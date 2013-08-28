//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// seqNode.cpp
//

#include "types.h"
#include "seqNode.h"

// TODO: add usage verify assertions to check that the head/tail entry does not exist
SeqNode * GrowSeqNode(SeqGraph * sgraph, SeqNode * orig_node, unsigned bases, bool skip_head_update, bool skip_tail_update) // TODO move this to SeqGraph
{
    // first, remove the orig_node from the graph -- XXX is this safe re: iterators???
    if (!skip_head_update) { sgraph->removeNodeHead(orig_node); }
    if (!skip_tail_update) { sgraph->removeNodeTail(orig_node); }
 
    // then, reallocate node
    assert( static_cast<unsigned>(orig_node->sequence.get_length()) + bases <= Sequence::MAX_BASES );
    uint16_t new_seq_length = orig_node->sequence.get_length() + bases;
    size_t new_size = SeqNode::ComputeNodeAllocatedBytes(new_seq_length);

#ifdef VELOUR_TBB
    SeqNode * new_node = tls_seqnode_allocator->ReallocateNodeMemory(orig_node, new_size);
#else
    SeqNode * new_node = g__SEQNODE_ALLOCATOR->ReallocateNodeMemory(orig_node, new_size);
#endif // VELOUR_TBB

#ifndef NDEBUG
    // fixup new allocation length
    new_node->sequence.alloc_length_ = Sequence::ComputeAllocatedBases(new_seq_length);
#endif

    // insert new node into sequence graph
    if (!skip_head_update) { sgraph->insertNodeHead(new_node); }
    if (!skip_tail_update) { sgraph->insertNodeTail(new_node); }

    // finally, return the reallocated node
    return new_node;
}

size_t SeqNode::GetNodeSerializedBytes() const
{
    size_t serialized_bytes = 0;

    serialized_bytes += sizeof(this->left_side);
    serialized_bytes += sizeof(this->right_side);
    serialized_bytes += sizeof(this->left_side_colors);
    serialized_bytes += sizeof(this->right_side_colors);
    serialized_bytes += sizeof(this->flags);
    serialized_bytes += sizeof(this->kmer_occurrences);

    serialized_bytes += this->sequence.GetSerializedBytes();

    return serialized_bytes;
}

void SeqNode::loadFromFile(FILE* file, int fileFormat)
{
	size_t retval;
    assert(file != NULL);

    switch (fileFormat) {
    /*case PREGRAPH:
        {
            long int dummy;
            fscanf(file, "NODE\t%li\n", &dummy);
            const int PREGRAPH_MAXLINE = 1000000; // FIXME constant
            char buf[PREGRAPH_MAXLINE];
            char *line = fgets(buf, PREGRAPH_MAXLINE, file);
            assert( line != NULL );
            char c;
            while ((c = *(line++)) != '\n') {
                Nucleotide n = BASE_MAP[int(c)];
                seq_append_base(&sg_node->sequence, n);  // FIXME: be more efficient
            }
            // FIXME: not creating connections?!
        }
        break; */
    case QUILT:
    case BUCKET:
        // TODO: any missing fields
        if (fread(&this->left_side, sizeof(this->left_side), 1, file) != 1) {
			assert( feof(file) && "fread() failed but not end-of-file!?" );
            return;
		}
        retval = fread(&this->right_side, sizeof(this->right_side), 1, file);
		assert(retval == 1 && "fread() error." );
        retval = fread(&this->left_side_colors, sizeof(this->left_side_colors), 1, file);
		assert(retval == 1 && "fread() error." );
        retval = fread(&this->right_side_colors, sizeof(this->right_side_colors), 1, file);
		assert(retval == 1 && "fread() error." );
        retval = fread(&this->flags, sizeof(this->flags), 1, file);
		assert(retval == 1 && "fread() error." );
        retval = fread(&this->kmer_occurrences, sizeof(this->kmer_occurrences), 1, file);
		assert(retval == 1 && "fread() error." );
        this->sequence.Load_BinaryFile(file);
        break;
    default:
        assert( false && "Invalid input format for sequence node." );
        fprintf(stderr, "Unrecoverable error.  Exiting...\n");
        exit(1);
    }
}

size_t SeqNode::emitToFile(FILE *file, int fileFormat, uintptr_t nodeIndex)
{
    size_t total_size = 0;
	size_t retval;
    assert(file != NULL);
    assert( !isNodeDead(this) );
    switch (fileFormat) {
    case PREGRAPH:
        fprintf(file, "NODE\t%"PRIuPTR"\t%hu\t%d\t%g\n", nodeIndex, (this->sequence.get_length()-g__FULLKMER_LENGTH+1), 0, this->getNodeKmerCoverage(g__FULLKMER_LENGTH)); // FIXME: PRIu64 when have offset to add to in future // FIXME: support non-16bit node length???
        //fprintf(file, "NODE\t%"PRIuPTR"\t%hu\n", nodeIndex, (this->sequence.get_length()-g__FULLKMER_LENGTH+1)); // FIXME: PRIu64 when have offset to add to in future // FIXME: support non-16bit node length???
        this->sequence.PrintToFile(file);
        fputc('\n', file);
        break;
    case QUILT:
    case BUCKET:
        // FIXME: use a variable-length encoding???
        // TODO: pseudo-bits?
        retval = fwrite(&this->left_side, sizeof(this->left_side), 1, file);
		assert(retval == 1 && "fwrite() error." );
        total_size += sizeof(this->left_side);
        retval = fwrite(&this->right_side, sizeof(this->right_side), 1, file);
		assert(retval == 1 && "fwrite() error." );
        total_size += sizeof(this->right_side);
        retval = fwrite(&this->left_side_colors, sizeof(this->left_side_colors), 1, file);
		assert(retval == 1 && "fwrite() error." );
        total_size += sizeof(this->left_side_colors);
        retval = fwrite(&this->right_side_colors, sizeof(this->right_side_colors), 1, file);
		assert(retval == 1 && "fwrite() error." );
        total_size += sizeof(this->right_side_colors);
        // TODO: bitconnections
        // TODO: unique
        retval = fwrite(&this->flags, sizeof(this->flags), 1, file); // TODO: clean flags?
		assert(retval == 1 && "fwrite() error." );
        total_size += sizeof(this->flags);
        retval = fwrite(&this->kmer_occurrences, sizeof(this->kmer_occurrences), 1, file);
		assert(retval == 1 && "fwrite() error." );
        total_size += sizeof(this->kmer_occurrences);
        total_size += this->sequence.Save_BinaryFile(file);
        break;
    default:
        assert( false && "Invalid output format for sequence node." );
        fprintf(stderr, "Unrecoverable error.  Exiting...\n");
        exit(1);
    }
    return total_size;
}

size_t SeqNode::loadFromMemory(char *memory)
{
    off_t offset = 0;

    memcpy(&this->left_side, (memory+offset), sizeof(this->left_side));
    offset += sizeof(this->left_side);

    memcpy(&this->right_side, (memory+offset), sizeof(this->right_side));
    offset += sizeof(this->right_side);

    memcpy(&this->left_side_colors, (memory+offset), sizeof(this->left_side_colors));
    offset += sizeof(this->left_side_colors);

    memcpy(&this->right_side_colors, (memory+offset), sizeof(this->right_side_colors));
    offset += sizeof(this->right_side_colors);

    memcpy(&this->flags, (memory+offset), sizeof(this->flags));
    offset += sizeof(this->flags);

    memcpy(&this->kmer_occurrences, (memory+offset), sizeof(this->kmer_occurrences));
    offset += sizeof(this->kmer_occurrences);

    offset += this->sequence.Load_Memory(static_cast<void *>(memory+offset));

    return offset;
}

size_t SeqNode::saveToMemory(char *memory, size_t space_remaining)
{
    off_t offset = 0;
    assert( !isNodeDead(this) );

    memcpy((memory+offset), &this->left_side, sizeof(this->left_side));
    offset += sizeof(this->left_side);

    memcpy((memory+offset), &this->right_side, sizeof(this->right_side));
    offset += sizeof(this->right_side);

    memcpy((memory+offset), &this->left_side_colors, sizeof(this->left_side_colors));
    offset += sizeof(this->left_side_colors);

    memcpy((memory+offset), &this->right_side_colors, sizeof(this->right_side_colors));
    offset += sizeof(this->right_side_colors);

    memcpy((memory+offset), &this->flags, sizeof(this->flags)); // TODO: clean flags?
    offset += sizeof(this->flags);

    memcpy((memory+offset), &this->kmer_occurrences, sizeof(this->kmer_occurrences));
    offset += sizeof(this->kmer_occurrences);

    offset += this->sequence.Save_Memory(static_cast<void *>(memory+offset), (space_remaining - offset));

    return offset;
}

void SeqNode::getNeighborCounts(unsigned *l_count, unsigned *r_count) // const
{
    unsigned left_count = 0, right_count = 0;
#ifdef VELOUR_TBB
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

