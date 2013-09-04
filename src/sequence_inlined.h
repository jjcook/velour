//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// sequence_inlined.h
//
//   the inlined member functions for class Sequence
//
//   ** IMPORTANT NOTE **
//     within 'xxx_Unsafe' functions cannot use length_ to derive allocated space
//
//   FIXME NOTE: word-sized manipulations assume sequence length allocations in multiples of words!

#ifndef VELOUR_SEQUENCE_INLINED_H
#define VELOUR_SEQUENCE_INLINED_H

namespace {

template<typename T>
static inline unsigned ComputeWordIndex(uint16_t bases)
{
    return (bases * 2) / (sizeof(T) * 8);
}

template<typename T>
static inline unsigned ComputeWordBitOffset(uint16_t bases)
{
    return (bases * 2) % (sizeof(T) * 8);
}

} // namespace: anonymous

extern Nucleotide BASE_MAP[256];
extern char CHAR_BASE_MAP[4];

// XXX initialize buffer to silence valgrind re uninitialized memory (SetBase reads before writes)
inline void Sequence::Init(uint16_t alloc_bases)
{
    bzero(str_, (alloc_bases >> 2)); // FIXME calculation
}

inline uint16_t Sequence::get_length(void) const
{
    return length_;
}

// NOTE: used by memory allocator to pseudo-init when marking unallocated node dead
inline void Sequence::zero_length()
{
    length_ = 0;
}
 
inline bool Sequence::CanAddBases(unsigned count) const
{
#ifdef VERIFY
    assert( count > 0 );
#endif
    return ((static_cast<unsigned>(length_) + count) <= MAX_BASES);
}
        
inline size_t Sequence::ComputeLengthInBytes(uint16_t bases)
{
    return (static_cast<unsigned>(bases) + 3) >> 2;  // round up
}

inline size_t Sequence::ComputeAllocatedBytes(uint16_t bases)
{
    size_t allocated_bytes = max(MIN_BYTES, round_up_power_of_2( ComputeLengthInBytes(bases) ));
    assert( allocated_bytes % sizeof(uintptr_t) == 0 ); // word-sized internal manipulations assume allocations in multiples of words!
    return allocated_bytes;
}

inline uint16_t Sequence::ComputeAllocatedBases(uint16_t bases)
{
    unsigned amount = ComputeAllocatedBytes(bases) << 2;
    return (amount < UINT16_MAX ? amount : UINT16_MAX);  // catch overflow of actual maximum allocation of UINT16_MAX+1 bases
}

inline size_t Sequence::GetLengthInBytes(void) const
{
    return ComputeLengthInBytes(length_);
}

inline uint16_t Sequence::GetAllocatedBases() const
{
    uint16_t bases = ComputeAllocatedBases(length_);
    assert( bases <= alloc_length_ );
    return bases;
}

inline size_t Sequence::GetAllocatedBytes(void) const
{
    return ComputeAllocatedBytes(length_);
}
        
inline Nucleotide Sequence::GetBase_Unsafe(const uint16_t index) const
{
    assert(index < length_);
    unsigned byte_index = index >> 2;
    unsigned bit_index = (index & 0x3) << 1;
    return (str_[byte_index] >> bit_index) & 0x3;
}

inline Nucleotide Sequence::GetBase(const uint16_t index) const
{
    assert(index < length_);
    return GetBase_Unsafe(index);
}

inline void Sequence::SetBase_Unsafe(const uint16_t index, Nucleotide base)
{
    unsigned byte_index = index >> 2;
	byte_t thebyte = str_[byte_index];

    unsigned bit_index = (index & 0x3) << 1;
	byte_t mask = ~(0x03 << bit_index);
	thebyte &= mask;                // TODO OPT: mask is necessary iff allocated space not zero'd

	thebyte |= (base << bit_index);
	
	str_[byte_index] = thebyte;
}
        
inline void Sequence::AppendBase_Unsafe(Nucleotide base)
{
    assert((base & 0x3) == base);
    unsigned index = length_;
    SetBase_Unsafe(index, base);
    ++ length_;
    assert( length_ <= alloc_length_ );
}

inline void Sequence::AppendBase(SeqGraph *sgraph, SeqNode **node, Nucleotide base)
{
    if (length_ == MAX_BASES) {
        fprintf(stderr, "ERROR: maximum sequence length!");
        assert( false );
        exit(EXIT_FAILURE);
    }

    if (length_ == GetAllocatedBases()) {
        // need to grow node
        *node = GrowSeqNode(sgraph, *node, 1, false, true);
        (*node)->sequence.AppendBase_Unsafe(base);
        return;
    }
    AppendBase_Unsafe(base);
}

inline void Sequence::PreparePrepend_Unsafe(uint16_t insert_length)
{
    assert( (static_cast<unsigned>(length_) + static_cast<unsigned>(insert_length)) <= alloc_length_ );

    if (insert_length == 0) return; // XXX: delete this?
    if (length_ == 0) return;       // XXX: delete this?

#ifdef VERIFY
    Sequence_StackAllocated memory1;
    Sequence *orig = new (&memory1) Sequence(*this, Sequence::MAX_BASES);
    assert( orig->get_length() == length_ );
#endif // VERIFY

    uintptr_t *words = reinterpret_cast<uintptr_t *>(str_); // FIXME FIXME alignment: this will break on SPARC! use memcpy?

    int sourceWordIndex = ComputeWordIndex<uintptr_t>(length_ - 1);
    int targetWordIndex = ComputeWordIndex<uintptr_t>(length_ + insert_length - 1);
    assert( targetWordIndex >= sourceWordIndex );
    assert((targetWordIndex * sizeof(uintptr_t)) < ComputeLengthInBytes(alloc_length_) );

    unsigned sourceBitOffset = ComputeWordBitOffset<uintptr_t>(length_);
    unsigned targetBitOffset = ComputeWordBitOffset<uintptr_t>(length_ + insert_length);
    if (sourceBitOffset == 0) { sourceBitOffset = sizeof(uintptr_t) * 8; }
    if (targetBitOffset == 0) { targetBitOffset = sizeof(uintptr_t) * 8; }

    // NOTE: doesn't zero the prepended placeholder bases
    if (targetBitOffset == sourceBitOffset) {
        assert( targetWordIndex > sourceWordIndex );
        for ( ; sourceWordIndex >= 0 ; --sourceWordIndex, --targetWordIndex ) {
            words[targetWordIndex] = words[sourceWordIndex];
        }
    } else if (targetBitOffset > sourceBitOffset) {
        unsigned shamt = targetBitOffset - sourceBitOffset;
        unsigned anti_shamt = (sizeof(uintptr_t) * 8) - shamt; // ok for shift below, is always < 64

        for ( ; sourceWordIndex > 0 ; --sourceWordIndex, --targetWordIndex ) {
            words[targetWordIndex] = (words[sourceWordIndex] << shamt) | (words[sourceWordIndex-1] >> anti_shamt);
        }
        words[targetWordIndex] = (words[sourceWordIndex] << shamt);
    } else if (targetBitOffset < sourceBitOffset) {
        assert( targetWordIndex > sourceWordIndex );
        unsigned shamt = sourceBitOffset - targetBitOffset;
        unsigned anti_shamt = (sizeof(uintptr_t) * 8) - shamt; // ok for shift below, is always < 64

        words[targetWordIndex] = (words[sourceWordIndex] >> shamt);
        --targetWordIndex;

        for ( ; sourceWordIndex > 0 ; --sourceWordIndex, --targetWordIndex ) {
            words[targetWordIndex] = (words[sourceWordIndex] << anti_shamt) | (words[sourceWordIndex-1] >> shamt);
        }
        words[targetWordIndex] = (words[sourceWordIndex] << anti_shamt);
    }

    /* // original code
    for (int i=length_ - 1 ; i >= 0 ; -- i) {
		Nucleotide base = GetBase_Unsafe(i);
		SetBase_Unsafe(i+insert_length, base);
	}*/

    length_ += insert_length; // set the new length in anticipation of being filled in

    // XXX: if SetBase does not load & mask, then need to zero the prepended space

    /* // clear the prepended space (inefficiently)
    for (int i=0 ; i < insert_length ; ++i) {
        SetBase_Unsafe(i, 0);
    }*/

#ifdef VERIFY
	for (int i=0 ; i < orig->get_length() ; ++i) {
		Nucleotide good_base = orig->GetBase_Unsafe(i);
        Nucleotide test_base = GetBase_Unsafe(i+insert_length);
        assert( test_base == good_base );
    }
#endif // VERIFY
}

inline void Sequence::PrependBase_Unsafe(Nucleotide base)
{
    assert((base & 0x3) == base);
	PreparePrepend_Unsafe(1);
	SetBase_Unsafe(0, base);
}

inline void Sequence::PrependBase(SeqGraph *sgraph, SeqNode **node, Nucleotide base)
{
    if (length_ == MAX_BASES) {
        fprintf(stderr, "ERROR: maximum sequence length!");
        assert( false );
        exit(EXIT_FAILURE);
    }

    if (length_ == GetAllocatedBases()) {
        // need to grow node
        *node = GrowSeqNode(sgraph, *node, 1, true, false);
        (*node)->sequence.PrependBase_Unsafe(base);
        return;
    }
    PrependBase_Unsafe(base);
}

inline void Sequence::PopFront(uint16_t pop_bases)
{
    assert( pop_bases <= length_ );

    if (pop_bases == length_) { length_ = 0; return; }
    if (pop_bases == 0) { return; }

#ifdef VERIFY
    Sequence_StackAllocated memory1;
    Sequence *orig = new (&memory1) Sequence(*this, Sequence::MAX_BASES);
    assert( orig->get_length() == length_ );
#endif // VERIFY

    uintptr_t *words = reinterpret_cast<uintptr_t *>(str_); // FIXME alignment -- use memcpy?

    int sourceWordIndex = ComputeWordIndex<uintptr_t>(pop_bases);
    int limitSourceWordIndex = ComputeWordIndex<uintptr_t>(length_ - 1); // breaks if we don't check 'pop_bases == length' above
    assert( limitSourceWordIndex >= sourceWordIndex );
    assert( (limitSourceWordIndex * sizeof(uintptr_t)) < ComputeLengthInBytes(alloc_length_) );

    unsigned shamt = ComputeWordBitOffset<uintptr_t>(pop_bases);
    unsigned anti_shamt = (sizeof(uintptr_t) * 8) - shamt;

    unsigned targetWordIndex = 0;

    for ( ; sourceWordIndex < limitSourceWordIndex ; ++sourceWordIndex, ++targetWordIndex ) {
        words[targetWordIndex] = (words[sourceWordIndex] >> shamt);
        if(anti_shamt < 64) words[targetWordIndex] |= (words[sourceWordIndex+1] << anti_shamt);
    }
    words[targetWordIndex] = (words[sourceWordIndex] >> shamt);

    length_ -= pop_bases;

    /* // original code
    for (uint16_t i=0 ; i < length_ - bases ; ++i) {
        Nucleotide base = GetBase_Unsafe(i+bases);
        SetBase_Unsafe(i, base);
    }*/

#ifdef VERIFY
	for (int i=0 ; i < length_ ; ++i) {
		Nucleotide good_base = orig->GetBase_Unsafe(i+pop_bases);
        Nucleotide test_base = GetBase_Unsafe(i);
        assert( test_base == good_base );
    }
#endif // VERIFY
}

inline void Sequence::InitializeWithKmer_Unsafe(Kmer kmer, unsigned kmer_length)
{
    Kmer *words = reinterpret_cast<Kmer *>(str_); // FIXME alignment -- use memcpy?
    words[0] = kmer;
    length_ = kmer_length;
}

inline void Sequence::InitializeWithKmer(SeqGraph *sgraph, SeqNode **node, Kmer kmer, unsigned kmer_length)
{
    if (length_ + kmer_length > GetAllocatedBases()) {
        // need to grow node
        *node = GrowSeqNode(sgraph, *node, kmer_length, true, true);
        (*node)->sequence.InitializeWithKmer_Unsafe(kmer, kmer_length);
        return;
    }
    InitializeWithKmer_Unsafe(kmer, kmer_length);
}

inline void Sequence::InitializeWithString(char *string)
{
    length_ = 0;
    //printf("input string: %s\n", string);
    unsigned alloc_bases = Sequence::MAX_BASES; // FIXME
    for (unsigned i=0; i < alloc_bases && string[i] != '\0' ; ++i) {
		AppendBase_Unsafe(BASE_MAP[static_cast<int>(string[i])]);
    }
    // FIXME TODO: catch case where strings are too large for sequence
    // FIXME TODO: catch case where character is not recognized???
    //puts("parsed string: "); DebugPrint();
}

inline Kmer Sequence::GetKmer(unsigned offset, unsigned kmer_length) const
{
    assert( offset <= length_ - kmer_length );

    const Kmer *words = reinterpret_cast<const Kmer *>(str_); // FIXME alignment -- use memcpy?

    unsigned wordIndex0 = ComputeWordIndex<Kmer>(offset);
    unsigned wordIndex1 = ComputeWordIndex<Kmer>(offset + kmer_length - 1);

    unsigned wordBitOffset = ComputeWordBitOffset<Kmer>(offset);

    Kmer kmerA;
    if (wordIndex0 == wordIndex1) {
        kmerA = maskKmer((words[wordIndex0] >> wordBitOffset), kmer_length); // XXX: shift up instead of mask?
    } else {
        assert( wordIndex0 + 1 == wordIndex1 );
        Kmer word0 = words[wordIndex0] >> wordBitOffset;  // upper bits implicitly zero'd
        Kmer word1 = words[wordIndex1] << ((sizeof(Kmer) * 8) - wordBitOffset);
        kmerA = maskKmer((word0 | word1), kmer_length);
    }

#ifdef VERIFY
    Kmer kmerB = 0;
    for (unsigned i=0 ; i < kmer_length ; ++ i) {
        unsigned shamt = i * 2;
        kmerB |= static_cast<Kmer>(GetBase_Unsafe(i+offset)) << shamt;
    }
    assert( kmerA == kmerB );
#endif // VERIFY

    return kmerA;
}
	
inline Kmer Sequence::GetHeadKmer(unsigned kmer_length) const
{
    return GetKmer(0, kmer_length);
}

inline Kmer Sequence::GetTailKmer(unsigned kmer_length) const
{
    return GetKmer((length_ - kmer_length), kmer_length);
}

inline Nucleotide Sequence::GetHeadKmerRightmostBase(unsigned kmer_length) const
{
    return GetBase(kmer_length - 1);
}

inline Nucleotide Sequence::GetTailKmerLeftmostBase(unsigned kmer_length) const
{
    return GetBase(length_ - kmer_length);
}

inline void Sequence::AppendString_Unsafe(char *string, unsigned string_length)
{
    for (unsigned i=0; i < string_length; ++ i) {
        assert( string[i] != '\0' );
		AppendBase_Unsafe(BASE_MAP[static_cast<int>(string[i])]);
	}
}

inline void Sequence::AppendString(SeqGraph *sgraph, SeqNode **node, char *string, unsigned string_length)
{
    if (length_ + string_length > MAX_BASES) {
        fprintf(stderr, "ERROR: maximum sequence length!\n");
        assert( false );
        exit(EXIT_FAILURE);
    }

    if (length_ + string_length > GetAllocatedBases()) {
        // need to grow node
        *node = GrowSeqNode(sgraph, *node, string_length, false, true);
        (*node)->sequence.AppendString_Unsafe(string, string_length);
        return;
    }
    AppendString_Unsafe(string, string_length);
}

inline void Sequence::ToString(char *string, unsigned limit) const
{
    if (limit-1 < length_) {
        fprintf(stderr, "ERROR: string limit too short!\n");
        assert( false );
        exit(EXIT_FAILURE);
    }

	for (unsigned i=0; i < length_; ++ i) {
		string[i] = CHAR_BASE_MAP[ GetBase_Unsafe(i) ];
	}
	string[length_] = '\0';
}

// TODO OPT: this is really inefficient, but should work.
// concatenates two sequences, without duplicating kmer_length-1 overlap
inline void Sequence::Concatenate_Unsafe(const Sequence *next_seq, bool append, bool sense_changed, unsigned kmer_length)
{
	if (append) {
		if (sense_changed) {
			for (int i=next_seq->length_ - kmer_length ; i >= 0 ; -- i) {
				Nucleotide base = next_seq->GetBase_Unsafe(i);
				AppendBase_Unsafe(COMPLEMENT(base));
			}
		} else {
			for (uint16_t i=kmer_length-1 ; i < next_seq->length_ ; ++ i) {
				Nucleotide base = next_seq->GetBase_Unsafe(i);
				AppendBase_Unsafe(base);
			}
		}
	} else {
		// first, make space for the sequence to be prepended
        unsigned next_length = next_seq->length_ - (kmer_length-1);
		PreparePrepend_Unsafe(next_length);
		if (sense_changed) {
			for (uint16_t i=0 ; i < next_length ; ++ i) {
				Nucleotide base = next_seq->GetBase_Unsafe(kmer_length-1+i);
				SetBase_Unsafe((next_length-1)-i, COMPLEMENT(base));
			}
		} else {
			for (uint16_t i=0 ; i < next_length ; ++ i) {
				Nucleotide base = next_seq->GetBase_Unsafe(i);
				SetBase_Unsafe(i, base);
			}
		}
	}
}

inline void Sequence::Concatenate(SeqGraph *sgraph, SeqNode **node, const Sequence *next_seq, bool append, bool sense_changed, unsigned kmer_length)
{
	unsigned next_length = next_seq->length_ - (kmer_length-1);

    if (length_ + next_length > MAX_BASES) {
        fprintf(stderr, "ERROR: maximum sequence length!\n");
        assert( false );
        exit(EXIT_FAILURE);
    }

    if (length_ + next_length > GetAllocatedBases()) {
        // need to grow node
        *node = GrowSeqNode(sgraph, *node, next_length, (!append), (append));
        (*node)->sequence.Concatenate_Unsafe(next_seq, append, sense_changed, kmer_length);
        return;
    }
    Concatenate_Unsafe(next_seq, append, sense_changed, kmer_length);
}

// TODO OPT: this is really inefficient, but should work.
// reverse complement a sequence in-place
inline void Sequence::ReverseComplement(void)
{
    unsigned i, k;
    for (i=0, k=(length_ - 1) ; i <= k ; ++i, --k) {
        Nucleotide ibase = GetBase_Unsafe(i);
        Nucleotide kbase = GetBase_Unsafe(k);
        SetBase_Unsafe(i, COMPLEMENT(kbase));
        SetBase_Unsafe(k, COMPLEMENT(ibase));
    }
}

//
// file io
//

inline void Sequence::Load_BinaryFile(FILE* fp)
{
    size_t retval;

    retval = fread(&length_, sizeof(length_), 1, fp);
    if (retval != 1) {
        fprintf(stderr, "ERROR: fread() failed.\n");
        assert( false );
        exit(EXIT_FAILURE);
    }
        
    assert( length_ <= alloc_length_ );

    size_t str_size = GetLengthInBytes();
    retval = fread(str_, str_size, 1, fp);  // NOTE: assumes sufficient space allocated
    if (retval != 1) {
        fprintf(stderr, "ERROR: fread() failed.\n");
        assert( false );
        exit(EXIT_FAILURE);
    }
}

inline size_t Sequence::Save_BinaryFile(FILE* fp) const
{
    size_t retval;

    retval = fwrite(&length_, sizeof(length_), 1, fp);
    if (retval != 1) {
        fprintf(stderr, "ERROR: fwrite() failed.\n");
        assert( false );
        exit(EXIT_FAILURE);
    }

    size_t str_size = GetLengthInBytes();
    retval = fwrite(str_, str_size, 1, fp);
    if (retval != 1) {
        fprintf(stderr, "ERROR: fwrite() failed.\n");
        assert( false );
        exit(EXIT_FAILURE);
    }
    return (sizeof(length_) + str_size);
}
        
inline size_t Sequence::Load_Memory(void *mem)
{
    memcpy(&length_, mem, sizeof(length_));
    assert( length_ <= alloc_length_ );
    size_t str_size = GetLengthInBytes();
    memcpy(str_, static_cast<char*>(mem)+sizeof(length_), str_size);

    return ( sizeof(length_) + str_size );
}
        
inline size_t Sequence::Save_Memory(void *mem, size_t bytes_remaining) const
{
    size_t str_size = GetLengthInBytes();
    size_t seq_size = sizeof(length_) + str_size;
    if (bytes_remaining < seq_size) {
        return 0;
    }

    memcpy(mem, &length_, sizeof(length_));
    memcpy(static_cast<char*>(mem)+sizeof(length_), str_, str_size);
    return seq_size;
}
        
inline size_t Sequence::GetSerializedBytes(void) const
{
    size_t str_size = GetLengthInBytes();
    size_t seq_size = sizeof(length_) + str_size;
    return seq_size;
}

#endif // VELOUR_SEQUENCE_INLINED_H
