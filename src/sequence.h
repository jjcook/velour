//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// sequence.h
//
//   variable length nucleotide sequence
//

// TODO: don't Init

#ifndef VELOUR_SEQUENCE_H
#define VELOUR_SEQUENCE_H

class SeqGraph;
class SeqNode;

class Sequence
{
    public:
        static const uint16_t MAX_BASES = UINT16_MAX; // NOTE: mutating members should check this
        static const size_t MIN_BYTES = STATIC_ROUND_UP_POWER_OF_2< sizeof(Kmer) >::RESULT;
        static const size_t MAX_BYTES = (MAX_BASES + 3) >> 2;

#ifdef NDEBUG
        Sequence(bool) /* : length_(MAX_BASES) */ {}   // used for stack allocation
        Sequence(const uint16_t alloc_bases) : length_(0) { /*Init(alloc_bases);*/ }
        Sequence(const Sequence& s, const uint16_t alloc_bases) :
            length_(s.length_) { /*Init(alloc_bases);*/ memcpy(str_, s.str_, s.GetLengthInBytes()); }
#else
        Sequence(bool) : length_(MAX_BASES), alloc_length_(MAX_BASES) {}    // used for stack allocation
        Sequence(const uint16_t alloc_bases) : length_(0), alloc_length_(alloc_bases) { /*Init(alloc_bases);*/ }
        Sequence(const Sequence& s, const uint16_t alloc_bases) :
            length_(s.length_), alloc_length_(alloc_bases)
            { /*Init(alloc_bases);*/ assert(alloc_length_ >= s.length_); memcpy(str_, s.str_, s.GetLengthInBytes()); }
#endif

        uint16_t get_length() const;
        void zero_length(); // NOTE: used by memory allocator to pseudo-init when marking unallocated node dead

        bool CanAddBases(unsigned count) const;

        static size_t ComputeLengthInBytes(uint16_t bases);
        static size_t ComputeAllocatedBytes(uint16_t bases);
        static uint16_t ComputeAllocatedBases(uint16_t bases);

        size_t GetLengthInBytes() const;
        size_t GetAllocatedBytes() const;
        uint16_t GetAllocatedBases() const;

        Nucleotide GetBase(const uint16_t index) const;

        void AppendBase(SeqGraph *sgraph, SeqNode **node, Nucleotide base);
        void AppendBase_Unsafe(Nucleotide base);

        void PrependBase(SeqGraph *sgraph, SeqNode **node, Nucleotide base);
        void PrependBase_Unsafe(Nucleotide base);

        void PopFront(uint16_t bases);

        void InitializeWithKmer(SeqGraph *sgraph, SeqNode **node, Kmer kmer, unsigned kmer_length);
        void InitializeWithKmer_Unsafe(Kmer kmer, unsigned kmer_length);

        void InitializeWithString(char *string);

        Kmer GetKmer(unsigned offset, unsigned kmer_length) const;
        Kmer GetHeadKmer(unsigned kmer_length) const;
        Kmer GetTailKmer(unsigned kmer_length) const;

        Nucleotide GetHeadKmerRightmostBase(unsigned kmer_length) const;
        Nucleotide GetTailKmerLeftmostBase(unsigned kmer_length) const;

        void AppendString(SeqGraph *sgraph, SeqNode **node, char *string, unsigned string_length);
        void AppendString_Unsafe(char *string, unsigned string_length);
        void ToString(char *string, unsigned limit) const;

        void Concatenate(SeqGraph *sgraph, SeqNode **node, const Sequence *next_seq, bool append, bool sense_changed, unsigned kmer_length);
        void Concatenate_Unsafe(const Sequence *next_seq, bool append, bool sense_changed, unsigned kmer_length);

        void ReverseComplement(void);

        void Load_BinaryFile(FILE* fp); // NOTE: assumes MAX_BASES space available
        size_t Save_BinaryFile(FILE* fp) const;

        size_t Load_Memory(void *mem); // NOTE: assumes MAX_BASES space available
        size_t Save_Memory(void *mem, size_t bytes_remaining) const;

        size_t GetSerializedBytes(void) const;

        void PrintToFile(FILE *file) const;
        void Print(void) const;
        void DebugPrint(void) const;


    private:
        void Init(uint16_t alloc_bases);
        Nucleotide GetBase_Unsafe(const uint16_t index) const;
        void SetBase_Unsafe(const uint16_t index, Nucleotide base);
        void PreparePrepend_Unsafe(uint16_t insert_length);

        // private declarations with no definition to disable compiler generated
        Sequence();
        void operator=(const Sequence&);

#ifndef NDEBUG
        friend SeqNode * GrowSeqNode(SeqGraph * sgraph, SeqNode * oldptr, unsigned bases, bool skip_head_update, bool skip_tail_update); // to support reallocation update of alloc_length_
#endif

    private:
        uint16_t length_;
#ifndef NDEBUG
        uint16_t alloc_length_;
#endif
        byte_t   str_[0];
} PACKED;

//
// use this class to stack allocate the variable length Sequence of maximal length
//
class Sequence_StackAllocated   // use instead of 'alloca()' since inlining it is dangerous
{
    public:
        Sequence_StackAllocated() : __dummy_sequence(false) {}

    private:
        Sequence __dummy_sequence;
        byte_t   y_[Sequence::MAX_BYTES];
} PACKED;

//
// NOTE: inlined Sequence functions are in sequence_inlined.h
//

template<typename T>
void sequence_process_kmers(Sequence *seq, unsigned kmer_length, T& functor)
{
    unsigned seq_length = seq->get_length();
    if (kmer_length >= seq_length) { return; }

    Kmer kmer = seq->GetHeadKmer(kmer_length);
    Kmer rc_kmer = reverseComplement(kmer, kmer_length);
    bool sense_reversed = rc_kmer < kmer;
    Kmer canonical_kmer = ( sense_reversed ? rc_kmer : kmer );

    functor(kmer, canonical_kmer, sense_reversed, 0, 0, kmer_length, true, (kmer_length == seq_length)); // is first

    unsigned double_kmer_length = kmer_length << 1;
    Kmer mask = (((Kmer)1) << double_kmer_length) - 1;

    for(unsigned i=kmer_length ; i < (seq_length - 1) ; ++i) {
        Nucleotide last_base = kmer & 0x3;
        Nucleotide base = seq->GetBase(i);
        kmer = KMER_APPEND(kmer, base, double_kmer_length);
        rc_kmer = KMER_PREPEND(rc_kmer, COMPLEMENT(base), double_kmer_length, mask);
        bool sense_reversed = rc_kmer < kmer;
        Kmer canonical_kmer = ( sense_reversed ? rc_kmer : kmer );

        functor(kmer, canonical_kmer, sense_reversed, base, last_base, kmer_length, false, false);
    }

    if (kmer_length < seq_length) {
        Nucleotide last_base = kmer & 0x3;
        Nucleotide base = seq->GetBase(seq_length-1);
        kmer = KMER_APPEND(kmer, base, double_kmer_length);
        rc_kmer = KMER_PREPEND(rc_kmer, COMPLEMENT(base), double_kmer_length, mask);
        bool sense_reversed = rc_kmer < kmer;
        Kmer canonical_kmer = ( sense_reversed ? rc_kmer : kmer );

        functor(kmer, canonical_kmer, sense_reversed, base, last_base, kmer_length, false, true); // is last
    }
}

#endif // VELOUR_SEQUENCE_H
