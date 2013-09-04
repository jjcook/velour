//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// kmer.h
//

#ifndef VELOUR_KMER_H
#define VELOUR_KMER_H


#if MAXKMERLENGTH <= 31

typedef uint64_t Kmer;

static inline Nucleotide getNucleotide(Kmer k, unsigned index)
{
    Nucleotide n = (k >> (index << 1)) & 0x3;
    return n;
}

static inline Kmer setNucleotide(Kmer k, unsigned index, Nucleotide n)
{
    Kmer mask = ~( (Kmer)0x3 << (index << 1) );
    Kmer insert = (Kmer)n << (index << 1);
    return (k & mask) | insert;
}

#  define convertKmerToUint64(k) (k)

#else
#  define LARGE_KMERS
class Kmer {
  public:
    static const unsigned kmer_bytes = MAXKMERLENGTH % 4 ? (MAXKMERLENGTH / 4) + 1 : (MAXKMERLENGTH / 4);
    static const unsigned kmer_words = kmer_bytes % 8 ? (kmer_bytes / 8) + 1 : (kmer_bytes / 8);

    Kmer() {}
    Kmer(const Kmer& k) { for(unsigned i=0; i < kmer_words; ++i) word[i] = k.word[i]; }
    Kmer(uint64_t j) { word[0] = j; for(unsigned i=1; i < kmer_words; ++i) word[i] = 0; }
  public:
    uint64_t word[kmer_words];

    void createMask(const unsigned bits)
    {
        unsigned words = bits >> 6; // div 64
        unsigned rem = bits & 0x3f; // mod 64
        unsigned i;
        for(i=0; i < words; ++i)
            word[i] = UINT64_MAX;
        word[i] = (((uint64_t)1) << rem) - 1;
        for(++i; i < kmer_words; ++i)
            word[i] = 0;
    }
} PACKED;

static inline uint64_t convertKmerToUint64(const Kmer& k)
{
#ifdef VERIFY
    for(unsigned i=Kmer::kmer_words-1; i >= 1; --i)
        assert( k.word[i] == 0 );
#endif
    return k.word[0];
}

static inline Nucleotide getNucleotide(Kmer& k, const unsigned index)
{
    unsigned offset_words = (index << 1) >> 6; // div 64
    unsigned offset_rem = (index << 1) & 0x3f; // mod 64

    return (k.word[offset_words] >> offset_rem) & 0x3;
}

static inline Kmer& setNucleotide(Kmer& k, const unsigned index, Nucleotide n)
{
    unsigned offset_words = (index << 1) >> 6; // div 64
    unsigned offset_rem = (index << 1) & 0x3f; // mod 64

    k.word[offset_words] &= ~((uint64_t)0x3 << offset_rem);
    k.word[offset_words] |= (uint64_t)n << offset_rem;

    return k;
}

/*static inline operator= (Kmer& lhs, const Kmer& rhs)
{
    for(int i=0; i < Kmer::kmer_words; ++i)
        lhs.word[i] = rhs.word[i];
}*/

static inline bool operator==(const Kmer& lhs, const Kmer& rhs)
{
    for(int i=Kmer::kmer_words-1; i >= 0; --i)
    {
        if (lhs.word[i] != rhs.word[i])
            return false;
    }
    return true;
}
static inline bool operator!=(const Kmer& lhs, const Kmer& rhs) { return !operator==(lhs, rhs); }
static inline bool operator< (const Kmer& lhs, const Kmer& rhs)
{
    for(int i=Kmer::kmer_words-1; i >= 0; --i)
    {
        if( lhs.word[i] < rhs.word[i] )
            return true;
        else if(lhs.word[i] > rhs.word[i])
            return false;
    }
    return false;
}
static inline bool operator> (const Kmer& lhs, const Kmer& rhs) { return  operator< (rhs, lhs);}
static inline bool operator<=(const Kmer& lhs, const Kmer& rhs) { return  !operator> (lhs, rhs);}
static inline bool operator>=(const Kmer& lhs, const Kmer& rhs) { return  !operator< (lhs, rhs);}

static inline bool operator< (const Kmer& lhs, const uint64_t rhs)
{
    return lhs.word[0] < rhs;
}

static inline Kmer operator<< (Kmer lhs, const unsigned shift)
{
    unsigned shift_words = shift >> 6; // div 64
    unsigned shift_rem = shift & 0x3f; // mod 64

    int i;
    for(i=Kmer::kmer_words-1; i - shift_words >= 1; --i)
    {
        lhs.word[i] = lhs.word[i-shift_words] << shift_rem;
        if(shift_rem > 0) lhs.word[i] |= lhs.word[i-shift_words-1] >> (64 - shift_rem);
    }
    lhs.word[i] = lhs.word[i-shift_words] << shift_rem;
    for(--i; i >= 0; --i)
        lhs.word[i] = 0;

    return lhs;
}

static inline Kmer operator>> (Kmer lhs, const unsigned shift)
{
    unsigned shift_words = shift >> 6; // div 64
    unsigned shift_rem = shift & 0x3f; // mod 64

    unsigned i;
    for(i=0; i + shift_words + 1 < Kmer::kmer_words; ++i)
    {
        lhs.word[i] = lhs.word[i+shift_words] >> shift_rem;
        if(shift_rem > 0) lhs.word[i] |= lhs.word[i+shift_words+1] << (64 - shift_rem);
    }
    lhs.word[i] = lhs.word[i+shift_words] >> shift_rem;
    for(++i; i < Kmer::kmer_words; ++i)
        lhs.word[i] = 0;

    return lhs;
}

static inline Kmer operator& (Kmer lhs, const Kmer& rhs)
{
    for(int i=Kmer::kmer_words-1; i >= 0; --i)
    {
        lhs.word[i] &= rhs.word[i];
    }
    return lhs;
}

static inline uint64_t operator& (Kmer lhs, const uint64_t mask)
{
    return lhs.word[0] & mask;
}

static inline Kmer operator| (Kmer lhs, const Kmer& rhs)
{
    for(int i=Kmer::kmer_words-1; i >= 0; --i)
    {
        lhs.word[i] |= rhs.word[i];
    }
    return lhs;
}

static inline Kmer operator^ (Kmer lhs, const Kmer& rhs)
{
    for(int i=Kmer::kmer_words-1; i >= 0; --i)
    {
        lhs.word[i] ^= rhs.word[i];
    }
    return lhs;
}

static inline Kmer operator~ (Kmer lhs)
{
    for(int i=Kmer::kmer_words-1; i >= 0; --i)
    {
        lhs.word[i] = ~lhs.word[i];
    }
    return lhs;
}
#endif // MAXKMERLENGTH


static inline Kmer maskKmer(Kmer k, unsigned kmer_length)
{
#ifdef LARGE_KMERS
    Kmer mask;
    mask.createMask(kmer_length << 1);
#else
    Kmer mask = ((((Kmer)1) << (kmer_length << 1)) - 1);
#endif
    return (k & mask);
}

// Reverses a word into its Watson-Crick reverse complement
static inline Kmer reverseComplement(Kmer k, unsigned kmer_length)
{
#ifdef LARGE_KMERS
    // XXX OPT: use bit-shuffling like as for the word-sized kmers
    for(int i=0; i < (kmer_length >> 1); ++i)
    {
        unsigned lhs_idx = kmer_length - 1 - i;
        unsigned rhs_idx = i;
        Nucleotide lhs = getNucleotide(k, lhs_idx);
        Nucleotide rhs = getNucleotide(k, rhs_idx);
        setNucleotide(k, lhs_idx, rhs);
        setNucleotide(k, rhs_idx, lhs);
    }
#else
    //k = ((k >>  1) & 0x5555555555555555ULL) | ((k <<  1) & 0xaaaaaaaaaaaaaaaaULL); // don't reverse the base itself!
    k = ((k >>  2) & 0x3333333333333333ULL) | ((k <<  2) & 0xccccccccccccccccULL);
    k = ((k >>  4) & 0x0f0f0f0f0f0f0f0fULL) | ((k <<  4) & 0xf0f0f0f0f0f0f0f0ULL);
    k = ((k >>  8) & 0x00ff00ff00ff00ffULL) | ((k <<  8) & 0xff00ff00ff00ff00ULL);
    k = ((k >> 16) & 0x0000ffff0000ffffULL) | ((k << 16) & 0xffff0000ffff0000ULL);
    k = ((k >> 32) & 0x00000000ffffffffULL) | ((k << 32) & 0xffffffff00000000ULL);

    k >>= 64 - (2 * kmer_length);
#endif

    Kmer rc = maskKmer(~k, kmer_length);
    return rc;
}

static inline Kmer canonicalKmer(Kmer k, unsigned kmer_length)
{
    Kmer rc = reverseComplement(k, kmer_length);
    return (k < rc ? k : rc);
}

//
// non-inline functions
//

void gdb_print_kmer(Kmer kmer, unsigned length);
void fprint_kmer(Kmer kmer, unsigned length, FILE *);

static inline void print_kmer(Kmer kmer, unsigned kmer_length)
{
    fprint_kmer(kmer, kmer_length, stdout);
}

#endif // VELOUR_KMER_H
