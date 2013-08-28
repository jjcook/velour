//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// kmer.h
//

#ifndef VELOUR_KMER_H
#define VELOUR_KMER_H

typedef uint64_t Kmer;

static inline Nucleotide getNucleotide(Kmer k, unsigned offset)
{
    Nucleotide n = (k >> (offset << 1)) & 0x3;
    return n;
}

static inline Kmer maskKmer(Kmer k, unsigned kmer_length)
{
    Kmer mask = ((((Kmer)1) << (kmer_length << 1)) - 1);
    return (k & mask);
}

// Reverses a word into its Watson-Crick reverse complement
static inline Kmer reverseComplement(Kmer k, unsigned kmer_length)
{
    //k = ((k >>  1) & 0x5555555555555555ULL) | ((k <<  1) & 0xaaaaaaaaaaaaaaaaULL); // don't reverse the base itself!
    k = ((k >>  2) & 0x3333333333333333ULL) | ((k <<  2) & 0xccccccccccccccccULL);
    k = ((k >>  4) & 0x0f0f0f0f0f0f0f0fULL) | ((k <<  4) & 0xf0f0f0f0f0f0f0f0ULL);
    k = ((k >>  8) & 0x00ff00ff00ff00ffULL) | ((k <<  8) & 0xff00ff00ff00ff00ULL);
    k = ((k >> 16) & 0x0000ffff0000ffffULL) | ((k << 16) & 0xffff0000ffff0000ULL);
    k = ((k >> 32) & 0x00000000ffffffffULL) | ((k << 32) & 0xffffffff00000000ULL);

    k >>= 64 - (2 * kmer_length);

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
