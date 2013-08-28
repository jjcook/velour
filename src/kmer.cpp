//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// kmer.cpp
//
// kmer utility functions
//

#include "types.h"
#include "kmer.h"

void gdb_print_kmer(Kmer kmer, unsigned kmer_length)
{
    fprint_kmer(kmer, kmer_length, stdout);
    putc('\n', stdout);
    fflush(stdout);
}

void fprint_kmer(Kmer kmer, unsigned kmer_length, FILE *file)
{
    for (unsigned i = 0 ; i < kmer_length ; ++ i) {
        Nucleotide base = (kmer >> (i << 1)) & 0x3;
        fputc(CHAR_BASE_MAP[base], file);
    }
}

