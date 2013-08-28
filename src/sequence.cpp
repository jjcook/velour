//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// sequence.cpp
//

#include "types.h"
#include "sequence.h"
        
// TODO OPT: use a larger table to make more efficient
void Sequence::PrintToFile(FILE *file) const
{
    for (unsigned i = 0 ; i < length_ ; ++ i) {
        Nucleotide base = GetBase_Unsafe(i);
        fputc(CHAR_BASE_MAP[base], file);
    }
}

void Sequence::DebugPrint(void) const
{
    PrintToFile(stdout);
    fputc('\n', stdout);
    fflush(stdout);
}

void Sequence::Print(void) const
{
    PrintToFile(stdout);
}
