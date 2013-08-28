#!/usr/bin/python
#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

import sys, string, random

def rc(line):
    output = ""
    for c in line:
        if c == "A":
            c2 = "T"
        elif c == "C":
            c2 = "G"
        elif c == "G":
            c2 = "C"
        elif c == "T":
            c2 = "A"
        elif c == "N":
            c2 = "N"
        else:
            sys.exit(1)
        output = c2 + output
    return output

def toint(line):
    output = 0
    for c in line:
        if c == "A":
            output = (output << 2) | 0x0
        elif c == "C":
            output = (output << 2) | 0x1
        elif c == "G":
            output = (output << 2) | 0x2
        elif c == "T":
            output = (output << 2) | 0x3
        else:
            sys.exit(1)
    return output

if len(sys.argv) == 2:
    #print "rc", rc(sys.argv[1])
    input_int = toint(sys.argv[1])
    input_hash = (input_int ^ (input_int >> 17)) & 0x0000000000ffffffL
    rc_int = toint(rc(sys.argv[1]))
    rc_hash = (rc_int ^ (rc_int >> 17)) & 0x0000000000ffffffL
    print "input: ", sys.argv[1], " ", hex(input_int), " ", hex(input_hash)
    print "   rc: ", rc(sys.argv[1]), " ", hex(rc_int), " ", hex(rc_hash)
    sys.exit(0)

