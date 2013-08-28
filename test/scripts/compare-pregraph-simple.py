#!/usr/bin/python

#
# compare-pregraph-simple.py
#
#   compare two pregraphs to equality modulo reverse complementation
#
#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

import sys, string, random, os

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
        else:
            print "ERROR: invalid base character %c.  exiting..." % (c)
            sys.exit(1)
        output = c2 + output
    return output

class Sequence:
    def __init__(self, i, string):
        self.i = i
        self.string = string
        self.length = len(string)
    def reverse_complement(self):
        return Sequence(self.i, rc(self.string))
    def equals(self, s):
        return self.string == s.string
    def print_me(self):
        print self.i, self.length, self.string

def get_sequence(file):
    lines = open(file).readlines()
    tokens = lines[0].split()
    num_strings = string.atoi(tokens[0], 0)
    array = []
    # print len(lines)
    for i in range(1, num_strings+1):
        array.append(Sequence(i, lines[i*2].strip()))
    return array

def find_and_remove(array, s):
    for i in range(len(array)):
        s1 = array[i]
        if s1.equals(s):
            # print array
            del array[i]
            # print array
            ##print "matched %d with %d" % (s1.i, s.i)
            return True
    return False

#if len(sys.argv) == 2:
#    print rc(sys.argv[1])
#    sys.exit(0)

if len(sys.argv) != 3:
    print "ERROR: Usage: %s <file1.pregraph> <file2.pregraph>" % (sys.argv[0])
    sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2]

#print "Matching %s with %s" % (file1, file2)

sequence_map = {}
sequences1 = get_sequence(file1)
for s in sequences1:
    l = s.length
    if not sequence_map.has_key(l):
        sequence_map[l] = []
    sequence_map[l].append(s)
    
sequences2 = get_sequence(file2)
#for s in sequences2:
#    s.print_me()
    
remaining_sequences = []
for s in sequences2:
    l = s.length
    if sequence_map.has_key(l):
        if find_and_remove(sequence_map[l], s) or find_and_remove(sequence_map[l], s.reverse_complement()):
            continue
    remaining_sequences.append(s)

if len(remaining_sequences) != 0:
    if os.getenv("BESILENT") is None:
        print "FAIL: Unmatched sequence(s) in file %s" % (file2)
    sys.exit(2)

for k in sequence_map.keys():
    if len(sequence_map[k]) != 0:
        if os.getenv("BESILENT") is None:
            print "FAIL: Unmatched sequence(s) in file %s" % (file1)
        sys.exit(2)

if os.getenv("BESILENT") is None:
    print "PASS: Pregraph files matched."
sys.exit(0)

# NOTE: the code below can be used to print what sequences did not match
#index_array = []
#for k in sequence_map.keys():
#    for s in sequence_map[k]:
#        index_array.append(s.i)
#index_array.sort()
#print "sequences left in " + file1
#for i in index_array:
#    sequences1[i-1].print_me()
#
#print "sequences left in " + file2
#for s in remaining_sequences:
#    s.print_me()
#
#for i in index_array:
#    tries = 5
#    s = sequences1[i-1]
#    l = len(s.string)
#    while tries > 0:
#        tries -= 1
#        index = random.randrange(3, l-8-3)
#        sub_sequence = s.string[index:index+8]
#        print "looking for subsequence %s in sequence %d" % (sub_sequence, i)
#        rc_sub_sequence = rc(sub_sequence)
#        for s2 in remaining_sequences:
#            if (s2.string.find(sub_sequence) != -1 or s2.string.find(rc_sub_sequence) != -1):
#                print "    found match in sequence %d" % s2.i
#                tries = 0
 
