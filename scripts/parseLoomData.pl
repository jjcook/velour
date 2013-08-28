#!/usr/bin/perl -w

#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

my $kmernodes = undef;
my $seqnodes = undef;
while(<>) {
  if( $_ =~ /^([0-9]+)\skmer\snodes\sbuilt/ )
  {
    $kmernodes = $1;
  }
  elsif( $_ =~ /^([0-9]+)\ssequence\snodes\sbuilt/ )
  {
    $seqnodes = $1;
    print "$kmernodes\t$seqnodes\n";
  }
}
