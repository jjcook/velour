#!/usr/bin/perl -w

#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

my $prenodes = 0;
my $nodes = 0;

while(<>) {
  if( $_ =~ /^([0-9]+)\suseful\spre-nodes/ )
  {
    $prenodes += $1;
  }
  elsif( $_ =~ /^([0-9]+)\suseful\snodes/ )
  {
    $nodes += $1;
  }
}
my $ratio = $nodes / $prenodes;
printf "\t%d\t%d\t%0.2f%%\n", $prenodes, $nodes, 100*$ratio;
