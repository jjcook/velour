#!/usr/bin/perl -w

#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

#108067 sequence nodes present before quilt file ssuis-21-hmetis-11-256/Partition-0256.quilt
#49482 sequence nodes loaded from quilt file ssuis-21-hmetis-11-256/Partition-0256.quilt
#157549 sequence nodes present after quilt f
#
#45325 sequence nodes dumped to pregraph

my $loaded_nodes = 0;
my $after_nodes = undef;
while(<>) {
  if( $_ =~ /^([0-9]+)\ssequence\snodes\sloaded/ )
  {
    $loaded_nodes += $1;
  }
  elsif( $_ =~ /^([0-9]+)\ssequence\snodes\spresent\sbefore/ )
  {
    print "$loaded_nodes\t$after_nodes\t$1\n" if $loaded_nodes > 0;
  }
  elsif( $_ =~ /^([0-9]+)\ssequence\snodes\spresent\safter/ )
  {
    $after_nodes = $1;
  }
  elsif( $_ =~ /^([0-9]+)\ssequence\snodes\sdumped/ )
  {
    print "$loaded_nodes\t$after_nodes\t$1\n";
  }
}
