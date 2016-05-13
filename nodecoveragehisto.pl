#!/usr/bin/perl -w

###
### Compute node coverage histograms
###

use warnings;
use strict;

my %cov_histo;
my %lwcc_histo;

while(<>) {
  if( $_ =~ /^([0-9]+)\t([0-9]+.[0-9]+)/ ) {
    my $contig_kmers = $1;
    my $contig_cov = $2;

    $cov_histo{ int($2 * 10.0) }++;
    $lwcc_histo{ int($contig_kmers * $contig_cov * 10.0) }++;
  }
  elsif( $_ =~ /^([0-9]+)\t([0-9]+)/ ) {
    my $contig_kmers = $1;
    my $contig_cov = $2;

    $cov_histo{ int($2 * 10.0) }++;
    $lwcc_histo{ int($contig_kmers * $contig_cov * 10.0) }++;
  }
}

foreach my $bucket ( sort { $a <=> $b } keys %cov_histo ) {
    printf "COV[%04.1f]:  %16d\n", $bucket / 10.0, $cov_histo{$bucket};
}

foreach my $bucket ( sort { $a <=> $b } keys %lwcc_histo ) {
    printf "LWCC[%04.1f]: %16d\n", $bucket / 10.0, $lwcc_histo{$bucket};
}
