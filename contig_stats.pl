#!/usr/bin/perl
#
# contig_stats.pl -- FASTA or PreGraph file formats
#
#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

#### Parameters ####

$ID = undef;
$MIN_SIZE = 0;
$TOP_N = undef;
$KMERSIZE = 1;

#### Statistics values -- iteratively updated

$num_contigs = 0;
$total_bases = 0;
$shortest_contig = undef;
$longest_contig = undef;

#### Statistics values -- computed 

$average_length = undef;
$median = undef;

$n50_length = undef;
$n50_count = undef;
$n50_percent = undef;

$n80_length = undef;
$n80_count = undef;
$n80_percent = undef;

#### Internal variables

@contig_lengths;

### get command line options
use Getopt::Std;
$::opt_i = $::opt_m = $::opt_t = $::opt_k = undef;
getopts('i:m:t:k:');

#$ID         = $::opt_i  if (defined $::opt_i);
$MIN_SIZE   = $::opt_m  if (defined $::opt_m);
#$TOP_N      = $::opt_t  if (defined $::opt_t);
$KMERSIZE   = $::opt_k  if (defined $::opt_k);

## go

$MIN_SIZE += ($KMERSIZE - 1);

$fasta_input_file = shift @ARGV;
$fasta_input_file = '-' if (!$fasta_input_file);
open(FASTAIN, $fasta_input_file) || die("Can't open fasta file: '$fasta_input_file'\n");

$ID = $fasta_input_file if (!$ID);

$contig = undef;
while($line = <FASTAIN>)
{
  chomp $line;
  if($line =~ /^>/)
  {
    process_contig($contig) if ($contig);
    $contig = undef;
  }
  elsif($line =~ /^NODE/)       # PREGRAPH
  {
    process_contig($contig) if ($contig);
    $contig = undef;
  }
  elsif($line =~ /^[0-9]+/)      # PREGRAPH
  {
    next;
  }
  else
  {
    $contig .= $line;
  }
}
if($contig) # handle end of file
{
  process_contig($contig);
}
computed_stats();
print_all_stats();
exit;


## helper routines

sub min
{
  my ($a, $b) = @_;
  return $a if (not defined $b);
  return $b if (not defined $a);
  if($a <= $b) { return $a; } else { return $b; }
}

sub max
{
  my ($a, $b) = @_;
  return $a if (not defined $b);
  return $b if (not defined $a);
  if($a >= $b) { return $a; } else { return $b; }
}

sub process_contig
{
  my ($contig) = @_;

  $len = length($contig);

  if($len < $MIN_SIZE) { return; }

  $len -= ($KMERSIZE - 1);

  $num_contigs++;
  $total_bases += $len;
  $shortest_contig = min($len, $shortest_contig);
  $longest_contig = max($len, $longest_contig);

  push @contig_lengths, $len;
}

sub computed_stats
{
  $average_length = $total_bases / $num_contigs;

  # compute median
  my @sorted_lengths = sort { $a <=> $b } @contig_lengths;
  if($num_contigs & 1)
  {
    $median = $sorted_lengths[$num_contigs >> 1];
  }
  else
  {
    my $sum = $sorted_lengths[$num_contigs >> 1] + $sorted_lengths[($num_contigs >> 1)-1];
    $median = ($sum & 1) ? sprintf("%.1f", $sum/2) : $sum >> 1;
  }

  ($n50_length, $n50_count) = compute_N(0.50);
  $n50_percent = $n50_count / $num_contigs;

  ($n80_length, $n80_count) = compute_N(0.80);
  $n80_percent = $n80_count / $num_contigs;
}

sub compute_N
{
  my ($percentage) = @_;

  my $nXX_length = $nXX_count = 0;

  my $base_coverage = int(($total_bases * $percentage) + 0.5);

  my @sorted_lengths = sort { $a <=> $b } @contig_lengths;

  while($base_coverage > 0)
  {
    $len = pop @sorted_lengths;
    $base_coverage -= $len;
    
    $nXX_count++;
    $nXX_length = $len;
  }	
					
  return ($nXX_length, $nXX_count);
}

sub print_all_stats
{
  #print "$ID\n";

  printf("num_contigs     = %10d\n", $num_contigs);
  printf("total_bases     = %10d\n", $total_bases);
  printf("shortest_contig = %10d\n", $shortest_contig);
  printf("longest_contig  = %10d\n", $longest_contig);

  printf("average_length  = %10d\n", $average_length);
  printf("median_length   = %10d\n", $median);

  printf("n50_length      = %10d\n", $n50_length);
  printf("n50_count       = %10d\n", $n50_count);
  printf("n50_percent     = %10d\n", $n50_percent*100);

  printf("n80_length      = %10d\n", $n80_length);
  printf("n80_count       = %10d\n", $n80_count);
  printf("n80_percent     = %10d\n", $n80_percent*100);
}

