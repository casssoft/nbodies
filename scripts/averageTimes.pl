#!/usr/bin/perl

use 5.010000;
use warnings;
use strict;

my $total = 0.0;
my $count = 0;

while (my $line = <STDIN>) {
  chomp($line);
  if ($line =~ m/real\s+(\d+)m(\d+)\.(\d+)s$/) {
    $total += 60.0 * $1 + $2 + ($3 / 1000.0);
    $count += 1;
  }
}

my $average;

if ($count == 0) {
  print "Count is 0\n";
  $average = 0;
} else {
  $average = $total / $count;
}

print "Average: $average\n";
