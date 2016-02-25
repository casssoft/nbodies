#!/usr/bin/perl

use 5.010000;
use warnings;
use strict;

my @output = (0.00248375, 0.00330892, -0.00866389, -0.596418, 0.0238128, -0.00878246);
my $elapsedYears = 6.52309e+06;
my $epsilon = 0.00001;
my $outputCounter = 0;

for (my $counter = 0; $counter < 3; ++$counter) {
  <STDIN>;
}

for (my $counter = 0; $counter < 3; ++$counter) {
  my $position = <STDIN>;
  my $expectedValue = $output[$outputCounter++];
  die("Expected $expectedValue, received $position") if ($position + $epsilon < $expectedValue 
    || $position - $epsilon > $expectedValue);
}

<STDIN>;

for (my $counter = 0; $counter < 3; ++$counter) {
  my $position = <STDIN>;
  my $expectedValue = $output[$outputCounter++];
  die("Expected $expectedValue, received $position") if ($position + $epsilon < $expectedValue 
    || $position - $epsilon > $expectedValue);
}
