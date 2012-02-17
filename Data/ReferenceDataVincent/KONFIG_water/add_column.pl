#!/usr/bin/perl

open IN, "CONFIG_water.exact" or die;

while (<IN>) {
   @line = split;
   $v1 = $line[0];
   $v2 = $line[1];
   $v3 = $line[2];
   print "1.0\t$v1\t$v2\t$v3\n";
}

close IN
