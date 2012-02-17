#!/usr/bin/perl

open IN, "CONFIG_water" or die;

while (<IN>) {
   @line = split;
   $v1 = $line[0] / 10.0;
   $v2 = $line[1] / 10.0;
   $v3 = $line[2] / 10.0;
   $v4 = $line[3] / 1.0;
   print "$v0\t$v1\t$v2\t$v3\t$v4\n";
}

close IN
