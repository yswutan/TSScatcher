#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "Usage: $0 <paraclu.input> <file.bed>\n";
my $infile = shift or die $usage;
my $bed = shift or die $usage;

my %store = ();

open(IN,'<',$infile) || die "Could not open $infile: $!\n";
while(<IN>){
   chomp;
   #chr1    -       16605   27
   my ($chr, $strand, $start, $count) = split(/\t/);
   $store{$strand}{$chr}{$start} = $count;
}
close(IN);

open(IN,'<',$bed) || die "Could not open $bed: $!\n";
while(<IN>){
   chomp;
   my ($chr, $start, $end, $junk, $count, $strand) = split('\t');
   my $max = 0;
   my @pos = ();
   for ($start .. $end){
      if (exists $store{$strand}{$chr}{$_}){
         my $count = $store{$strand}{$chr}{$_};
         if ($count > $max){
            $max = $count;
            @pos = ();
            $pos[0] = $_;
         } elsif ($count == $max){
            push(@pos, $_);
         }
      }
   }
   my $peak = join(',', @pos);
   print join("\t", $_, $max, $peak),"\n";
}
close(IN);

exit(0);
