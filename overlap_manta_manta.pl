#! /usr/bin/perl

use strict;
use Data::Dumper;

my $manta = shift; # reference
my $manta2 = shift; # sample
my $dist = shift;

my %svs;
open MANTA, $manta;
while(<MANTA>) {
    chomp;
    next if $_ =~ /^#/;
    my ($chr, $pos1, $id, $ref, $alt, $qual, $filter, $info, $format, $genotype) = split/\t/, $_;
    my $chr2 = $chr;
    my $pos2 = 0;
    ($chr2, $pos2) = ($1, $2) if $alt =~ /(\w+):(\d+)/;
    $pos2 = $1 if $info =~ /END=(\d+);/;
    $pos2 = $1 if $info =~ /;END=(\d+)/;

    my ( $pos1a, $pos1b ) = ( $pos1, $pos1 );
    $pos1a -= $1 if $info =~ /CIPOS=-*(\d+),/;
    $pos1b += $1 if $info =~ /CIPOS=-*\d+,(\d+)/;
    
    my ( $pos2a, $pos2b ) = ( $pos2, $pos2 );
    $pos2a -= $1 if $info =~ /CIEND=-*(\d+),/;
    $pos2b += $1 if $info =~ /CIEND=-*\d+,(\d+)/;
    
    my $ori = "XX";
    if ($id =~ /BND/) {
	$ori = "TH" if $alt =~ /^\w+\[/;
        $ori = "HT" if $alt =~ /^\]/;
        $ori = "TT" if $alt =~ /^\w+\]/;
	$ori = "HH" if $alt =~ /^\[/;
    } else {
        $ori = "TH" if $id =~ /DEL/;
        $ori = "TH" if $id =~ /INS/;
        $ori = "HT" if $id =~ /DUP/;
        $ori = "TT" if $id =~ /INV/ and $info =~ /INV3/;
        $ori = "HH" if $id =~ /INV/ and $info =~ /INV5/;
    }
    
    my $type = $1 if $info =~ /SVTYPE=(\w+)/;
    $type = "BND" if $type eq "INV";
    $type = "BND" if $type eq "TRA";

    $type = "X";
    
    if ($pos1>$pos2 and $chr eq $chr2) {
      ($pos1a, $pos1b, $pos2a, $pos2b) = ($pos2a, $pos2b, $pos1a, $pos1b);
      $ori = reverse($ori);
    }
    if ($chr2 < $chr) {
      ($chr, $chr2, $pos1a, $pos1b, $pos2a, $pos2b) = ($chr2, $chr, $pos2a, $pos2b, $pos1a, $pos1b);
      $ori=reverse($ori);
    }

    $svs{$type}{$chr}{$chr2}{$ori}{$pos1a}{$pos1b}{$pos2a}{$pos2b} = $id;
}
close MANTA;


open MANTA2, $manta2;
while(<MANTA2>) {
  chomp;
  print $_ ."\n" if $_ =~ /^#/;
  next if $_ =~ /^#/;
  my @line = split/\t/, $_;
  my ($chr, $pos1, $id, $ref, $alt, $qual, $filter, $info, $format, $genotype) = split/\t/, $_;
  my $chr2 = $chr;
  my $pos2 = 0;
  ($chr2, $pos2) = ($1, $2) if $alt =~ /(\w+):(\d+)/;
  $pos2 = $1 if $info =~ /END=(\d+);/;
  $pos2 = $1 if $info =~ /;END=(\d+)/;

  my ( $pos1a, $pos1b ) = ( $pos1, $pos1 );
  $pos1a -= $1 if $info =~ /CIPOS=-*(\d+),/;
  $pos1b += $1 if $info =~ /CIPOS=-*\d+,(\d+)/;
  
  my ( $pos2a, $pos2b ) = ( $pos2, $pos2 );
  $pos2a -= $1 if $info =~ /CIEND=-*(\d+),/;
  $pos2b += $1 if $info =~ /CIEND=-*\d+,(\d+)/;
  
  my $ori = "XX";
  if ($id =~ /BND/) {
      $ori = "TH" if $alt =~ /^\w+\[/;
      $ori = "HT" if $alt =~ /^\]/;
      $ori = "TT" if $alt =~ /^\w+\]/;
      $ori = "HH" if $alt =~ /^\[/;
  } else {
      $ori = "TH" if $id =~ /DEL/;
      $ori = "TH" if $id =~ /INS/;
      $ori = "HT" if $id =~ /DUP/;
      $ori = "TT" if $id =~ /INV/ and $info =~ /INV3/;
      $ori = "HH" if $id =~ /INV/ and $info =~ /INV5/;
  }
  
  my $type = $1 if $info =~ /SVTYPE=(\w+)/;
  $type = "BND" if $type eq "INV";
  $type = "BND" if $type eq "TRA";

  $type = "X";
  
  if ($pos1>$pos2 and $chr eq $chr2) {
    ($pos1a, $pos1b, $pos2a, $pos2b) = ($pos2a, $pos2b, $pos1a, $pos1b);
    $ori = reverse($ori);
  }
  if ($chr2 < $chr) {
    ($chr, $chr2, $pos1a, $pos1b, $pos2a, $pos2b) = ($chr2, $chr, $pos2a, $pos2b, $pos1a, $pos1b);
    $ori=reverse($ori);
  }
  
  my @manta_ids;
  
  $type = "X";

  foreach my $p1a (keys %{$svs{$type}{$chr1}{$chr2}{$ori}}) {
    foreach my $p1b (keys %{$svs{$type}{$chr1}{$chr2}{$ori}{$p1a}}) {
      next unless $p1a <= $pos1b and $p1b >= $pos1a;
      foreach my $p2a (keys %{$svs{$type}{$chr1}{$chr2}{$ori}{$p1a}{$p1b}}) {
	foreach my $p2b (keys %{$svs{$type}{$chr1}{$chr2}{$ori}{$p1a}{$p1b}{$p2a}}) {
	  next unless $p2a <= $pos2b and $p2b >= $pos2a;
	  if ( ($chr1=~ /^$chr2$/) and ($pos2-$pos1<1000) and ($alt !~ 'INS')) {
	    if ( $p1a <= ( $pos2b - $dist ) and $p2b >= ( $pos1a + $dist ) ) {
	      $line[6] = "RefOverlap" if $line[6] eq "PASS";
	      $line[6] .= ",RefOverlap" if $line[6] ne "PASS" and $line[6] !~ /RefOverlap/;
	      push @manta_ids, $svs{$type}{$chr1}{$chr2}{$ori}{$p1a}{$p1b}{$p2a}{$p2b};
	    }
	  } else {
	    $line[6] = "RefOverlap" if $line[6] eq "PASS";
	    $line[6] .= ",RefOverlap" if $line[6] ne "PASS" and $line[6] !~ /RefOverlap/;
	    push @manta_ids, $svs{$type}{$chr1}{$chr2}{$ori}{$p1a}{$p1b}{$p2a}{$p2b};
	  }
	}
      }
    }
  }
  $line[7] .= ";OVERLAPIDS=" . join(",", @manta_ids) if scalar( @manta_ids ) > 0;
  print join("\t", @line) . "\n";
}
close MANTA2;