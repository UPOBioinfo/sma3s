#!/usr/bin/perl
use strict;

my %F;
my (@tsv) = `ls set_*_UniInv_source_go_goslim_source_summary.tsv`;
foreach my $tsv (@tsv) {
  my $a1;

  open in, $tsv;
  while (<in>) {
    chomp;

    next if (/^$/ || /^#UniProt Keyword categories/ || /^#GO Slim/);
    if (/^#/) {
      $a1 = $_;
    } else {
      my ($annot, $freq, $a3) = split/\t/;
      if ($annot =~ /^GO:/) {
        $F{$a1}{"$annot\t$freq"} += $a3;        
      } else {
        $F{$a1}{$annot} += $freq;
      }
    }
  }
  close in;
}

# Results
foreach my $group (sort keys %F) {
  print "$group\n";
  foreach my $a (sort {$F{$group}{$b} <=> $F{$group}{$a}} keys %{$F{$group}}) {
    print "$a\t$F{$group}{$a}\n";    
  }
  print "\n";
}
