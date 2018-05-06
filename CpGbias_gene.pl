#!/usr/bin/perl -w

=pod
Author: Yong Li
Writing: April 13, 2015

This script is to calculate CpG bias for every gene including introns.

The CpG bias of a sequence is defined as the ratio of the observed frequency of CpG dinucleotides divided by the expected frequency of CpG dinucleotides, where the expected number of CpG dinucleotides is the product of the frequency of C and G nucleotides in a given sequence. When no Cs or no Gs are observed, the CpG bias is arbi- trarily set to one.

Foret, S., Kucharski, R., Pittelkow, Y., Lockett, G.A., and Maleszka, R. (2009). Epigenetic regulation of the honey bee transcriptome: unravelling the nature of methylated genes. BMC Genomics 10, 472.

=cut

use Bio::SeqIO;

my $thisScript = $0; $thisScript = $1 if ($0 =~ /\\([^\\]+)$/ || $0 =~ /\/([^\/]+)$/);

if ( @ARGV != 2 || $ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "--help" ) {
	print "\n$0:\n";
    print "Usage: ./$thisScript genes.fa out_bias.txt.\n";
	print "This script is to calculate CpG bias for every gene.\n\n";
    exit(1);
}

print "\nRunning $thisScript for $ARGV[0] now ...\n";

my $seqio = Bio::SeqIO->new(-file => "<$ARGV[0]", -format => 'fasta', verbose => -1);
open( OUT, ">$ARGV[1]" ) or die "$!\n";

while ( my $seqi = $seqio->next_seq ) {
    my $seq = $seqi->seq;
	my $cg = $c = $g = 0;
	my $slen = $seqi->length;
	while($seq =~ /CG/ig){$cg++;}
	while($seq =~ /C/ig){$c++;}
	while($seq =~ /G/ig){$g++;}
	if ($c == 0 || $g == 0){
		my $cpgbias = 1;
	}
	else{
		my $cpgbias = ($cg * $slen)/($c * $g);
	}
	print OUT "$cpgbias\n";
}
close(OUT);

print "    ... $thisScript running finished.\n\n";

