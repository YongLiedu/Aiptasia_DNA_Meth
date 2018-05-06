#!/usr/bin/perl -w

=pod
Author: Yong Li
Writing: April 13, 2015

This script is to calculate CpG bias for every gene including introns&UTRs in gff3 file.

The CpG bias of a sequence is defined as the ratio of the observed frequency of CpG dinucleotides divided by the expected frequency of CpG dinucleotides, where the expected number of CpG dinucleotides is the product of the frequency of C and G nucleotides in a given sequence. When no Cs or no Gs are observed, the CpG bias is arbi- trarily set to one.

Foret, S., Kucharski, R., Pittelkow, Y., Lockett, G.A., and Maleszka, R. (2009). Epigenetic regulation of the honey bee transcriptome: unravelling the nature of methylated genes. BMC Genomics 10, 472.

=cut

use Bio::SeqIO;
use Parallel::ForkManager;

my $thisScript = $0; $thisScript = $1 if ($0 =~ /\\([^\\]+)$/ || $0 =~ /\/([^\/]+)$/);

if ( @ARGV != 3 || $ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "--help" ) {
	print "\n$0:\n";
    print "Usage: ./$thisScript genome.fa gff_file.gff3 out_bias.txt.\n";
	print "This script is to calculate CpG bias for every gene.\n\n";
    exit(1);
}

print "\nRunning $thisScript for $ARGV[0] now ...\n";

my $MAX_PROCESSES = 30;
my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

system("grep '\tgene\t' $ARGV[1] > tmpXxTdErf.tmp");
open( GFF, "<tmpXxTdErf.tmp" ) or die "$!\n";
my @gff = <GFF>;
close(GFF);
unlink("tmpXxTdErf.tmp");

unlink($ARGV[2]);

Fork:
foreach my $in (@gff){
	my $pid = $pm->start and next Fork;
	next if $in =~ /^#/;
	next if $in =~ /^$/;
	my @arr = split /\t/, $in;
	next if ($arr[2] ne 'gene');
	$arr[0] =~ s/size/|size/;
	my $seqio = Bio::SeqIO->new(-file => "<$ARGV[0]", -format => 'fasta', -verbose => -1);
	while ( my $seqi = $seqio->next_seq ) {
		#print $arr[0] . "\t" . $seqi->id . "\n";
		if ($arr[0] eq $seqi->id){
			my $seq = $seqi->subseq($arr[3],$arr[4]);
			if ($arr[6] eq '-'){
				$seq =~ tr/AGCT/TCGA/;
				$seq = reverse($seq);
			}
			my $cg = $c = $g = 0;
			my $slen = length($seq);
			while($seq =~ /CG/ig){$cg++;}
			while($seq =~ /C/ig){$c++;}
			while($seq =~ /G/ig){$g++;}
			my $cpgbias;
			if ($c == 0 || $g == 0){
				$cpgbias = 1;
			}
			else{
				$cpgbias = ($cg * $slen)/($c * $g);
			}
			my @arr2 = split /\;/, $arr[8];
			$arr2[0] =~ s/ID=//;
			open( OUT, ">>$ARGV[2]" ) or die "$!\n";
			print OUT $arr2[0] . "\t$cpgbias\n";
			close(OUT);
			last;
		}
	}
	$pm->finish;
}
$pm->wait_all_children;

print "    ... $thisScript running finished.\n\n";

