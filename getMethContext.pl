#!/usr/bin/perl -w

=pod
Author: Yong Li
Writing: Apr 16, 2015

This script is to return the base combinations (CpG, CHG, CHH) of the methylated sites.

=cut

use Bio::SeqIO;
use Parallel::ForkManager;

my $thisScript = $0; $thisScript = $1 if ($0 =~ /\\([^\\]+)$/ || $0 =~ /\/([^\/]+)$/);

if ( @ARGV != 3 || $ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "--help" ) {
	print "\n$0:\n";
    print "Usage: ./$thisScript meth.sites.txt genome.fa out_file.\n";
	print "This script is to return the base combinations (CpG, CHG, CHH) of the methylated sites.\n\n";
    exit(1);
}

print "\nRunning $thisScript for $ARGV[0] now ...\n";

my $MAX_PROCESSES = 30;
my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

open( POS, "<$ARGV[1]" ) or die "$!\n";
my @pos = <POS>;
chomp @pos;
close(POS);

unlink($ARGV[2]) if (-e $ARGV[2]);

Fork:
foreach my $in (@pos){
	my $pid = $pm->start and next Fork;
	my @arr = split /\t/, $in;
	$arr[0] =~ s/size/|size/;
	my $seqio = Bio::SeqIO->new(-file => "<$ARGV[0]", -format => 'fasta', -verbose => -1);
	while ( my $seqi = $seqio->next_seq ) {
		if ($arr[0] eq $seqi->id){
			my $seq = $seqi->subseq($arr[1]-2,$arr[1]+2);
			$seq = uc($seq);
			#  plus strand:  ATmCG TC
			# minus strand:  TA GCmAG
			if (substr($seq,2,1) eq 'C'){
				$seq = substr($seq,2);
			}
			elsif(substr($seq,2,1) eq 'G'){
				$seq = substr($seq,0,3);
				$seq = tr/AGCT/TCGA/;
				$seq = reverse($seq);
			}
			else{
				print "warning: the methylated base is not C!\n";
				next;
			}
			my $mcontext = "";
			if ($seq =~ /^CG/){
				$mcontext = 'CpG';
			}
			elsif($seq =~ /^C[^G]G/){
				$mcontext = 'CHG';
			}
			else{
				$mcontext = 'CHH';
			}
			open( OUT, ">>$ARGV[2]" ) or die "$!\n";
			print OUT "$arr[0]\t$arr[1]\t$mcontext\n";
			close(OUT);
			last;
		}
	}
	$pm->finish;
}
$pm->wait_all_children;

print "    ... $thisScript running finished.\n\n";

