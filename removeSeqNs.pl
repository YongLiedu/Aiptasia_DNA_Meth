#!/usr/bin/perl -w

=pod
Author: Yong Li
Writing: Mar 15, 2014

This script is for ...
=cut

my $thisScript = $0; $thisScript = $1 if ($0 =~ /\\([^\\]+)$/ || $0 =~ /\/([^\/]+)$/);

if ( @ARGV != 2 || $ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "--help" ) {
	print "\n$0:\n";
    print "Usage: ./$thisScript in.fasta out.fasta.\n";
	print "This script is to remove sequences with only Ns from fasta file.\n\n";
    exit(1);
}

print "\nRunning $thisScript for $ARGV[0] now ...\n";

open( IN,  "<$ARGV[0]" ) or die "$!\n";
open( OUT, ">$ARGV[1]" ) or die "$!\n";

my $o = 0;
my $cache = "";

my $in = <IN>;
print OUT "$in";

while ( defined( my $in = <IN> ) ) {
    chomp($in);
    if ($in =~ /^>/){
		$o = 0;
		$cache = $in . "\n";
	} elsif ($o){
		print OUT "$in\n";
	} else {
		if ($in =~ /[ACGT]/i){
			$o = 1;
			print OUT "$cache" if $cache ne "";
			print OUT "$in\n";
		}else{
			$cache .= $in . "\n";
		}
	}
}

close(IN);
close(OUT);

print "    ... $thisScript running finished.\n\n";

