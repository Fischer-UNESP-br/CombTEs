# Developed by Carlos Fischer - 14.06.2017
# This script is used to split a chromosome into subsequences with an overlap region between two subsequences
# Usage: perl splitOverlap.pl <fileNameChromosome>

use Bio::SeqIO;
use Bio::SeqUtils;


my $chromo = $ARGV[0];

my $length  = 10000; # length in NTs of the subsequences
my $overlap =  1000; # length in NTs of the overlap region

my $fileIn  = Bio::SeqIO->new(-file => $chromo, -format => 'Fasta');
my $fileOut = Bio::SeqIO->new(-file => ">$chromo\_subseqs.fasta", -format => 'Fasta');

while (my $seqIN = $fileIn->next_seq() ) {
	my $seqNT = $seqIN->seq();
	my $seqLength = length($seqNT);

	my $numSub = 0;
	my $from = 0;
	while ($from < $seqLength) {
		$numSub++;
		my $idSeq = "Subsequence-".$numSub;
		my $subSeq = substr($seqNT,$from,$length);

		$seqOut = Bio::Seq->new(-seq => $subSeq, -id  => $idSeq);
		$fileOut->write_seq($seqOut);

		my $fromAnt = $from;
		$from = $fromAnt + $length - $overlap;
	}
}

