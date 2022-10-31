# Developed by Carlos Fischer - 20 Oct 2020
# Updated: 02 Aug 2022

# To extract and format the predictions from the output file of RPS-Blast, excluding redundancies when 2 or more predictions of different domain types are in the same region - the best one (with the lowest e-value) is mantained in the analyses.

# Redundancy may occur when there is an overlapping region (total or partial) between any 2 PREDICTIONS; the redundancy is verified by (i) calculating (based on the "$maxPercOutOverlap" percentage) the maximum size ("$maxSizeOutOverlap") allowed for the region outside the overlapping region to consider that there is an overlap between the 2 predictions and (ii) verifying for both predictions if the regions outside the overlapping region is smaller than a predefined value ("$maxRegionOutOverlap").

# This script can be used to analyse both a whole long sequence and a set of small sequences.

# Usage: perl extractRPSB.pl output_file_rpsblast

# input:  file with results from RPS-Blast ("output_file_rpsblast")
# output: two files - one with ALL predictions and other with predictions without redundancies, for ALL considered TE types.

# It considers the use of the "RPSTBLASTN 2.10.1+" version (or other versions that maintain the same output format).

# THIS SCRIPT WOULD BE RUN JUST ONCE
# UNLESS the user changes the values of some parameters used here.

# This script uses 4 parameters:
# - 1 of them can be changed in "ParamsGeneral.pm":
#	- @TEtypes: the used TE types.

# - 3 others can be changed in "ParamsRpsblast.pm":
#	- $maxPercOutOverlap: percentage used to calculate the maximum size ("$maxSizeOutOverlap") allowed for the region outside the 		overlapping region to consider that there is an overlap between 2 predictions (default = 0.20 = 20%);
#	- $maxRegionOutOverlap: (predefined) maximum length allowed for the region(s) outside the overlapping region between 2 		predictions (default = 100 - a lower value would increase the number of predictions to be considered in the following  		analyses);
#	- %domains: conserved domains considered in the analyses (after running RPS-Blast) - defined in "ParamsRpsblast.pm" package.

###########################################################################################

use strict;
use warnings;

use Cwd qw(getcwd);
use lib getcwd(); 

use ParamsGeneral  qw(@TEtypes);
use ParamsRpsblast qw($maxPercOutOverlap $maxRegionOutOverlap %domains);

###########################################################################################

sub excludeRedund { # to exclude redundancies - for predictions inside a same sequence

    (my @predictions) = @_;

    my $qttR = scalar(@predictions);
    my @auxPreds = ();
    for (my $i = 0; $i < $qttR; $i++) {
#				 FROM--fff---TO--ttt---LENGTH
	if ($predictions[$i] =~ /FROM--(\d+)---TO/) {
	    push (@auxPreds, {line => $predictions[$i], from => $1});
	}
    }
    my @sortedPredicts = sort{$a->{from} <=> $b->{from}} @auxPreds;

    for (my $i = 0; $i < $qttR-1; $i++) { # comparing predictions all against all
	for (my $j = $i+1; $j < $qttR; $j++) {
	    if ($sortedPredicts[$i]->{line} =~ /-->EXCLUD/) { last; }
	    if ($sortedPredicts[$j]->{line} !~ /-->EXCLUD/) {
		my $resp = verifyOverlap ($sortedPredicts[$i]->{line},$sortedPredicts[$j]->{line});
		if    ($resp eq "II") { $sortedPredicts[$i]->{line} .= "-->EXCLUD"; }
		elsif ($resp eq "JJ") { $sortedPredicts[$j]->{line} .= "-->EXCLUD"; }
		else		      { last; } # stop comparison when no overlap is identified
	    }
	} # FOR ($j = $i+1; ...
    } # FOR ($i = 0; ...

    my @sortPredsNoRedund = ();
    for (my $i = 0; $i < $qttR; $i++) {
	if ($sortedPredicts[$i]->{line} !~ /-->EXCLUD/) { push (@sortPredsNoRedund, $sortedPredicts[$i]->{line}); }
    }

    return @sortPredsNoRedund;
} # END of SUB "excludeRedund"


sub verifyOverlap { # to verify overlap between 2 predictions and, if so, retrieve the ID of the one to be disregarded; this prediction would be the one with the highest e-value or, for the same e-values, the shortest prediction

    (my $seq1, my $seq2) = @_;

    my ($from1, $to1, $length1, $eval1, $sense1);
    my ($from2, $to2, $length2, $eval2, $sense2);

#	     $seq= DOMAIN--dom---FROM--ff---TO--tt---LENGTH--ll---EVALUE--ev---SCORE--sc---SENSE--dr
    if ($seq1 =~ /^DOMAIN--.*---FROM--(\d+)---TO--(\d+)---LENGTH--(\d+)---EVALUE--(.*)---SCORE--.*---SENSE--(.*)/) {
	    $from1   = $1;
	    $to1     = $2;
	    $length1 = $3;
	    $eval1   = $4;
	    $sense1  = $5;
    }
    if ($seq2 =~ /^DOMAIN--.*---FROM--(\d+)---TO--(\d+)---LENGTH--(\d+)---EVALUE--(.*)---SCORE--.*---SENSE--(.*)/) {
	    $from2   = $1;
	    $to2     = $2;
	    $length2 = $3;
	    $eval2   = $4;
	    $sense2  = $5;
    }

    my $respOverl = "NONE"; # initially, it considers that there is NOT overlap between seq_1 and seq_2; otherwise, it retrieves the ID of the prediction to be disregarded
    if ($sense1 eq $sense2) {
	my $hasOverlap = "no";
#     testing if the predictions are apart from each other
	if ( ($to1 > $from2) and ($to2 > $from1) ) {
#	  calculating the maximum size ("$maxSizeOutOverlap" - according to the length of the smallest prediction) outside the 		  overlapping region to consider that there is an overlap between 2 predictions:
	    my $smallestLeng = $length1;
	    if ($length2 < $length1) { $smallestLeng = $length2; }
	    my $maxSizeOutOverlap = $maxPercOutOverlap * $smallestLeng;

#	  IF: (seq_2 inside seq_1) OR (seq_1 inside seq_2) --->> overlap between predictions
	    if ( ( ($from1 <= $from2) and ($to2 <= $to1) ) or ( ($from2 <= $from1) and ($to1 <= $to2) ) ) { $hasOverlap = "yes"; }
#	  testing the maximum allowed region outside the overlap
	    elsif ( ( (abs($from1-$from2) <= $maxSizeOutOverlap) and (abs($from1-$from2) <= $maxRegionOutOverlap) ) or
		    ( (abs($to1-$to2)     <= $maxSizeOutOverlap) and (abs($to1-$to2)     <= $maxRegionOutOverlap) )
		  ) # --->> if so, overlap between predictions
		  { $hasOverlap = "yes"; }
	} # IF ( ($to1 > $from2) and ($to2 > $from1) )

	if ($hasOverlap eq "yes") { # there is overlap between seq_1 and seq_2
	    if    ($eval1 < $eval2) { $respOverl = "JJ"; }
	    elsif ($eval1 > $eval2) { $respOverl = "II"; }
	    else {
		if ($length1 >= $length2) { $respOverl = "JJ"; }	
		else			  { $respOverl = "II"; }
	    }
	} # IF ($hasOverlap eq "yes")
    } # IF ($sense1 eq $sense2)

    return "$respOverl";
} # END of SUB "verifyOverlap"

###########################################################################################


unless(@ARGV) { die "\nUSAGE: perl $0 file_from_RpsBlast\n\n"; }

my $fileRpsb = $ARGV[0];
open (INPRED, $fileRpsb) or die "Can't open $fileRpsb!\n";

my $allPred = "ConservedDomains_ALL.pred";
open (ALLPRED,">$allPred") or die "Can't open $allPred";
print ALLPRED "ALL predictions of RPSBlastn for conserved domains, from file \"$fileRpsb\".\n\n";

my $predNoRedund = "ConservedDomains_NoRedund.pred";
open (PREDNO, ">$predNoRedund") or die "Can't open $predNoRedund";
print PREDNO "Predictions of RPSBlastn for conserved domains of interest, from file \"$fileRpsb\".\n\n";

my $currDomsInterest = "";
foreach my $type (@TEtypes) {
	if ($currDomsInterest ne "") { $currDomsInterest .= "|"; }
	$currDomsInterest .= $domains{$type};
}

my ($countDom, $idSeq);
my @sortedDomsInter = (); # for each ID and its domains of interest

my $line = readline(INPRED); # 1st line
while (not eof INPRED) {
	$countDom = 0;
	$line = readline(INPRED); chomp $line; # 2nd line
	if ($line =~ /Query: (.*)/)  { $idSeq = $1; }

	$line = readline(INPRED); # 3th line: "# Database: "
	$line = readline(INPRED); chomp $line; # 4th line must be "0 hits found" OR "Fields: ..."
	if ($line !~ /# Fields:/) {
		if ($line !~ /# 0 hits found/) { die ("PROBLEMS in the RPSBLAST's output file: $line!!! - maybe you did not use the \"RPSTBLASTN 2.10.1+\" (or a compatible) version!!") }
	}
	else { $line = readline(INPRED); }

	my @allDoms = ();	# for ALL found domains
	my @domsInterest = ();	# for the current domains of interest
	$line = readline(INPRED);
	while ( ($line !~ /# RPSTBLASTN/) and (not eof INPRED) ) {
		$countDom++;

		my @array1  = split '\t',$line;
		my $nameDom = $array1[1];
		my $frame   = $array1[2];
		my $evalue  = $array1[3];
		my $score   = $array1[4];
		my $fromDom = $array1[5];
		my $toDom   = $array1[6];

		my $sense   = "Direct";
		if ($frame < 0) { # identified in reverse sense
			$sense   = "Reverse";
			my $aux     = $fromDom;
			$fromDom = $toDom;
			$toDom   = $aux;
		}
		my $length = $toDom - $fromDom + 1;

		my @array2 = split ', ',$nameDom;
		my $dom = $array2[1];
		my $domDescr = "DOMAIN--$dom---FROM--$fromDom---TO--$toDom---LENGTH--$length---EVALUE--$evalue---SCORE--$score---SENSE--$sense";
		push (@allDoms, {line => $domDescr, from => $fromDom});

#	      selecting the domains of interest;
#		it also tests the possible prediction of the "integrase" ("rve") domain (when considered in the analyses)
#		if this domain is NOT being considered in the analyses, this fact has no influence here
		if ( ($dom eq "rve") or ($dom =~ /$currDomsInterest/) ) { push (@domsInterest, $domDescr); }

		$line = readline(INPRED);
	} # WHILE ( ($line !~ /# RPSTBLASTN ...

	my $seqID = ">>>SEQUENCE: $idSeq";
	print ALLPRED "$seqID\n";

	if ($countDom == 0) {
		print ALLPRED "NO DOMAIN\n";
		print ALLPRED "###\n";
	}
	else {
		my @sortedAllDoms = sort{$a->{from} <=> $b->{from}} @allDoms;
		foreach my $lineAux (@sortedAllDoms) { print ALLPRED $lineAux->{line}."\n"; }
		print ALLPRED "###\n";
		my $qttDomInter = scalar(@domsInterest);
		if ($qttDomInter != 0) {
#		    excluding redundancies for predictions inside a same sequence
			my @sortedInter = excludeRedund (@domsInterest);
			push (@sortedDomsInter, $seqID);
			foreach my $lineAux (@sortedInter) { push (@sortedDomsInter, $lineAux); }
			push (@sortedDomsInter, "###");
		}
	}
} # WHILE (not eof INPRED)
close (INPRED);

foreach $line (@sortedDomsInter) { print PREDNO "$line\n"; }

close (ALLPRED);
close (PREDNO);

