# Developed by Carlos Fischer - 20 Oct 2020
# Updated: 02 Aug 2022

# To extract and format the predictions from the output files of HMMER and RepeatMasker, excluding redundancies when 2 or more predictions are in the same region - the best one (with the lowest e-value OR the higher SWscore) is mantained in the analyses.
# Interesting for HMMER: the script excludes redundancies when using 2 or more pHMMs.

# Redundancy may occur when there is an overlapping region (total or partial) between any 2 PREDICTIONS; the redundancy is verified by (i) calculating (based on the "$maxPercOutOverlap" percentage) the maximum size ("$maxSizeOutOverlap") allowed for the region outside the overlapping region to consider that there is an overlap between the 2 predictions and (ii) verifying for both predictions if the regions outside the overlapping region is smaller than a predefined value ("$maxRegionOutOverlap").

# This script can be used to analyse either a whole long sequence or a set of small sequences.

# Usage: perl extractHmmerRM.pl TOOL TE-type file-tetype-tool

# input ("file-tetype-tool"): file with results from HMMER or RepeatMasker for each TE type
# output: two files - one with ALL predictions (no matter the prediction's length) and other with predictions without redundancies (for the last file, the prediction's length is limited to "$minLenPred"), for EACH TE type


# THIS SCRIPT WOULD BE RUN JUST ONCE FOR EACH TE TYPE
# UNLESS the user changes the values of some parameters used here.

# This script uses 4 parameters:
# - 1 of them can be changed in "ParamsGeneral.pm":
#	- @TEtypes: the used TE types.

# - 3 others can be changed in "ParamsHmmerRM.pm":
#	- $maxPercOutOverlap: percentage used to calculate the maximum size ("$maxSizeOutOverlap") allowed for the region outside the 		overlapping region to consider that there is an overlap between 2 predictions (default = 0.20 = 20%);
#	- $maxRegionOutOverlap: (predefined) maximum length allowed for the region(s) outside the overlapping region between 2 		predictions (default = 100 - a lower value would increase the number of predictions to be considered in the following 		analyses);
#	- "$minLenPred": minimum length for a prediction to be included in the analyses (default = 20).

###########################################################################################

use strict;
use warnings;

use Cwd qw(getcwd);
use lib getcwd(); 

use ParamsGeneral qw(@TEtypes);
use ParamsHmmerRM qw($maxPercOutOverlap $maxRegionOutOverlap $minLenPred);

###########################################################################################


sub excludeRedund { # to exclude redundancies - for predictions inside a same sequence

    (my $method, my @predictions) = @_;

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
		my $resp = verifyOverlap ($method, 0, $sortedPredicts[$i]->{line},$sortedPredicts[$j]->{line});
		if    ($resp eq "II") { $sortedPredicts[$i]->{line} .= "-->EXCLUD"; }
		elsif ($resp eq "JJ") { $sortedPredicts[$j]->{line} .= "-->EXCLUD"; }
	    }
	} # FOR ($j = $i+1; ...
    } # FOR ($i = 0; ...

    my @sortPredsNoRedund = ();
    for (my $i = 0; $i < $qttR; $i++) {
	if ($sortedPredicts[$i]->{line} !~ /-->EXCLUD/) { push (@sortPredsNoRedund, $sortedPredicts[$i]->{line}); }
    }

    return @sortPredsNoRedund;
} # END of SUB "excludeRedund"

###############################################

sub verifyOverlap { # to verify overlap between 2 predictions and, if so, retrieve the ID of the one to be disregarded; this prediction would be:
# - for HMMER:        the one with the highest e-value or, for the same e-values, the shortest prediction
# - for RepeatMasker: the one with the lowest RMscore or,  for the same RMscores, the shortest prediction

    (my $method, my $offset, my $seq1, my $seq2) = @_;

    my ($from1, $to1, $length1, $evalSWs1, $sense1);
    my ($from2, $to2, $length2, $evalSWs2, $sense2);

# For HMMER:
# $seq= PREDIC---FROM--ff---TO--tt---LENGTH--ll---EVALUE--ev---SCORE--sc---SENSE--dr---TETYPE--type
# For RepeatMasker:
# $seq= PREDIC---FROM--ff---TO--tt---LENGTH--ll---RMSCORE--sc---SENSE--dr---MATCHINGREPEAT--match---TETYPE--type
    if (
	($seq1 =~ /^PREDIC---FROM--(.*)---TO--(.*)---LENGTH--(.*)---EVALUE--(.*)---SCORE--.*---SENSE--(.*)---TETYPE/)    or
	($seq1 =~ /^PREDIC---FROM--(.*)---TO--(.*)---LENGTH--(.*)---RMSCORE--(.*)---SENSE--(.*)---MATCHINGREPEAT/)
       ) {
	$from1    = $1;
	$to1      = $2;
	$length1  = $3;
	$evalSWs1 = $4;
	$sense1   = $5;
    }
    if (
	($seq2 =~ /^PREDIC---FROM--(.*)---TO--(.*)---LENGTH--(.*)---EVALUE--(.*)---SCORE--.*---SENSE--(.*)---TETYPE/)    or
	($seq2 =~ /^PREDIC---FROM--(.*)---TO--(.*)---LENGTH--(.*)---RMSCORE--(.*)---SENSE--(.*)---MATCHINGREPEAT/)
       ) {
	$from2    = $1 + $offset;
	$to2      = $2 + $offset;
	$length2  = $3;
	$evalSWs2 = $4;
	$sense2   = $5;
    }

    my $respOverl = "NONE"; # initially, it considers that there is NOT overlap between seq_1 and seq_2; otherwise, it retrieves the ID of the prediction to be disregarded
    if ($sense1 eq $sense2) {
	my $hasOverlap = "no";
#     testing if the predictions are apart from each other
	if ( ($to1 > $from2) and ($to2 > $from1) ) {
#	  calculating the maximum size ("$maxSizeOutOverlap" - according to the length of the smallest prediction) outside the 		  overlapping region to consider that there is an overlap between 2 predictions:
	    my $referLength = $length1;
	    if ($length2 < $length1) { $referLength = $length2; }
	    my $maxSizeOutOverlap = $maxPercOutOverlap * $referLength;

#	  IF: (seq_2 inside seq_1) OR (seq_1 inside seq_2) --->> overlap between predictions
	    if ( ( ($from1 <= $from2) and ($to2 <= $to1) ) or ( ($from2 <= $from1) and ($to1 <= $to2) ) ) { $hasOverlap = "yes"; }
#	  testing the maximum allowed region outside the overlap
	    elsif ( ( (abs($from1-$from2) <= $maxSizeOutOverlap) and (abs($from1-$from2) <= $maxRegionOutOverlap) ) or
		    ( (abs($to1-$to2)     <= $maxSizeOutOverlap) and (abs($to1-$to2)     <= $maxRegionOutOverlap) )
		  ) # --->> if so, overlap between predictions
		  { $hasOverlap = "yes"; }
	} # IF ( ($to1 > $from2) and ($to2 > $from1) )

	if ($hasOverlap eq "yes") { # there is overlap between seq_1 and seq_2
	    if ($method eq "HMMER") {
		if    ($evalSWs1 < $evalSWs2) { $respOverl = "JJ"; }
		elsif ($evalSWs1 > $evalSWs2) { $respOverl = "II"; }
		else {
		    if ($length1 >= $length2) { $respOverl = "JJ"; }	
		    else		      { $respOverl = "II"; }
		}
	    }
	    elsif ($method eq "RepeatMasker") {
		if    ($evalSWs1 > $evalSWs2) { $respOverl = "JJ"; }
		elsif ($evalSWs1 < $evalSWs2) { $respOverl = "II"; }
		else {
		    if ($length1 >= $length2) { $respOverl = "JJ"; }	
		    else		      { $respOverl = "II"; }
		}
	    }
	} # IF ($hasOverlap eq "yes")
    } # IF ($sense1 eq $sense2)

    return ($respOverl);
} # END of SUB "verifyOverlap"

###########################################################################################


    unless(@ARGV) { die "\nUSAGE: perl $0 TOOL TETYPE file_for_tetype\n\n"; }

    my $tool    = $ARGV[0];
    my $typeIN  = $ARGV[1];
    my $fileIN  = $ARGV[2];

    if ( ($tool ne "HMMER") and ($tool ne "RepeatMasker") ) { die "\nWrong TOOL!!! MUST use \"HMMER\" or \"RepeatMasker\".\n\n"; }

    my $OKtype = 0;
    foreach my $type (@TEtypes) { if ($type eq $typeIN) { $OKtype = 1; } }
    if ($OKtype == 0) { die "\nWrong TE type!!!\n\n"; }
    open (INPRED, $fileIN) or die "\nCan't open \"$fileIN\"!!!\n\n";

    my $allPred = "$typeIN\_$tool\_ALL.pred";
    open (ALLPRED,">$allPred") or die "Can't open $allPred!\n";
    print ALLPRED "ALL predictions of $tool for \"$typeIN\", from file \"$fileIN\".\n\n";

    my $predNoRedund = "$typeIN\_$tool.pred";
    open (PREDNO, ">$predNoRedund") or die "Can't open $predNoRedund";
    print PREDNO "Predictions of $tool for \"$typeIN\", from file \"$fileIN\".\n\n";

    my @allPredicts = ();

    readline(INPRED); readline(INPRED);
    if ($tool eq "RepeatMasker") { readline(INPRED); }

    my ($line, $idSeq);
    while (not eof INPRED) {
	$line = readline(INPRED); chomp $line;
	if ($line !~ /^#/) {
	    my @auxSplit = split ' ',$line;
	    my ($from, $to, $strand, $evalue, $score, $sense, $length, $lineOUT, $matchRM);
	    if ($tool eq "HMMER") {
		$idSeq  = $auxSplit[0];
		$from   = $auxSplit[6];
		$to     = $auxSplit[7];
		$strand = $auxSplit[11];
		$evalue = $auxSplit[12];
		$score  = $auxSplit[13];

		$sense  = "Direct";
		if ($strand eq "-") { # identified in reverse sense
			$sense = "Reverse";
			my $aux= $from;
			$from  = $to;
			$to    = $aux;
		}
		$length = $to - $from + 1;

		$lineOUT = "PREDIC---FROM--$from---TO--$to---LENGTH--$length---EVALUE--$evalue---SCORE--$score---SENSE--$sense---TETYPE--$typeIN";
	    }
	    elsif ($tool eq "RepeatMasker") {
		$score   = $auxSplit[0];
		$idSeq   = $auxSplit[4];
		$from    = $auxSplit[5];
		$to      = $auxSplit[6];
		$strand  = $auxSplit[8];
		$matchRM = $auxSplit[9];

		$sense   = "Direct";
		if ($strand eq "C") { $sense = "Reverse"; } # identified in reverse sense
		$length = $to - $from + 1;

		$lineOUT = "PREDIC---FROM--$from---TO--$to---LENGTH--$length---RMSCORE--$score---SENSE--$sense---MATCHINGREPEAT--$matchRM---TETYPE--$typeIN";
	    }
	    push (@allPredicts, {line => $lineOUT, from => $from, idSeq => $idSeq, leng => $length});
	} # IF ($line !~ /^#/)

    } # WHILE (not eof INPRED)
    close (INPRED);

    my $qttPred = scalar(@allPredicts);
    if ($qttPred == 0) { die "\nInput file has no predictions!!!\n\n"; }

    my @sortedAllPredicts = sort{$a->{idSeq} cmp $b->{idSeq}} @allPredicts;
    $qttPred = scalar(@sortedAllPredicts);

    my $indexAllPredicts = 0;
    while ($indexAllPredicts < $qttPred) {
	my $currID = $sortedAllPredicts[$indexAllPredicts]->{idSeq};
	print ALLPRED ">>>SEQUENCE: $currID\n";
	my @eachSeq = ();
	push (@eachSeq, $sortedAllPredicts[$indexAllPredicts]);
	$indexAllPredicts++;
	while ( ($indexAllPredicts < $qttPred) and ($currID eq $sortedAllPredicts[$indexAllPredicts]->{idSeq}) ) {
	    push (@eachSeq, $sortedAllPredicts[$indexAllPredicts]);
	    $indexAllPredicts++;
	}
	my @sortedEachSeq = sort{$a->{from} <=> $b->{from}} @eachSeq;
	my @filteredSortedEachSeq = ();
	foreach my $lineAux (@sortedEachSeq) {
	    print ALLPRED $lineAux->{line}."\n";
	    if ($lineAux->{leng} >= $minLenPred) { push (@filteredSortedEachSeq, $lineAux->{line}); }
	}
	print ALLPRED "###\n";

	my @excludedFilteredSortedEachSeq = excludeRedund ($tool, @filteredSortedEachSeq); # excluding redundancies
	my $qttExclEachSeq = scalar(@excludedFilteredSortedEachSeq);
	if ($qttExclEachSeq != 0) {
	    print PREDNO ">>>SEQUENCE: $currID\n";
	    foreach my $lineAux (@excludedFilteredSortedEachSeq) { print PREDNO "$lineAux\n"; } 
	    print PREDNO "###\n";
	}
    } # WHILE ($indexAllPredicts < $qttPred)

    close (ALLPRED);
    close (PREDNO);

