# Developed by Carlos Fischer - 20 Oct 2020
# Updated: 29 Jan 2021

# To extract and format the predictions from the output files of RPS-Blast (and excluding redundancies).

# Redundancy may occur when there is an overlapping region (total or partial) between 2 PREDICTIONS; the redundancy is verified by (i) calculating (based on the "$maxPercOutOverlap" percentage) the maximum size allowed for the region outside the overlapping region to consider the overlap as (possibly) acceptable and (ii) verifying for both predictions if the regions outside the overlapping region is smaller than a predefined value ("$maxRegionOutOverlap").

# This script can be used when analysing a whole long sequence or subsequences of a long sequence.

# Usage: perl extractRPSB.pl output_file_rpsblast

# input:  file with results from RPS-Blast
# output: two files - one with ALL predictions and other with predictions without redundancies, for ALL considered TE superfamilies

# It considers the use of "RPSTBLASTN 2.10.1+" (and other versions that maintain the same output format)

# THIS SCRIPT WOULD BE RUN JUST ONCE
# UNLESS the user changes the values of some parameters used here.

# This script uses 8 parameters:
# - 4 of them can be changed in "ParamsGeneral.pm":
#	- @superfamilies: the used superfamilies;
#	- $usingSubseqs: when using subsequences of a long sequence ($usingSubseqs = "yes") or directly on a whole long sequence 		($usingSubseqs = "no");
#	- $overlapBetwSubseqs: length of the original overlapping region between subsequent subsequences;
#	- $limitForSubseq: equals to "$lengthSubseq" (length of the subsequences) minus "$overlapBetwSubseqs" (length of the 		overlapping region between subsequences).
# - 4 others can be changed in "ParamsRpsblast.pm":
#	- $varExtr: to allow small variation between respective extremities for cases of overlap between 2 predictions (default=30);
#	- $maxPercOutOverlap: percentage used to calculate the maximum size allowed for the region outside the overlapping region, to 		consider the overlap as (possibly) acceptable (default = 0.20 = 20%);
#	- $maxRegionOutOverlap: (predefined) maximum length allowed for the region(s) outside the overlapping region between 2 		predictions (default = 100 - a lower value would increase the number of predictions to be considered in the following 		analyses);
#	- %domains: conserved domains considered in the analyses (after running RPS-Blast), for:
#		- Bel:   DUF1759, Peptidase_A17, DUF1758, RT_pepA17;
#		- Copia: Retrotran_gag, gag_pre-integrs, RVT_2, RNase_HI_RT_Ty1;
#		- Gypsy: Retrotrans_gag, gag-asp_proteas, retropepsin_like, RP_Saci_like, RVP_2, RT_LTR, RVT_3, RNase_HI_RT_Ty3, 			RNase_HI_like;
#		- and:   rve (Integrase), for Bel, Copia, and Gypsy.


# For the case of subsequences of a long sequence, the IDs of each subsequence must be in the format "(XX)Subsequence-NumSubseq": "XX" means anything, or even nothing, and "NumSubseq" is a number in ascending order for each subsequence.
# Also available, the "splitOverlap.pl" script splits a long sequence in small ones with an overlapping region between 2 subsequent subsequences.

###########################################################################################

use strict;
use warnings;

use Cwd qw(getcwd);
use lib getcwd(); 

use ParamsGeneral  qw(@superfamilies $usingSubseqs $overlapBetwSubseqs $limitForSubseq);
use ParamsRpsblast qw($varExtr $maxPercOutOverlap $maxRegionOutOverlap %domains);

###########################################################################################

sub excludeRedund { # to exclude redundancies

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
	    if ($sortedPredicts[$i]->{line} =~ /-->OUT/) { last; }
	    if ($sortedPredicts[$j]->{line} !~ /-->OUT/) {
		my $resp = verifyOverlap (0, $sortedPredicts[$i]->{line},$sortedPredicts[$j]->{line});
		if    ($resp eq "II") { $sortedPredicts[$i]->{line} .= "-->OUT"; }
		elsif ($resp eq "JJ") { $sortedPredicts[$j]->{line} .= "-->OUT"; }
	    }
	} # FOR ($j = $i+1; ...
    } # FOR ($i = 0; ...

    my @sortPredsNoRedund = ();
    for (my $i = 0; $i < $qttR; $i++) {
	if ($sortedPredicts[$i]->{line} !~ /-->OUT/) { push (@sortPredsNoRedund, $sortedPredicts[$i]->{line}); }
    }

    return @sortPredsNoRedund;
} # END of SUB "excludeRedund"


sub verifyOverlap { # to verify overlap between 2 predictions and, if so, retrieve the ID of the one to be disregarded; this prediction would be the one with the highest e-value or, for the same e-values, the shortest prediction

    (my $offset, my $seq1, my $seq2) = @_;

    my ($from1, $to1, $length1, $eval1, $sense1);
    my ($from2, $to2, $length2, $eval2, $sense2);

# $seq= DOMAIN--dom---FROM--ff---TO--tt---LENGTH--ll---EVALUE--ev---SCORE--sc---SENSE--dr
    if ($seq1 =~ /^DOMAIN--.*---FROM--(.*)---TO--(.*)---LENGTH--(.*)---EVALUE--(.*)---SCORE--.*---SENSE--(.*)/) {
	    $from1   = $1;
	    $to1     = $2;
	    $length1 = $3;
	    $eval1   = $4;
	    $sense1  = $5;
    }
    if ($seq2 =~ /^DOMAIN--.*---FROM--(.*)---TO--(.*)---LENGTH--(.*)---EVALUE--(.*)---SCORE--.*---SENSE--(.*)/) {
	    $from2   = $1 + $offset;
	    $to2     = $2 + $offset;
	    $length2 = $3;
	    $eval2   = $4;
	    $sense2  = $5;
    }

    my $respOverl = "NONE"; # define that there is NOT overlap between seq_1 and seq_2; otherwise, it retrieves the ID of the prediction to be disregarded
    if ($sense1 eq $sense2) {
	my $hasOverlap = "no";
#     testing if the predictions are apart from each other
	if ( ($to1 > $from2) and ($to2 > $from1) ) {
#	  calculating the maximum size ("$maxSizeOutOverlap" - according to the length of the smallest prediction) outside the 		  overlapping region to consider the overlap as (possibly) acceptable:
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
print ALLPRED "All predictions of RPSBlastn for conserved domains, from file \"$fileRpsb\".\n\n";

my $predNoRedund = "ConservedDomains_NoRedund.pred";
open (PREDNO, ">$predNoRedund") or die "Can't open $predNoRedund";
print PREDNO "Predictions of RPSBlastn for conserved domains of interest, from file \"$fileRpsb\".\n\n";

my $domsInterestNow = "";
my $firstSp = "yes";
foreach my $spfam (@superfamilies) {
	if ($firstSp eq "no") { $domsInterestNow .= "|"; }
	$domsInterestNow .= $domains{$spfam};
	$firstSp = "no";
}

my ($countDom, $numSubseq);
my @sortedDomsInter = (); # for each ID and its domains of interest (for ALL subseqs together when analysing subseqs of a long sequence)

my $line = readline(INPRED); # 1st line
while (not eof INPRED) {
	$countDom = 0;
	$line = readline(INPRED); chomp $line; # 2nd line
	if    ($line =~ /Subseq-(\d+)/) { $numSubseq = $1; } # for subsequences
	elsif ($line =~ /Query: (.*)/)  { $numSubseq = $1; } # for a whole sequence

	$line = readline(INPRED); # 3th line: "# Database: "
	$line = readline(INPRED); chomp $line; # 4th line must be "0 hits found" OR "Fields: ..."
	if ($line !~ /# Fields:/) {
		if ($line !~ /# 0 hits found/) { die ("PROBLEMS in RPSBLAST file: $line!!!") }
	}
	else { $line = readline(INPRED); }

	my @allDoms = ();	# for ALL found domains: for just a subsequence OR for the whole long sequence
	my @domsInterest = ();	# for the domains of interest: for just a subsequence OR for the whole long sequence
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

#	      selecting domains of interest
		if ( ($dom eq "rve") or ($dom =~ $domsInterestNow) ) { push (@domsInterest, $domDescr); }

		$line = readline(INPRED);
	} # WHILE ( ($line !~ /# RPSTBLASTN ...

	my $seqID;
	if ($countDom != 0) {
		if ($usingSubseqs eq "yes") { $seqID = ">>>SUBSEQUENCE_$numSubseq"; }
		else			    { $seqID = ">>>SEQUENCE: $numSubseq"; }

		print ALLPRED "$seqID\n";
		my @sortedAllDoms = sort{$a->{from} <=> $b->{from}} @allDoms;
		foreach my $lineAux (@sortedAllDoms) { print ALLPRED $lineAux->{line}."\n"; }
		print ALLPRED "###\n";

		my $qttDomInter = scalar(@domsInterest);
		if ($qttDomInter != 0) {
#		    excluding redundancies for predictions inside a subsequence
			my @sortedInter = excludeRedund (@domsInterest);
			if ($usingSubseqs eq "yes") { push (@sortedDomsInter, $seqID); }
			foreach my $lineAux (@sortedInter) { push (@sortedDomsInter, $lineAux); }
		}
	}
} # WHILE (not eof INPRED)
close (INPRED);

if ($usingSubseqs eq "yes") { # ONLY for subsequences from a long sequence
##    excluding (possible) redundancies for predictions into the overlapping region between subsequent subsequences
	my $qtt = scalar(@sortedDomsInter);
	my $fstIndx = 0;
	while ($fstIndx < $qtt) {
		my $predInOverlReg = "no";
		my $secIndx = $fstIndx + 1;
		$line = $sortedDomsInter[$secIndx];
#				>>>SUBSEQUENCE_numSubseq
		while ( ($line !~ /SUBSEQUENCE_(\d+)/) and ($secIndx < $qtt) ) {
			# firstly, check if any domain prediction is into overlapping region
			if ($line =~ /TO--(\d+)---LENGTH/) { if ($1 > $limitForSubseq) { $predInOverlReg = "yes"; } }
			$secIndx++;
			if ($secIndx < $qtt) { $line = $sortedDomsInter[$secIndx]; }
		}

		my $trdIndx = $secIndx + 1;
		if ( ($secIndx < $qtt) and ($predInOverlReg eq "yes") ) {
			$line = $sortedDomsInter[$trdIndx];
			while ( ($line !~ /SUBSEQUENCE_(\d+)/) and ($trdIndx < $qtt) ) {
				$trdIndx++;
				if ($trdIndx < $qtt) { $line = $sortedDomsInter[$trdIndx]; }
			}
		}

		if ($predInOverlReg eq "yes") { # checking possible overlap between predictions of subsequent subsequences
		   for (my $i = ($fstIndx+1); $i < $secIndx; $i++) { # compare predictions from 1 subseq to all preds from next subseq
		      if ($sortedDomsInter[$i] =~ /FROM--(\d+)---TO--(\d+)---LENGTH/) {
			my $fromI = $1;
			my $toI   = $2;
			if ($toI > $limitForSubseq) {
			   my $toJ;
			   for (my $j = ($secIndx+1); $j < $trdIndx; $j++) {
				if ($sortedDomsInter[$i] =~ /-->OUT/) { last; }
				if ($sortedDomsInter[$j] !~ /-->OUT/) { # for COMPLETE or PARTIAL overlap
				    if ($sortedDomsInter[$j] =~ /TO--(\d+)---LENGTH/) { $toJ = $1; }
				    if ( ($fromI >= ($limitForSubseq-$varExtr)) or ($toJ <= ($overlapBetwSubseqs+$varExtr)) ) {
					my $resp = verifyOverlap ($limitForSubseq, $sortedDomsInter[$i],$sortedDomsInter[$j]);
					if    ($resp eq "II") { $sortedDomsInter[$i] .= "-->OUT"; }
					elsif ($resp eq "JJ") { $sortedDomsInter[$j] .= "-->OUT"; }
				    }
				} # IF ($sortedDomsInter[$j] !~ /-->OUT/)
			   }
			} # IF ($toI > $limitForSubseq)
		      } # IF ($sortedDomsInter[$i] ...
		   } # FOR (my $i = $fstIndx; $i < $secIndx-1; $i++)
		} # IF ($predInOverlReg eq "yes")
		$fstIndx = $secIndx;

	} # WHILE ($fstIndx < $qtt)
##    END of excluding redundancy for subsequent subsequences

##    mapping all predictions
	my $startSubseq;
	my @candidSuperfam = (); # ALL preds: mapped and without redundancies (with each ID and its domains)
	for (my $i = 0; $i < $qtt; $i++) {
		$line = $sortedDomsInter[$i];
		if ($line =~ /SUBSEQUENCE_(\d+)/) {
			my $numbSubseq  = $1;
			$startSubseq = ($numbSubseq - 1) * $limitForSubseq;
		}
		elsif ( ($line !~ /-->OUT/) and ($line =~ /DOMAIN--(.*)---FROM--(\d+)---TO--(\d+)---LENGTH--(.*)/) ) {
			my $domMap   = $1;
			my $fromMap  = $2 + $startSubseq;
			my $toMap    = $3 + $startSubseq;
			my $restLine = $4;
			my $linePrint = "DOMAIN--$domMap---FROM--$fromMap---TO--$toMap---LENGTH--$restLine";
			push (@candidSuperfam, {line => $linePrint, from => $fromMap});
		}
	}

	my @auxSorted = sort{$a->{from} <=> $b->{from}} @candidSuperfam;
	foreach my $lineAux (@auxSorted) { print PREDNO $lineAux->{line}."\n"; }

} # IF ($usingSubseqs eq "yes")
else { # for a whole sequence ($usingSubseqs eq "no")
	foreach $line (@sortedDomsInter) { print PREDNO "$line\n"; }
}

close (ALLPRED);
close (PREDNO);

