# Developed by Carlos Fischer - 20 Oct 2020

# To extract the predictions from the output files of HMMER and RepeatMasker, excluding redundancies between predictions.
# This script can be used to extract info of interest when analysing a whole long sequence or subsequences of a long sequence.

# Usage: perl extractHmmerRM.pl TOOL SUPERFAM file-superfam-tool

# input ("file-superfam-tool"): file with results from HMMER or RepeatMasker for each superfamily
# output: two files - one with ALL predictions and other with predictions without redundancies, for EACH superfamily


# THIS SCRIPT WOULD BE RUN JUST ONCE FOR EACH SUPERFAMILY
# UNLESS the user changes the values of some parameters used here.

# This script uses 7 parameters:
# - 4 of them can be changed in "ParamsGeneral.pm":
#	- @superfamilies: the used superfamilies;
#	- $usingSubseqs: when using subsequences of a long sequence ($usingSubseqs = "yes") or directly on a whole long sequence 		($usingSubseqs = "no");
# 	- $overlap: length of the overlap region;
#	- $limit: equals to "$lengthSubseq" (length of the subsequences) minus "$overlap" (length of the overlap region) (see 		"ParamsGeneral.pm");
# - 3 others can be changed in "ParamsHmmerRM.pm":
#	- $varExtr: to allow small variation between respective extremities for cases of overlapping between 2 predictions (default= 		50).
#	- $minOverlPred: minimum value to consider that there is overlapping between two predictions (default = 100 - a higher value 		would increase the number of predictions to be considered in the following analysis);
#	- "$minLenPred": minimum length for a prediction to be included in the analysis (default = 100).


# For the case of subsequences of a long sequence, the IDs of each subsequence must be in the format "(XX)Subsequence-NumSubseq": "XX" means anything, or even nothing, and "NumSubseq" is a number in ascending order for each subsequence.
# Also available, there is a script ("splitOverlap.pl") that splits a long sequence in small ones with an overlapping region between 2 subsequent subsequences.

###########################################################################################

use strict;
use warnings;

use Cwd qw(getcwd);
use lib getcwd(); 

use ParamsGeneral qw(@superfamilies $usingSubseqs $overlap $limit);
use ParamsHmmerRM qw($varExtr $minOverlPred $minLenPred);

###########################################################################################

sub excludeRedund { # to exclude redundant predictions

	(my $method, my @predictions) = @_;

	my $qttR = scalar(@predictions);
	my @auxPreds = ();
	for (my $i = 0; $i < $qttR; $i++) {
#					 FROM--fff---TO--ttt---LENGTH
		if ($predictions[$i] =~ /FROM--(\d+)---TO--(\d+)---LENGTH/) {
			push (@auxPreds, {line => $predictions[$i], from => $1, to => $2});
		}
	}
	my @sortedPredicts = sort{$a->{from} <=> $b->{from}} @auxPreds;

	for (my $i = 0; $i < $qttR-1; $i++) { # comparing predictions all against all
		for (my $j = $i+1; $j < $qttR; $j++) {
			if ($sortedPredicts[$i]->{line} =~ /-->OUT/) { last; }
			if ($sortedPredicts[$j]->{line} !~ /-->OUT/) {
				my $resp = verifyOverlap ($method, 0, $sortedPredicts[$i]->{line},$sortedPredicts[$j]->{line});
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


sub verifyOverlap { # to verify overlap between 2 predictions and, if so, retrieve an ID of the one to be disregarded; this prediction would be:
# - for HMMER:        the one with the highest e-value or, for the same e-values, the shortest prediction
# - for RepeatMasker: the one with the lowest SWscore or,  for the same SWscores, the shortest prediction.

	(my $method, my $offset, my $seq1, my $seq2) = @_;

	my ($from1, $to1, $length1, $evalSWs1, $sense1);
	my ($from2, $to2, $length2, $evalSWs2, $sense2);

# For HMMER:
# $seq= PREDIC---FROM--ff---TO--tt---LENGTH--ll---EVALUE--ev---SCORE--sc---SENSE--dr---SUPERFAM--spfam
# For RepeatMasker:
# $seq= PREDIC---FROM--ff---TO--tt---LENGTH--ll---SWSCORE--sc---SENSE--dr---MATCHINGREPEAT--match---SUPERFAM--spfam
	if (
	     ($seq1 =~ /^PREDIC---FROM--(.*)---TO--(.*)---LENGTH--(.*)---EVALUE--(.*)---SCORE--.*---SENSE--(.*)---SUPERFAM/)    or
	     ($seq1 =~ /^PREDIC---FROM--(.*)---TO--(.*)---LENGTH--(.*)---SWSCORE--(.*)---SENSE--(.*)---MATCHINGREPEAT/)
	   ) {
		$from1    = $1;
		$to1      = $2;
		$length1  = $3;
		$evalSWs1 = $4;
		$sense1   = $5;
	}
	if (
	     ($seq2 =~ /^PREDIC---FROM--(.*)---TO--(.*)---LENGTH--(.*)---EVALUE--(.*)---SCORE--.*---SENSE--(.*)---SUPERFAM/)    or
	     ($seq2 =~ /^PREDIC---FROM--(.*)---TO--(.*)---LENGTH--(.*)---SWSCORE--(.*)---SENSE--(.*)---MATCHINGREPEAT/)
	   ) {
		$from2    = $1 + $offset;
		$to2      = $2 + $offset;
		$length2  = $3;
		$evalSWs2 = $4;
		$sense2   = $5;
	}

	my $respOverl = "NONE"; # there is NOT overlap between seq_1 and seq_2 OR retrieve an ID of the prediction to be disregarded
	if ($sense1 eq $sense2) {
	    if ( ($from2 <= $to1) or ($to2 >= $from1) ) {
		if (
		    (abs($from1-$from2) <= $varExtr) or (abs($to1-$to2) <= $varExtr) or # small variation in corresponding extremities
		    (($from1 < $from2) and ($to2 <  $to1)) or				# seq_2 inside seq_1
		    (($from2 < $from1) and ($to1 <  $to2)) or				# seq_1 inside seq_2
		    (($from1 < $from2) and (($to1 - $from2) >= $minOverlPred)) or	# overlap >= $minOverlPred
		    (($to2   < $to1)   and (($to2 - $from1) >= $minOverlPred))
		   ) { # there is overlap between seq_1 and seq_2
			if ($method eq "HMMER") {
				if    ($evalSWs1 < $evalSWs2) { $respOverl = "JJ"; }
				elsif ($evalSWs1 > $evalSWs2) { $respOverl = "II"; }
				else {
					if ($length1 >= $length2) { $respOverl = "JJ"; }	
					else			  { $respOverl = "II"; }
				} 
			}
			elsif ($method eq "RepeatMasker") {
				if ( ($length1 >= $minLenPred) and ($length2 >= $minLenPred) ) { 
					if    ($evalSWs1 > $evalSWs2) { $respOverl = "JJ"; }
					elsif ($evalSWs1 < $evalSWs2) { $respOverl = "II"; }
					else {
						if ($length1 >= $length2) { $respOverl = "JJ"; }	
						else			  { $respOverl = "II"; }
					}
				}
				else {
					if ($length1 >= $length2) { $respOverl = "JJ"; }	
					else			  { $respOverl = "II"; }
				}
			}
		} # IF ( (abs($from1-$from2) <= $varExtr) or ....
	    } # IF ( ($from2 < $to1) or ($to2 > $from1) )
	} # IF ($sense1 eq $sense2)

	return "$respOverl";
} # END of SUB "verifyOverlap"

###########################################################################################


unless(@ARGV) { die "\nUSAGE: perl $0 TOOL SUPERFAM file_for_superfam\n\n"; }

my $tool    = $ARGV[0];
my $spfamIN = $ARGV[1];
my $fileIN  = $ARGV[2];

if ( ($tool ne "HMMER") and ($tool ne "RepeatMasker") ) { die "\nWrong TOOL!!! Use \"HMMER\" or \"RepeatMasker\".\n\n"; }

my $OKspfam = 0;
foreach my $spfam (@superfamilies) { if ($spfam eq $spfamIN) { $OKspfam = 1; } }
if ($OKspfam == 0) { die "\nWrong SUPERFAMILY!!!\n\n"; }
open (INPRED, $fileIN) or die "\nCan't open \"$fileIN\"!!!\n\n";

my $outPred = "$spfamIN\_$tool\_ALL.pred";
open (OUTPRED,">$outPred") or die "Can't open $outPred!\n";
print OUTPRED "All predictions of $tool for superfamily \"$spfamIN\", from file \"$fileIN\".\n\n";

my $predNoRedund = "$spfamIN\_$tool.pred";
open (PREDNO, ">$predNoRedund") or die "Can't open $predNoRedund";
print PREDNO "Predictions of $tool for superfamily \"$spfamIN\", from file \"$fileIN\".\n";
print PREDNO "Minimum value to consider overlapping between two predictions: $minOverlPred.\n\n";

my @allPredicts   = ();

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
				my $aux   = $from;
				$from  = $to;
				$to    = $aux;
			}
			$length = $to - $from + 1;

			$lineOUT = "PREDIC---FROM--$from---TO--$to---LENGTH--$length---EVALUE--$evalue---SCORE--$score---SENSE--$sense---SUPERFAM--$spfamIN";
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

			$lineOUT = "PREDIC---FROM--$from---TO--$to---LENGTH--$length---SWSCORE--$score---SENSE--$sense---MATCHINGREPEAT--$matchRM---SUPERFAM--$spfamIN";
		}
		my $numSubseq;
		if ($idSeq =~ /Subseq-(\d+)/) { $numSubseq = $1; }	# for subsequences of a long sequence
		elsif ($usingSubseqs eq "no") { $numSubseq = $idSeq; }	# for a whole long sequence
		else { die "ERROR: formatting problems of the input file!!!";}

		if    ($usingSubseqs eq "yes") { push (@allPredicts, {line => $lineOUT, numSubseq => $numSubseq}); }
		elsif ($usingSubseqs eq "no")  { push (@allPredicts, {line => $lineOUT, from => $from}); }
	} # IF ($line !~ /^#/)

} # WHILE (not eof INPRED)
close (INPRED);

my $qttPred = scalar(@allPredicts);
if ($qttPred == 0) { die "\nInput file is empty!!!\n\n"; }

if ($usingSubseqs eq "yes") { # ONLY for subsequences from a long sequence
	my @sortedAllPredicts = sort{$a->{numSubseq} <=> $b->{numSubseq}}@allPredicts;

	$qttPred = scalar(@sortedAllPredicts);
	my @predSubAll = ();
	my @predExclud = ();
	my $countPred  = 0;

	my @arrayAux = (); # for each ID and its predictions (for ALL subseqs together, when analysing subseqs of a long sequence)
	my $currNumSubseq = 0;
	foreach my $lineSorted (@sortedAllPredicts) {
		$countPred++;

		if ($lineSorted->{numSubseq} != $currNumSubseq) {
			if ($currNumSubseq != 0) {
				my @predSubAllSorted = sort{$a->{from} <=> $b->{from}} @predSubAll;
				foreach my $lineAux (@predSubAllSorted) { print OUTPRED $lineAux->{line}."\n"; }
				print OUTPRED "###\n";
				@predSubAll = ();
#			    excluding redundancies between predictions of interest
				my @sortedStructAux = excludeRedund ($tool, @predExclud);
				foreach my $lineAux (@sortedStructAux) { push (@arrayAux, $lineAux); } 
				push (@arrayAux, "###");
				@predExclud = ();
			}
			$currNumSubseq = $lineSorted->{numSubseq};
			print OUTPRED    ">>>SUBSEQUENCE_$currNumSubseq\n";
			push (@arrayAux, ">>>SUBSEQUENCE_$currNumSubseq");
		}

		my $lineOUT = $lineSorted->{line};
		if ($lineOUT =~ /FROM--(\d+)---TO/) { push (@predSubAll, {line => $lineOUT, from => $1}); }
		push (@predExclud, $lineOUT);

		if ($countPred == $qttPred) {
			my @predSubAllSorted = sort{$a->{from} <=> $b->{from}} @predSubAll;
			foreach my $lineAux (@predSubAllSorted) { print OUTPRED $lineAux->{line}."\n"; }
			print OUTPRED "###\n";
			my @sortedStructAux = excludeRedund ($tool, @predExclud);
			foreach my $lineAux (@sortedStructAux) { push (@arrayAux, $lineAux); } 
			push (@arrayAux, "###");
		}
 	} # FOREACH $lineSorted (@sortedAllPredicts)

##    excluding redundancy between subsequences (predictions into an overlapping region)
	my $qtt = scalar(@arrayAux);
	my $fstIndx = 0;
	while ($fstIndx < $qtt) {
		my $predInOverlReg = "no";
		my $secIndx = $fstIndx + 1;
		$line = $arrayAux[$secIndx];
#				>>>SUBSEQUENCE_numSubseq
		while ( ($line !~ /SUBSEQUENCE_(\d+)/) and ($secIndx < $qtt) ) {
			# firstly, check if any pred is into overlapping region
			if ($line =~ /TO--(\d+)---LENGTH/) { if ($1 > $limit) { $predInOverlReg = "yes"; } }
			$secIndx++;
			if ($secIndx < $qtt) { $line = $arrayAux[$secIndx]; }
		}

		my $trdIndx = $secIndx + 1;
		if ( ($secIndx < $qtt) and ($predInOverlReg eq "yes") ) {
			$line = $arrayAux[$trdIndx];
			while ( ($line !~ /SUBSEQUENCE_(\d+)/) and ($trdIndx < $qtt) ) {
				$trdIndx++;
				if ($trdIndx < $qtt) { $line = $arrayAux[$trdIndx]; }
			}
		}

		if ($predInOverlReg eq "yes") { # checking possible overlapping between predictions of subsequent subsequences
		   for (my $i = ($fstIndx+1); $i < $secIndx; $i++) { # compare predictions from 1 subseq to all preds from next subseq
		      if ($arrayAux[$i] =~ /FROM--(\d+)---TO--(\d+)---LENGTH/) {
			my $fromI = $1;
			my $toI   = $2;
			if ($toI > $limit) {
			   my $toJ;
			   for (my $j = ($secIndx+1); $j < $trdIndx-1; $j++) {
				if ($arrayAux[$i] =~ /-->OUT/) { last; }
				if ($arrayAux[$j] !~ /-->OUT/) { # for COMPLETE or PARTIAL overlapping
				    if ($arrayAux[$j] =~ /TO--(\d+)---LENGTH/) { $toJ = $1; }
				    if ( ($fromI >= ($limit-$varExtr)) or ($toJ <= ($overlap+$varExtr)) ) {
					my $resp = verifyOverlap ($tool, $limit, $arrayAux[$i], $arrayAux[$j]);
					if    ($resp eq "II") { $arrayAux[$i] .= "-->OUT"; }
					elsif ($resp eq "JJ") { $arrayAux[$j] .= "-->OUT"; }
				    }
				} # IF ($arrayAux[$j] !~ /-->OUT/)
			   }
			} # IF ($toI > $limit)
		      } # IF ($arrayAux[$i] ...
		   } # FOR (my $i = $fstIndx; $i < $secIndx-1; $i++)
		} # IF ($predInOverlReg eq "yes")
		$fstIndx = $secIndx;

	} # WHILE ($fstIndx < $qtt)
##    END of excluding redundancy between subsequences

##    mapping each prediction of each superfamily
	my $startSubseq;
	my @candidSuperfam = ();
	for (my $i = 0; $i < $qtt; $i++) {
		$line = $arrayAux[$i];
		if ($line =~ /SUBSEQUENCE_(\d+)/) {
			my $numbSubseq  = $1;
			$startSubseq = ($numbSubseq - 1) * $limit;
		}
		elsif ( ($line !~ /-->OUT/) and ($line =~ /FROM--(\d+)---TO--(\d+)---LENGTH--(.*)/) ) {
			my $fromMap  = $1 + $startSubseq;
			my $toMap    = $2 + $startSubseq;
			my $restLine = $3;
			my $linePrint = "PREDIC---FROM--$fromMap---TO--$toMap---LENGTH--$restLine";
			push (@candidSuperfam, {line => $linePrint, from => $fromMap});
		}
	} # FOR (my $i = 0; $i < $qtt; $i++)

	my @auxSorted = sort{$a->{from} <=> $b->{from}} @candidSuperfam;
	foreach my $lineAux (@auxSorted) { print PREDNO $lineAux->{line}."\n"; }

} # IF ($usingSubseqs eq "yes")
else { # for a whole sequence ($usingSubseqs eq "no")
	my @sortedAllPredicts = sort{$a->{from} <=> $b->{from}} @allPredicts;
	print OUTPRED ">>>SEQUENCE: $idSeq\n";

	my @sortedAllPreds = ();
	foreach my $lineAux (@sortedAllPredicts) {
		print OUTPRED $lineAux->{line}."\n";
		push (@sortedAllPreds, $lineAux->{line});
	}
	print OUTPRED "###\n";

	my @sortedStructAux = excludeRedund ($tool, @sortedAllPreds); # excluding redundancies between predictions of interest
	foreach my $lineAux (@sortedStructAux) { print PREDNO "$lineAux\n"; } 
}

close (OUTPRED);
close (PREDNO);

