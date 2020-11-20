# Developed by Carlos Fischer - 20 Oct 2020

# To extract the predictions from the output files of RPS-Blast, excluding redundancies between predictions
# This script can be used to extract info of interest when analysing a whole long sequence or subsequences of a long sequence.

# Usage: perl extractRPSB.pl output_file_rpsblast

# input:  file with results from RPS-Blast
# output: two files - one with ALL predictions and other with predictions without redundancies, for ALL considered superfamilies

# It considers the use of "RPSTBLASTN 2.10.1+" (and other version which maintains the same output format)


# THIS SCRIPT WOULD BE RUN JUST ONCE
# UNLESS the user changes the values of some parameters used here.

# This script uses 7 parameters:
# - 4 of them can be changed in "ParamsGeneral.pm":
#	- @superfamilies: the used superfamilies;
#	- $usingSubseqs: when using subsequences of a long sequence ($usingSubseqs = "yes") or directly on a whole long sequence 		($usingSubseqs = "no");
# 	- $overlap: length of the overlap region;
#	- $limit: equals to "$lengthSubseq" (length of the subsequences) and "$overlap" (length of the overlap region) (see 		"ParamsGeneral.pm");
# - 3 others can be changed in "ParamsRpsblast.pm":
#	- $varExtr: to allow small variation between respective extremities for cases of overlapping between 2 predictions (default= 		50).
#	- $minOverlPred: minimum value to consider that there is overlapping between two predictions (default = 100 - a higher value 		would increase the number of predictions to be considered in the following analysis);
#	- %domains: conserved domains considered in the analysis (after running RPS-Blast), for:
#		- Bel:   DUF1759, Peptidase_A17, DUF1758, RT_pepA17;
#		- Copia: Retrotran_gag, gag_pre-integrs, RVT_2, RNase_HI_RT_Ty1;
#		- Gypsy: Retrotrans_gag, gag-asp_proteas, retropepsin_like, RP_Saci_like, RVP_2, RT_LTR, RVT_3, RNase_HI_RT_Ty3, 			RNase_HI_like;
#		- and:   rve (Integrase), for Bel, Copia, and Gypsy.


# For the case of subsequences of a long sequence, the IDs of each subsequence must be in the format "(XX)Subsequence-NumSubseq": "XX" means anything, or even nothing, and "NumSubseq" is a number in ascending order for each subsequence.
# Also available, there is a script ("splitOverlap.pl") that splits a long sequence in small ones with an overlapping region between 2 subsequent subsequences.

###########################################################################################

use strict;
use warnings;

use Cwd qw(getcwd);
use lib getcwd(); 

use ParamsGeneral  qw(@superfamilies $usingSubseqs $overlap $limit);
use ParamsRpsblast qw($varExtr $minOverlPred %domains);

###########################################################################################

sub excludeRedund { # to exclude redundant predictions

	(my @predictions) = @_;

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


sub verifyOverlap { # to verify overlap between 2 predictions and, if so, retrieve an ID of the one to be disregarded; this prediction would be the one with the highest e-value or, for the same e-values, the shortest prediction

	(my $offset, my $seq1, my $seq2) = @_;

	my ($from1, $to1, $length1, $eval1, $sense1);
	my ($from2, $to2, $length2, $eval2, $sense2);

#	$seq= DOMAIN--dom---FROM--ff---TO--tt---LENGTH--ll---EVALUE--ev---SCORE--sc---SENSE--dr
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
			if    ($eval1 < $eval2) { $respOverl = "JJ"; }
			elsif ($eval1 > $eval2) { $respOverl = "II"; }
			else {
				if ($length1 >= $length2) { $respOverl = "JJ"; }	
				else			  { $respOverl = "II"; }
			}
		}
	    } # IF ( ($from2 < $to1) or ($to2 > $from1) )
	} # IF ($sense1 eq $sense2)

	return "$respOverl";
} # END of SUB "verifyOverlap"

###########################################################################################


unless(@ARGV) { die "\nUSAGE: perl $0 file_from_RpsBlast\n\n"; }

my $fileRpsb = $ARGV[0];
open (INPRED, $fileRpsb) or die "Can't open $fileRpsb!\n";

my $allPred = "ConservedDomains_ALL-3.pred";
open (ALLPRED,">$allPred") or die "Can't open $allPred";
print ALLPRED "All predictions of RPSBlastn for conserved domains, from file \"$fileRpsb\".\n\n";

my $predNoRedund = "ConservedDomains_NoRedund-3.pred";
open (PREDNO, ">$predNoRedund") or die "Can't open $predNoRedund";
print PREDNO "Predictions of RPSBlastn for conserved domains of interest, from file \"$fileRpsb\".\n";
print PREDNO "Minimum value to consider overlapping between two predictions: $minOverlPred.\n\n";

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

	$line = readline(INPRED); # 3th line - # Database: 
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

#		selecting domains of interest
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
#		    excluding redundancies between predictions of interest
			my @sortedInter = excludeRedund (@domsInterest);
			if ($usingSubseqs eq "yes") { push (@sortedDomsInter, $seqID); }
			foreach my $lineAux (@sortedInter) { push (@sortedDomsInter, $lineAux); }
		}
	}
} # WHILE (not eof INPRED)
close (INPRED);

if ($usingSubseqs eq "yes") { # ONLY for subsequences from a long sequence
##    excluding redundancy between subsequences (predictions into an overlapping region)
	my $qtt = scalar(@sortedDomsInter);
	my $fstIndx = 0;
	while ($fstIndx < $qtt) {
		my $predInOverlReg = "no";
		my $secIndx = $fstIndx + 1;
		$line = $sortedDomsInter[$secIndx];
#				>>>SUBSEQUENCE_numSubseq
		while ( ($line !~ /SUBSEQUENCE_(\d+)/) and ($secIndx < $qtt) ) {
			# firstly, check if any domain prediction is into overlapping region
			if ($line =~ /TO--(\d+)---LENGTH/) { if ($1 > $limit) { $predInOverlReg = "yes"; } }
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

		if ($predInOverlReg eq "yes") { # checking possible overlapping between predictions of subsequent subsequences
		   for (my $i = ($fstIndx+1); $i < $secIndx; $i++) { # compare predictions from 1 subseq to all preds from next subseq
		      if ($sortedDomsInter[$i] =~ /FROM--(\d+)---TO--(\d+)---LENGTH/) {
			my $fromI = $1;
			my $toI   = $2;
			if ($toI > $limit) {
			   my $toJ;
			   for (my $j = ($secIndx+1); $j < $trdIndx; $j++) {
				if ($sortedDomsInter[$i] =~ /-->OUT/) { last; }
				if ($sortedDomsInter[$j] !~ /-->OUT/) { # for COMPLETE or PARTIAL overlapping
				    if ($sortedDomsInter[$j] =~ /TO--(\d+)---LENGTH/) { $toJ = $1; }
				    if ( ($fromI >= ($limit-$varExtr)) or ($toJ <= ($overlap+$varExtr)) ) {
					my $resp = verifyOverlap ($limit, $sortedDomsInter[$i],$sortedDomsInter[$j]);
					if    ($resp eq "II") { $sortedDomsInter[$i] .= "-->OUT"; }
					elsif ($resp eq "JJ") { $sortedDomsInter[$j] .= "-->OUT"; }
				    }
				} # IF ($sortedDomsInter[$j] !~ /-->OUT/)
			   }
			} # IF ($toI > $limit)
		      } # IF ($sortedDomsInter[$i] ...
		   } # FOR (my $i = $fstIndx; $i < $secIndx-1; $i++)
		} # IF ($predInOverlReg eq "yes")
		$fstIndx = $secIndx;

	} # WHILE ($fstIndx < $qtt)
##    END of excluding redundancy between subsequences

##    mapping all prediction
	my $startSubseq;
	my @candidSuperfam = (); # ALL preds: mapped and without redundancies (with each ID and its domains)
	for (my $i = 0; $i < $qtt; $i++) {
		$line = $sortedDomsInter[$i];
		if ($line =~ /SUBSEQUENCE_(\d+)/) {
			my $numbSubseq  = $1;
			$startSubseq = ($numbSubseq - 1) * $limit;
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

##    sorting predictions of "@candidSuperfam"
	my @auxSorted = sort{$a->{from} <=> $b->{from}} @candidSuperfam;
	foreach my $lineAux (@auxSorted) { print PREDNO $lineAux->{line}."\n"; }

} # IF ($usingSubseqs eq "yes")
else { # for a whole sequence ($usingSubseqs eq "no")
	foreach $line (@sortedDomsInter) { print PREDNO "$line\n"; }
}

close (ALLPRED);
close (PREDNO);

