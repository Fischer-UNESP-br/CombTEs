# Developed by Carlos Fischer - 20 Oct 2020
# Updated: 25 Aug 2022 

# CombTEs is used to produce the final combinations considering TE candidates from tools used to search for TEs, for all "TE types" ("class", "order", "superfamily", etc.) of interest, assigning a TE classification for each final combination (called here the "FINAL GROUP").
# CombTEs.pl launches its associated scripts that generate a separate file for each used tool with its TE candidates.

# Usage: perl CombTEs.pl


# CombTEs can deal with TE candidates from any tool, since its final results (CANDIDATEs) are in the format:
# CANDIDATE_1 - FROM: pos_INIT - TO: pos_FINAL - LENGTH: length - SENSE: Direct/Reverse - CLASSIFICATION: TEtype_A
# PREDIC---FROM--pos_from---TO--anyPos---LENGTH--length_1---METRICS--metr---SENSE--D/R---TETYPE--TEtype_A
# PREDIC---FROM--anyPos---TO--anyPos---LENGTH--length_2---METRICS--metr---SENSE--D/R---TETYPE--TEtype_A
# ...
# PREDIC---FROM--anyPos---TO--anyPos---LENGTH--length_N---METRICS--metr---SENSE--D/R---TETYPE--TEtype_A
# (one blank line) (PAY ATTENTION HERE)*********
# CANDIDATE_2 - FROM: pos_INIT - TO: pos_FINAL - LENGTH: length - SENSE: Direct/Reverse - CLASSIFICATION: TEtype_B
# PREDIC---FROM--pos_from---TO--anyPos---LENGTH--length_1---METRICS--metr---SENSE--D/R---TETYPE--TEtype_B
# ...
# ...

# "CANDIDATE_1", "CANDIDATE_2", etc. must be sorted by their "FROM: " positions ("pos_INIT") inside each query sequence; the value of "METRICS" ("metr") will be used when filtering the predictions before producing the final candidates of each used tool.

# The predictions ("PREDIC") that make up each "CANDIDATE" must be sorted using their "FROM--" positions; then, the "FROM--" position ("pos_from") of the first prediction (i.e., the lowest value among the "FROM--" positions of all "PREDIC") of a same CANDIDATE will be the "FROM: " position ("pos_INIT") of that CANDIDATE. For the "TO: " position ("pos_FINAL") of a CANDIDATE, select the highest value among the "TO--" position ("anyPos") of all "PREDIC".


# WARNING:
# This script was implemented in a way that it can deal directly with TE predictions from the RepeatMasker, HMMER, and RPS-Blast programs (because they were used to test CombTEs and, then, maintained here). However, it can deal with results from ANY tools, and NO MATTER if those 3 programs are or are not considered in the analyses.

# If any one of those 3 programs is not considered, just exclude it in the "ParamsGeneral.pm" package (file), changing the "@tools" parameter - see instructions in that package.

# If another program is used, just include it in that package ("ParamsGeneral.pm"), changing the following 4 parameters there: "@tools", "%outFileNames", "%filterTools", and "%goodMetrTools" (see explanations/instructions there).
# ALSO: include a new "test condition" here in this script - search for "OBS_A" and "OBS_B" below.

# This script uses 6 parameters, which can be changed in the "ParamsGeneral.pm" package:

###########################################################################################

use strict;
use warnings;

use Cwd qw(getcwd);
use lib getcwd();

use ParamsGeneral qw(@TEtypes @tools @outFileNames %filterTools %goodMetrTools $distBetweenCands);

    my @classTEtype = (@TEtypes, "inconclusive"); # for the final classification of a TE candidate - the "inconclusive" classification refers to when the "integrase" domain is the only one domain found in a sequence (for example, when using RPS-Blast in the analyses)

    my @combClassif = (@TEtypes, "inconclusive"); # for the combinations between "types" of "@TEtypes"- a region can receive a classification like (for superfamily level) "Bel" or "Bel/Copia" or "Copia/Gypsy" -  used to produce the extra output info

###########################################################################################

    my $usedRpsblast = "no";
    foreach my $tool (@tools) { if ($tool eq "RpsBlast") { $usedRpsblast = "yes"; } }

    my $finalFile = "FinalClassification.txt";
    open (CLASS, ">$finalFile") or die "Can't open $finalFile";

    my $qttTools = scalar(@tools);
    print CLASS "Final classification for \"$tools[0]\"";
    for (my $i = 1; $i < $qttTools; $i++) { print CLASS ", \"$tools[$i]\""; }
    print CLASS ".\n";

    print CLASS "Used thresholds: **$tools[0] = $filterTools{$tools[0]}**";
    for (my $i = 1; $i < $qttTools; $i++) { print CLASS ", **$tools[$i] = $filterTools{$tools[$i]}**"; }
    print CLASS ".\n";
    print CLASS "Maximum distance from a tool's candidate to the next one to consider them inside the same FINAL GROUP: $distBetweenCands.\n\n";

    my $finalFileTab = "FinalClassificationCSV.csv";
    open (CLASSTAB, ">$finalFileTab") or die "Can't open $finalFileTab";

    print CLASSTAB "ID SEQ,GROUP,FROM,TO,LENGTH,SENSE,TE CLASSIF,IDENTIFIED BY,";
    if ($usedRpsblast eq "yes") { print CLASSTAB "CONSERVED DOMAINS,"; }
    print CLASSTAB "GOOD METRICS: **$tools[0] = $goodMetrTools{$tools[0]}**";
    for (my $i = 1; $i < $qttTools; $i++) { print CLASSTAB " / **$tools[$i] = $goodMetrTools{$tools[$i]}**"; }
    print CLASSTAB "\n";

    my %filehandles = ();
    for (my $i = 0; $i < $qttTools; $i++) {
	my $fileName = $outFileNames[$i];
	my $fileFinders = "FinalClassification_$fileName.txt";
	open (my $fh, ">", "$fileFinders") or die "Can't open $fileFinders";
	$filehandles{$fileName} = $fh;
	print $fh "This file shows the groups produced by $fileName tool(s)\n\n";
    }

########### inserting ALL candidates (of ALL TE types) from ALL tools in "@allCandid" for FINAL combination ###########
    my @seqNames; # for the IDs of the sequences with at least one candidate
    my %hashCandids = ();
    foreach my $tool (@tools) {
	my @candidTEtype = ();

# To combine **TE candidates** from TE search tools, CombTEs needs those candidates to be in an appropriate format (as shown in the beginning of this script), for ALL TE types together.

# WHEN using HMMER OR RepeatMasker: the script "finalCandidsHmmerRM.pl" produces the **TE candidates** for each tool (the "finalCandidates_TOOL.txt" file). In the way "finalCandidsHmmerRM.pl" was implemented, it needs the **initial predictions** from each tool to be in another appropriate format, what can be obtained using "extractHmmerRM.pl", which produces the "TEtype_TOOL.pred" file.
# WHEN using RPS-Blast: the same as for HMMER/RepeatMasker but using "finalCandidsRpsBlast.pl" (to produce the "finalCandidates_RpsBlast.txt" file) and "extractRPSB.pl" (to produce the "ConservedDomains_NoRedund.pred" file, with the identified domain predictions).

####
# The important thing here is the format of the file with the **TE candidates** from each used tool. The names of the output files and used scripts, of course, do not matter.
####

	if ( ($tool eq "HMMER") or ($tool eq "RepeatMasker") )	{ @candidTEtype = `perl finalCandidsHmmerRM.pl $tool`; }
	if ($tool eq "RpsBlast")				{ @candidTEtype = `perl finalCandidsRpsBlast.pl`; }

# OBS_A : for a different tool, just replicate the last line above, like:
#	if ($tool eq "otherTOOL1") { @candidTEtype = `perl scriptOtherTOOL1`; }
#	if ($tool eq "otherTOOL2") { @candidTEtype = `python scriptOtherTOOL2`; }
# just exclude the "comment"("#") in the beginning of the line above.
# For each different tool ("otherTOOL1", "otherTOOL2", etc.), the related script ("scriptOtherTOOL1", "scriptOtherTOOL2", etc.) should produce the TE candidates from that tool.
# Just use a similar line (command) for each different tool - above, "perl scriptOtherTOOL1" would be used to run this new script, when written in Perl; for a script implemented in another language, just replace "perl scriptOtherTOOL1" with the appropriate way to run it (for example, "python scriptOtherTOOL2").

	my ($from, $idSeqHash);
	foreach my $line (@candidTEtype) {
	    if ($line =~ /SEQUENCE: (.*)/) {
		$idSeqHash = $1;
		if ( !(exists ($hashCandids{$idSeqHash})) ) {
		    push (@seqNames, $idSeqHash);
		    $hashCandids{$idSeqHash} = "";
	    	}
	    }
	    elsif ($line =~ /CANDIDATE|###/) {
		if ($hashCandids{$idSeqHash} ne "") { $hashCandids{$idSeqHash} .= "###"; } # to separate the candidates
		chomp $line;
		if ($line !~ /###/) { $hashCandids{$idSeqHash} .= "    $line\t\t\t$tool\n"; }
	    }
	    elsif ($line =~ /PREDIC|DOMAIN/) { # patterns used in the output file from HMMER/RepeatMasker or RPS-Blast when using the cited scripts (above) to parse the initial results from these tools (the "extractHmmerRM.pl" and "extractRPSB.pl" scripts)
		$hashCandids{$idSeqHash} .= "\t$line";
	    }
	} # FOREACH my $line (@candidTEtype)
    } # FOREACH my $tool (@tools)

####
# Here, after the FOREACH above ("foreach my $tool (@tools)"), for each analised sequence (with "$idSeqHash" ID), the related "$hashCandids{$idSeqHash}" hash will contain ALL the candidates from ALL the used tools
####

    my @allCandid = (); # to be used to sort (using the start position of) the candidates inside a same sequence ("$idSeqHash"). The candidates are separated by a "###"
    my %sameSeq = ();
    foreach my $idSeqHash (@seqNames) {
	$sameSeq{"general"} = "no";
	for (my $i = 0; $i < $qttTools; $i++) {
		my $fileName = $outFileNames[$i];
		$sameSeq{$fileName} = "no";
	}

	my @arraySplit = split '###', $hashCandids{$idSeqHash};
	my $qttSplit = scalar(@arraySplit);
	for (my $j = 0; $j < $qttSplit; $j++) {
	    my $lineCand = $arraySplit[$j];
	    if ($lineCand =~ /CANDIDATE_\d+ - FROM: (\d+) - TO/) {
		my $startCandid = $1;
		push (@allCandid, {line => $lineCand, from => $startCandid} );
	    }
	}

	my @auxSorted = sort{ $a->{from} <=> $b->{from} } @allCandid;
	@allCandid = ();
	my @sortedAllCandid = ();
	foreach my $lineAux (@auxSorted) { push (@sortedAllCandid, $lineAux->{line}); }
	my $qttCandid = scalar(@sortedAllCandid);

########### producing the combined groups ###########
	my %classif = ();
	foreach my $clte (@classTEtype) { $classif{$clte} = 0; }

	my %toolFinder = ();
	foreach my $tool (@tools) { $toolFinder{$tool} = 0; }

	my %typeClassFinal = ();
	for (my $i = 0; $i < $qttTools; $i++) {
		my $fileName = $outFileNames[$i];
		foreach my $comb (@combClassif) { $typeClassFinal{$fileName}{$comb} = 0; }
	}

	my $group      = 0;
	my @arrayGroup = ();
	my $fromGroup  = 0;
	my $toGroup    = 0;
	my $newGroup   = "yes";
	my $i = 0;
	my ($from1, $to1, $sense1, $senseGroup, $from2, $sense2);
	while ($i < $qttCandid) {
		push (@arrayGroup, $sortedAllCandid[$i]);

### the format of each final TE candidate is:
#	CANDIDATE_numb - FROM: ff - TO: tt - LENGTH: ll - SENSE: d/r - CLASSIFICATION: type \n 
#	PREDIC---FROM--ff--- ...	(OR)
#	DOMAIN--dom---FROM-- ...
		if ($sortedAllCandid[$i] =~ /FROM: (\d+) - TO: (\d+) - LENGTH: \d+ - SENSE: (\w+) - CLASSIFICATION: (\w+)\t\t\t(\w+)\n\t(PREDIC|DOMAIN)/) {
			$from1  = $1;
			$to1    = $2;
			$sense1 = $3;
			my $class1 = $4;
			$classif{$class1} = 1;
			my $tool1  = $5;
			$toolFinder{$tool1} = 1;
		}

		if ($newGroup eq "yes") {
			$fromGroup  = $from1;
			$senseGroup = $sense1;
			$newGroup   = "no";
		}

		if ($to1 > $toGroup) { $toGroup = $to1; }

		if ( ($i+1) < $qttCandid ) { # there is a next candidate
			if ($sortedAllCandid[$i+1] =~ /FROM: (\d+) - TO: .* - SENSE: (.*) - CLASSIFICATION/) {
				$from2  = $1;
				$sense2 = $2;
			}
		}
		if ( ($i+1) == $qttCandid ) { $newGroup   = "yes"; }
		elsif ( (($from2-$toGroup) > $distBetweenCands) or ($sense2 ne $senseGroup) ) { $newGroup   = "yes"; }

		if ($newGroup eq "yes") {
# if (current CANDIDATE is the last one) OR (distance between next candidate and current GROUP is higher than "$distBetweenCands") OR (they are in different senses), then a new GROUP must be formed: write the current GROUP in the output files.
			my $finalClass = "";
			my $finalTools = "";
			$group++;
			my $lenGroup = $toGroup - $fromGroup + 1;

			foreach my $clte (@TEtypes) {
				if ($classif{$clte} == 1) {
					if ($finalClass eq "")  { $finalClass = $clte; }
					else			{ $finalClass .= "/$clte"; }
					$classif{$clte} = 0;
				}
			}
			if ($classif{"inconclusive"} == 1) {
				if ($finalClass eq "") { $finalClass = "inconclusive"; }
				$classif{"inconclusive"} = 0;
			}

			my $numTools = 0;
			foreach my $tool (@tools) {
				if ($toolFinder{$tool} == 1) {
					if ($finalTools eq "")  { $finalTools = $tool; }
					else			{ $finalTools .= " / $tool"; }
					$numTools++;
					$toolFinder{$tool} = 0;
				}
			}

			my $indx = $numTools - 1;
			my $fileName = $outFileNames[$indx];
			my $FILEtoWRITE = $filehandles{$fileName};

			my $idPrint = ">>>SEQUENCE: $idSeqHash\n";
			if ($sameSeq{"general"} eq "no") {
				print CLASS "$idPrint";
				print CLASSTAB "\n$idSeqHash";
				$sameSeq{"general"} = "yes";
			}
			if ($sameSeq{$fileName} eq "no") {
				print $FILEtoWRITE $idPrint;
				$sameSeq{$fileName} = "yes";
			}

			my $descrGroup = "GROUP $group - FROM: $fromGroup - TO: $toGroup - LENGTH: $lenGroup - SENSE: $senseGroup - CLASS: $finalClass - TOOL(S): $finalTools\n";
			print CLASS $descrGroup;
			print $FILEtoWRITE $descrGroup;

			my $descrGroupTab = ",$group,$fromGroup,$toGroup,$lenGroup,$senseGroup,$finalClass,$finalTools,";

			my %toolsGoodMetr = (); # specification of the tools with "at least one GOOD metrics"
			foreach my $tool (@tools) { $toolsGoodMetr{$tool} = 0; }

			my $foundDomains = "";
			my $qttArrayGroup = scalar(@arrayGroup);
			for (my $j = 0; $j < $qttArrayGroup; $j++) {
			   print CLASS "$arrayGroup[$j]\n";
			   print $FILEtoWRITE "$arrayGroup[$j]\n";

			   my @arraySplit = split '\n',$arrayGroup[$j];
			   my $qttSepara = scalar(@arraySplit);
			   for (my $jj = 0; $jj < $qttSepara; $jj++) {
				if    ($arraySplit[$jj] =~ /PREDIC---.*---EVALUE--(.*)---SCORE/) {
					if ($1 <= $goodMetrTools{'HMMER'}) { $toolsGoodMetr{"HMMER"} = 1; }
				}
				elsif ($arraySplit[$jj] =~ /RMSCORE--(\d+)---SENSE/) {
					if ($1 >= $goodMetrTools{'RepeatMasker'}) { $toolsGoodMetr{"RepeatMasker"} = 1; }
				}
				elsif ($arraySplit[$jj] =~ /DOMAIN--(.*)---FROM--.*---EVALUE--(.*)---SCORE/) {
					if ($foundDomains eq "") { $foundDomains = $1; }
					else			 { $foundDomains .= " / $1"; }
					if ($2 <= $goodMetrTools{'RpsBlast'}) { $toolsGoodMetr{"RpsBlast"} = 1; }
				}

#				elsif ($arraySplit[$jj] =~ /"regex"/) {
#					if ($1 >= $goodMetrTools{'otherTOOL'}) { $toolsGoodMetr{"otherTOOL"} = 1;}
#				}
# OBS_B : IF A DIFFERENT TOOL ("otherTOOL") was considered in the searches for TEs, a new test ("ELSIF") MUST BE INCLUDE HERE; this new test should filter (using "$goodMetrTools{'otherTOOL'}") the candidates based on the metrics reported by this new "otherTOOL".
# Just delete the comments ("#") in the beginning of the 3 lines above and prepare the condition test as for HMMER, for example - it is not necessary to do anything else here, unless you know what you are doing!!!

			   } # FOR (my $jj = 0; $jj < $qttSepara; $jj++)
			} # FOR (my $j = 0; $j < $qttArrayGroup; $j++)

			print CLASS "\n";
			print $FILEtoWRITE "\n";

			my $toolsWithGoodMetrics = "";
			foreach my $tool (@tools) {
				if ($toolsGoodMetr{$tool} == 1) {
					if ($toolsWithGoodMetrics eq "") { $toolsWithGoodMetrics = $tool; }
					else { $toolsWithGoodMetrics .= " / $tool"; }
				}
			}

			if ($foundDomains eq "") { $foundDomains = "No domain"; }

			if ($usedRpsblast eq "yes") { $descrGroupTab .= "$foundDomains,"; }
			$descrGroupTab .= "$toolsWithGoodMetrics\n";
			print CLASSTAB $descrGroupTab;

			@arrayGroup = ();
			$fromGroup  = 0;
			$toGroup    = 0;

			my $combInCombClassif = "no";
			foreach my $comb (@combClassif) { if ($finalClass eq $comb) { $combInCombClassif = "yes"; last; } }
			if ($combInCombClassif eq "no") {
				$typeClassFinal{$fileName}{$finalClass} = 0;
				push (@combClassif, $finalClass); # to include a new "TEtype combination" in "@combClassif"
			}
			$typeClassFinal{$fileName}{$finalClass}++;
		} # IF ($newGroup eq "yes")

		$i++;
	} # WHILE ($i < $qttCandid)
########### END of producing the combined groups ###########

###  extra output info:
	print CLASS "##########\n\n";
	print CLASS "SUMMARY OF THE GENERATED GROUPS:\n\n";

	for (my $i = 0; $i < $qttTools; $i++) {
		my $fileName = $outFileNames[$i];
		my $qttGroup = 0;

		foreach my $comb (@combClassif) {
			if ( exists $typeClassFinal{$fileName}{$comb} ) { $qttGroup += $typeClassFinal{$fileName}{$comb}; }
		}
		print CLASS "Groups generated by $fileName tool(s): $qttGroup\n";
		print CLASS "   Final classification:\n";
		foreach my $comb (@combClassif) {
			if ( (exists $typeClassFinal{$fileName}{$comb}) and ($typeClassFinal{$fileName}{$comb} != 0) ) {
				print CLASS "      - $comb: $typeClassFinal{$fileName}{$comb}\n";
			}
		}
		print CLASS "\n";
	}
	print CLASS "######################################################################################################\n\n";

    } # FOREACH my $idSeqHash (@seqNames)

    for (my $i = 0; $i < $qttTools; $i++) {
	my $fileName = $outFileNames[$i];
	close ($filehandles{$fileName});
    }
    close (CLASS);
    close (CLASSTAB);

