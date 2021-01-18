# Developed by Carlos Fischer - 20 Oct 2020

# Used to combine the final candidates of all used tools, for all superfamilies of interest, and produce their final classification.

# Usage: perl CombTEs.pl

# This script was originally aimed to process results and produce TE candidates of HMMER and RepeatMasker programs.

# However, it can deal with TE candidates from other tools as well, since their final results (CANDIDATEs) are in the format:
# CANDIDATE_1 - FROM: pos_INIT - TO: pos_FINAL - LENGTH: length - SENSE: Direct/Reverse - CLASSIFICATION: superfam_A
# PREDIC---FROM--pos_from---TO--anyPos---LENGTH--length_1---METRICS--metr---SENSE--D/R---SUPERFAM--spfam_A
# PREDIC---FROM--anyPos---TO--anyPos---LENGTH--length_2---METRICS--metr---SENSE--D/R---SUPERFAM--spfam_A
# ...
# PREDIC---FROM--anyPos---TO--anyPos---LENGTH--length_N---METRICS--metr---SENSE--D/R---SUPERFAM--spfam_A
# (one blank line) (ATTENTION HERE)*********
# CANDIDATE_2 - FROM: pos_INIT - TO: pos_FINAL - LENGTH: length - SENSE: Direct/Reverse - CLASSIFICATION: superfam_B
# PREDIC---FROM--pos_from---TO--anyPos---LENGTH--length_1---METRICS--metr---SENSE--D/R---SUPERFAM--superfam_B
# ...
# ...
# "CANDIDATE_1", "CANDIDATE_2", etc. must be sorted by their "FROM" positions ("pos_INIT") inside the query sequence and the value of "METRICS" ("metr") would be used when filtering the predictions before producing the final candidates for each tool.

# The predictions ("PREDIC") that make up each "CANDIDATE" must be sorted by their "FROM" positions; then, the "FROM" position ("pos_from") of the first prediction (i.e., the lowest value among the "FROM" positions of all "PREDIC") of a CANDIDATE will be the "FROM" position ("pos_INIT") of that CANDIDATE. For the "TO" position ("pos_FINAL") of a CANDIDATE, select the highest value among the "TO" position ("anyPos") of all "PREDIC".

# WARNING: two things must be done for CombTEs works fine with a new tool:
# i) change the following parameters in the "ParamsGeneral.pm" file: "@otherTools", "%outFileNamesOthers", "%filterOtherTools", and "%goodMetrOtherTools" (see explanation in "ParamsGeneral.pm");
# ii) include a new "test condition" here - search for "OBS_2" below.

# For example, in the case of RPS-Blast program, the “finalCandidsRpsBlast.pl” script (also provided here) can be used to generate the final candidates in such a format.

# This script launches the "finalCandidsHmmerRM.pl" script, which generates the final candidates of HMMER and RepeatMasker, for all considered superfamilies (see more about the input files for "finalCandidsHmmerRM.pl" in such a script).

# If RPS-Blast is used in the searches for LTR-RT, delete the comments in front of the lines .... (search for "OBS_1" below).

# If a different tool is used to search for LTR-RT, a new "test condition" must be inserted at ... (search for "OBS_2 (a AND b)" below).

# This script uses 11 parameters, which can be changed in "ParamsGeneral.pm" package:


###########################################################################################

use strict;
use warnings;

use Cwd qw(getcwd);
use lib getcwd(); 

use ParamsGeneral qw(@superfamilies @combClassif @tools @outFileNames @otherTools %outFileNamesOthers %filterTools %filterOtherTools %goodMetrTools %goodMetrOtherTools $distBetweenCands);


my @classSpfam  = (@superfamilies, "inconclusive");

my $usedRpsblast = "no";
foreach my $tool (@otherTools) {
	push (@tools, $tool);
	push (@outFileNames, $outFileNamesOthers{$tool});
	$filterTools{$tool} = $filterOtherTools{$tool};
	$goodMetrTools{$tool} = $goodMetrOtherTools{$tool};
	if ($tool eq "RpsBlast") { $usedRpsblast = "yes"; };
}


###########################################################################################

	my $finalFile = "FinalClassification.txt";
	open (CLASS, ">$finalFile") or die "Can't open $finalFile";

	my $qttOtherTools = scalar(@otherTools);
	print CLASS "Final classification for \"HMMER\", \"RepeatMasker\"";
	for (my $i = 0; $i < $qttOtherTools; $i++) { print CLASS ", \"$otherTools[$i]\""; }
	print CLASS ".\n";

	print CLASS "Used filters: **HMMER = $filterTools{'HMMER'}**, **RepeatMasker = $filterTools{'RepeatMasker'}**";
	for (my $i = 0; $i < $qttOtherTools; $i++) {
		print CLASS ", **$otherTools[$i] = $filterTools{$otherTools[$i]}**";
	}
	print CLASS ".\n";
	print CLASS "Maximum distance from a candidate to the next one to consider them inside the same FINAL GROUP: $distBetweenCands.\n\n";

	my $finalFileTab = "FinalClassificationCSV.txt";
	open (CLASSTAB, ">$finalFileTab") or die "Can't open $finalFileTab";
	if ($usedRpsblast eq "no") {
		print CLASSTAB "Group,From,To,Length,Sense,TE classification,Identified by,Tools with good metrics ($goodMetrTools{'HMMER'} / $goodMetrTools{'RepeatMasker'}";
	}
	else {
		print CLASSTAB "Group,From,To,Length,Sense,TE classification,Identified by,Conserved Domains,Tools with good metrics ($goodMetrTools{'HMMER'} / $goodMetrTools{'RepeatMasker'}";
	}
	for (my $i = 0; $i < $qttOtherTools; $i++) { print CLASSTAB " / $goodMetrTools{$otherTools[$i]}"; }
	print CLASSTAB ")\n";

	my %filehandles = ();
	foreach my $fileName (@outFileNames) {
		my $fileFinders = "FinalClassification_$fileName.txt";
		open (my $fh, ">", "$fileFinders") or die "Can't open $fileFinders";
		$filehandles{$fileName} = $fh;
		print $fh "This file shows the groups produced by $fileName tool(s)\n\n";
	}


## insertion of ALL candidates (of ALL superfamilies) from ALL tools in "@allCandid" for FINAL classification
	my @allCandid = ();
	foreach my $tool (@tools) {
		my @candidSpfam = ();
		if ( ($tool eq "HMMER") or ($tool eq "RepeatMasker") )	{ @candidSpfam = `perl finalCandidsHmmerRM.pl $tool`; }

# WHEN using RPS-Blast: CombTEs assumes that "extractRPSB.pl" was used to extract and format the predictions from RPS-Blast and, then, its output file, with the conserved domains, was named "ConservedDomains_NoRedund.pred", to be used by "finalCandidsRpsBlast.pl":
		elsif ($tool eq "RpsBlast") { @candidSpfam = `perl finalCandidsRpsBlast.pl`; }
# OBS_2 (a):
# WHEN a different tool ("otherTOOL") is used to search for LTR-RT, the user must provide the correspondent script to process the initial candidates of this new tool; the name of the new script will replace "scriptOtherTOOL" (if the new script is written in Perl) in the next line (and then, DELETE the comment ("#") in front of the next line):
#		elsif ($tool eq "otherTOOL") { @candidSpfam = `perl scriptOtherTOOL`; }

		my ($from, $lineCand);
		foreach my $line (@candidSpfam) {
			if ($line =~ /FROM: (\d+) - TO/) {
				$from = $1;
				chomp $line;
				$lineCand = "\t$line\t\t$tool\n";
			}
			elsif ($line =~ /PREDIC|DOMAIN/) { $lineCand .= "\t\t$line"; }
			elsif ($line =~ /###/) {
				chomp $lineCand;
				push (@allCandid, {line => $lineCand, from => $from});
			}
		}
	} # FOREACH my $tool (@tools)

##    sorting the candidates
	my @auxSorted = sort{ $a->{from} <=> $b->{from} } @allCandid;
	my @sortedAllCandid = ();
	foreach my $lineAux (@auxSorted) { push (@sortedAllCandid, $lineAux->{line}); }

########### producing the groups of the final classification
	my %classif = ();
	foreach my $clsf (@classSpfam) { $classif{$clsf} = 0; }

	my %toolFinder = ();
	foreach my $tool (@tools) { $toolFinder{$tool} = 0; }

	my %typeClassFinal = ();
	foreach my $fileName (@outFileNames) {
		foreach my $comb (@combClassif) { $typeClassFinal{$fileName.$comb} = 0; }
	}

	my $qttCandid = scalar(@sortedAllCandid);

	my $group      = 0;
	my @arrayGroup = ();
	my $fromGroup  = 0;
	my $toGroup    = 0;
	my $newGroup   = "yes";
	my $i = 0;
	my ($from1, $to1, $sense1, $senseGroup, $from2, $sense2);
	while ($i < $qttCandid) {
		push (@arrayGroup, $sortedAllCandid[$i]);

# from the input files:
#	CANDIDATE_num - FROM: ff - TO: tt - LENGTH: ll - SENSE: d/r - CLASSIFICATION: Sf \n 
#	PREDIC---FROM--ff--- ...	(OR)
#	DOMAIN--dom---FROM-- ...
		if ($sortedAllCandid[$i] =~ /FROM: (\d+) - TO: (\d+) - LENGTH: \d+ - SENSE: (\w+) - CLASSIFICATION: (\w+)\t\t(\w+)\n\t\t(PREDIC|DOMAIN)/) {
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

		if ( (($from2-$toGroup) > $distBetweenCands) or ($sense2 ne $senseGroup) or (($i+1) == $qttCandid) ) {
# if (distance between next candidate and current GROUP is higher than "$distBetweenCands") OR (they are in different senses) OR (current candidate is the last one), then a new GROUP must be formed: write the current GROUP in the output files.
			my $finalClass = "";
			my $finalTools = "";
			$group++;
			my $lenGroup = $toGroup - $fromGroup + 1;

			foreach my $clsf (@superfamilies) {
				if ($classif{$clsf} == 1) {
					if ($finalClass eq "")  { $finalClass = $clsf; }
					else			{ $finalClass .= "/$clsf"; }
				}
			}
			if ($classif{"inconclusive"} == 1) {
				if ($finalClass eq "") { $finalClass = "inconclusive"; }
			}
			foreach my $clsf (@classSpfam) { $classif{$clsf} = 0; }

			my $numTools = 0;
			foreach my $tool (@tools) {
				if ($toolFinder{$tool} == 1) {
					if ($finalTools eq "")  { $finalTools = $tool; }
					else			{ $finalTools .= " / $tool"; }
					$numTools++;
				}
				$toolFinder{$tool} = 0;
			}
			my $indx = $numTools - 1;
			my $file = $outFileNames[$indx];
			my $FILEtoWRITE = $filehandles{$file};

			my $descrGroup = "GROUP $group - FROM: $fromGroup - TO: $toGroup - LENGTH: $lenGroup - SENSE: $senseGroup - CLASS: $finalClass - TOOL(S): $finalTools\n";
			print CLASS $descrGroup;
			print $FILEtoWRITE $descrGroup;


######## for the CSV file
			my $descrGroupTab = "$group,$fromGroup,$toGroup,$lenGroup,$senseGroup,$finalClass,$finalTools,";

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
				elsif ($arraySplit[$jj] =~ /SWSCORE--(\d+)---SENSE/) {
					if ($1 >= $goodMetrTools{'RepeatMasker'}) { $toolsGoodMetr{"RepeatMasker"} = 1;}
				}
# OBS_1: IF RPS-BLAST IS USED IN THE SEARCHES FOR LTR-RT, DELETE THE COMMENTS ("#") IN FRONT OF THE FOLLOWING 5 LINES:
#				elsif ($arraySplit[$jj] =~ /DOMAIN--(.*)---FROM--.*---EVALUE--(.*)---SCORE/) {
#					if ($foundDomains eq "") { $foundDomains = $1; }
#					else			 { $foundDomains .= " / $1"; }
#					if ($2 <= $goodMetrTools{'RpsBlast'}) { $toolsGoodMetr{"RpsBlast"} = 1; }
#				} # END OF "ELSIF" SPECIFIC FOR RPS-Blast.

# OBS_2 (b): IF A DIFFERENT TOOL IS CONSIDERED IN THE SEARCHES FOR LTR-RT, A NEW TEST CONDITION ("ELSIF") MUST BE INSERTED HERE; this new condition should filter the candidates based on the metrics used by the new tool.

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

			if ($usedRpsblast eq "no") {
				$descrGroupTab .= "$toolsWithGoodMetrics\n";
			}
			else {
				$descrGroupTab .= "$foundDomains,$toolsWithGoodMetrics\n";
			}
			print CLASSTAB $descrGroupTab;

			@arrayGroup = ();
			$fromGroup  = 0;
			$toGroup    = 0;
			$newGroup   = "yes";

			$typeClassFinal{$file.$finalClass}++;

		} # IF ( (($from2-$toGroup) > $distBetweenCands) or ...

		$i++;
	} # WHILE ($i < $qttCandid)
########### END of producing the groups of the final classification


## extra output info:
	print CLASS "#####################################################\n\n";

	print CLASS "SUMMARY OF THE GENERATED GROUPS:\n\n";

	foreach my $fileName (@outFileNames) {
		my $qttGroup = 0;
		foreach my $comb (@combClassif) { $qttGroup += $typeClassFinal{$fileName.$comb}; }
		print CLASS "Groups generated by $fileName tool(s): $qttGroup\n";
		print CLASS "   Final classification:\n";
		foreach my $comb (@combClassif) { print CLASS "      - $comb: $typeClassFinal{$fileName.$comb}\n"; }
		print CLASS "\n";
	}

## closing all output files
	foreach my $fileName (@outFileNames) { close ($filehandles{$fileName}); }
	close (CLASS);
	close (CLASSTAB);

