# Developed by Carlos Fischer - 20 Oct 2020

# General parameters used in all scripts.

###########################################################################################

package ParamsGeneral;
# ParamsGeneral.pm

use strict;
use warnings;

use Exporter qw(import);

# When using also other tools, use the following "@EXPORT_OK":
our @EXPORT_OK = qw(@superfamilies @combClassif $usingSubseqs $overlap $limit @tools @outFileNames @otherTools %outFileNamesOthers %filterTools %filterOtherTools %goodMetrTools %goodMetrOtherTools $distBetweenCands);


our @superfamilies = ("Bel", "Copia", "Gypsy");
our @combClassif = (@superfamilies, "inconclusive", "Bel/Copia", "Bel/Gypsy", "Copia/Gypsy", "Bel/Copia/Gypsy"); # if a new superfamily is included in the analyses, insert ALL possible combinations with the other superfamilies in this array.



our $usingSubseqs = "yes"; # when using subsequences of a long sequence (OR: use $usingSubseqs = "no" for a long sequence and also for general small sequences)
my  $lengthSubseq = 10000; # length of the subsequences.
our $overlap      = 1000;  # length of the overlap region.
our $limit        = $lengthSubseq - $overlap;



our @tools = ("HMMER", "RepeatMasker");
our @outFileNames  = ("ONE", "TWO"); # when using only HMMER and RepeatMasker.

# IF "RPS-Blast" is to be used in the analyses, search for "RPS-Blast" below and do the requested action.

our @otherTools = (); # if a new tool is to be used in the analyses, insert its name in this array - DO NOT COMMENT THIS LINE.
#@otherTools = ("RpsBlast"); # if "RPS-Blast" is used: delete the comment ("#") in front of this line.
#@otherTools = ("OTHERTOOL_1"); # and so on.
our %outFileNamesOthers = (); # when using OTHER tools.
#$outFileNamesOthers{'RpsBlast'} = ("THREE"); # if "RPS-Blast" is used: delete the comment ("#") in front of this line.
#$outFileNamesOthers{'OTHERTOOL_1'} = ("FOUR"); and so on.



# change HERE the filter values for the specific tool:
our %filterTools = ();
$filterTools{'HMMER'} = 1.0e-05;    # to filter HMMER's predictions based on e-values (default = 1.0e-05).
$filterTools{'RepeatMasker'} = 225; # to filter RepeatMasker's predictions based on SWscores (default = 225).
our %filterOtherTools = (); # if a new tool will be used in the analyses, insert its filter value in this hash - DO NOT COMMENT THIS LINE.
#$filterOtherTools{'RpsBlast'} = 1.0e-05; # if "RPS-Blast" is used: delete the comment ("#") in front of this line.
#$filterOtherTools{'otherTOOL_1'} = "filter_value"; # filter considering the used metrics of other tool.



# to set the "good metrics" for HMMER and RepeatMasker:
# (the term "good metrics" for a tool is used here to define the value of the metrics used by that tool for which the user would consider a prediction as a correct one; for example, for "$goodMetrTools{'HMMER'} = 1.0e-50", the sequences predicted by HMMER with E-values equal or lower than "1.0e-50" would be considered correct ones by that user):
our %goodMetrTools = ();
$goodMetrTools{'HMMER'}        = 1.0e-50;
$goodMetrTools{'RepeatMasker'} = 500;
our %goodMetrOtherTools = (); # if a new tool will be used in the analyses, insert its "good metrics" value in this hash - DO NOT COMMENT THIS LINE.
#$goodMetrOtherTools{'RpsBlast'} = 1.0e-20; # if "RPS-Blast" is used: delete the comment ("#") in front of this line.
#$goodMetrOtherTools{'otherTOOL_1'} = "filterValue";


our $distBetweenCands = 500; # maximum distance from a candidate to the next one to consider them inside a SAME FINAL GROUP.


