# Developed by Carlos Fischer - 20 Oct 2020
# Updated: 02 Aug 2022

# General parameters used in all scripts.

###########################################################################################

package ParamsGeneral;
# ParamsGeneral.pm

use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw(@TEtypes @tools @outFileNames %filterTools %goodMetrTools $distBetweenCands);


# ------------------------------------
# TE TYPES - to specify the "TE types" considered in the analyses:
our @TEtypes = ("Bel", "Copia", "Gypsy"); # in this case, the "superfamily level" would be the "TE type" considered in the analyses
#our @TEtypes = ("Ale", "Alesia", "Ivana"); # here, the "TE type" would be "lineages"

# ------------------------------------
# TOOL - to specify the programs used in the searches for TEs:
our @tools = ("HMMER", "RepeatMasker", "RpsBlast"); # in this case, these 3 programs would have been used to search for TEs
#our @tools = ("otherTool1", "otherTool2", "otherTool3", "otherTool4"); # when using 4 programs

our @outFileNames  = ("ONE", "TWO", "THREE", "FOUR", "FIVE"); # if you use more than 5 programs, include new names in "@outFileNames"


# ------------------------------------
# FILTER VALUES - include/change HERE the filter values for a specific tool:
our %filterTools = ();
$filterTools{'HMMER'} = 1.0e-05;    # to filter HMMER's predictions based on e-values (default = 1.0e-05)
$filterTools{'RepeatMasker'} = 230; # to filter RepeatMasker's predictions based on SWscores (default = 225)
$filterTools{'RpsBlast'} = 1.0e-05; # to filter RpsBlast's predictions based on e-values (default = 1.0e-05)
#$filterTools{'otherTool'} = xxxx;  # to filter otherTool's predictions based on the used metrics of that tool


# ------------------------------------
# "GOOD METRICS" - to set the "good metrics":
# (the term "good metrics" for a tool is used here to define the value of the metrics used by that tool for which the user would consider a prediction as (possibly) a correct one; for example, for "$goodMetrTools{'HMMER'} = 1.0e-20", the sequences predicted by HMMER with E-values equal or lower than "1.0e-20" would be considered correct ones by that user):
our %goodMetrTools = ();
$goodMetrTools{'HMMER'}        = 1.0e-20;
$goodMetrTools{'RepeatMasker'} = 400;
$goodMetrTools{'RpsBlast'}     = 1.0e-20;
#$goodMetrTools{'otherTool'}    = xxxx; # the "good metrics" for the "otherTool" tool


# ------------------------------------
# DISTANCE BETWEEN CANDIDATES:
our $distBetweenCands = 500; # maximum distance (default = 500 nt) from a candidate to the next one to consider them inside a SAME FINAL GROUP


