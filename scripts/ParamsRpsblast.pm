# Developed by Carlos Fischer - 20 Oct 2020

# Parameters used in "extractRPSB.pl" and "finalCandidsRpsBlast.pl" scripts.

###########################################################################################

package ParamsRpsblast;
# ParamsRpsblast.pm

use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw($varExtr $minOverlPred %domains $distPredsRPSB);



####### used in "extractRPSB.pl":
our $varExtr = 50; # to allow small variation between respective extremities for cases of overlapping between 2 predictions (default= 50).
our $minOverlPred = 100; # minimum value to consider that there is overlapping between two predictions (default = 100).
#   $usingSubseqs = "yes"(default) - CHANGE IT IN "ParamsGeneral.pm".

####### Conserved domains considered in the analyses (after running RPS-Blast):
our %domains = ();
$domains{'Bel'}   = 'DUF1758|DUF1759|Peptidase_A17|RT_pepA17';
$domains{'Copia'} = 'Retrotran_gag|gag_pre-integrs|RVT_2|RNase_HI_RT_Ty1';
$domains{'Gypsy'} = 'Retrotrans_gag|gag-asp_proteas|retropepsin_like|RP_Saci_like|RVP_2|RT_LTR|RVT_3|RNase_HI_RT_Ty3|RNase_HI_like';



####### used in "finalCandidsRpsBlast.pl":
#   $filterRpsBlast = 1.0e-05 (default) - to filter RPS-Blast's predictions based on e-values - CHANGE IT IN "ParamsGeneral.pm".
our $distPredsRPSB  = 1300; # maximum distance between 2 predictions to consider them inside the SAME candidate (default = 1300);

