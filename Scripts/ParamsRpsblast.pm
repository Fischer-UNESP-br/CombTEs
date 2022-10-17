# Developed by Carlos Fischer - 20 Oct 2020
# Updated: 02 Aug 2022

# Parameters used in "extractRPSB.pl" and "finalCandidsRpsBlast.pl" scripts.

###########################################################################################

package ParamsRpsblast;
# ParamsRpsblast.pm

use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw($maxPercOutOverlap $maxRegionOutOverlap %domains $distPredsRPSB);


####### used in "extractRPSB.pl":
our $maxPercOutOverlap = 0.20; # percentage used to calculate the maximum size allowed for the region outside the overlapping region, to consider the overlap as (possibly) acceptable (default = 0.20 == 20% - a lower value would increase the number of predictions to be considered in the following analyses).
our $maxRegionOutOverlap = 100; # (predefined) maximum length allowed for region(s) outside the overlapping region between 2 predictions (default = 100 - a lower value would increase the number of predictions to be considered in the following analyses).


####### Conserved domains considered in the analyses (after running RPS-Blast):
our %domains = ();
$domains{'Bel'}   = 'DUF1758|DUF1759|Peptidase_A17|RT_pepA17';
$domains{'Copia'} = 'Retrotran_gag|gag_pre-integrs|RVT_2|RNase_HI_RT_Ty1';
$domains{'Gypsy'} = 'Retrotrans_gag|gag-asp_proteas|retropepsin_like|RP_Saci_like|RVP_2|RT_LTR|RVT_3|RNase_HI_RT_Ty3|RNase_HI_like';
# and also:   rve (Integrase), for Bel, Copia, and Gypsy.


####### used in "finalCandidsRpsBlast.pl":
#   $filterRpsBlast = 1.0e-05 (default) - to filter RPS-Blast's predictions based on e-values - CHANGE IT IN "ParamsGeneral.pm".
our $distPredsRPSB  = 300; # maximum distance between 2 predictions to consider them inside the SAME candidate (default = 300);

