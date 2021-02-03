# Developed by Carlos Fischer - 20 Oct 2020

# Parameters used in "extractHmmerRM.pl" and "finalCandidsHmmerRM.pl" scripts.

###########################################################################################

package ParamsHmmerRM;
# ParamsHmmerRM.pm

use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw($varExtr $maxPercOutOverlap $maxRegionOutOverlap $maxPercCandidOutOverlap $distPredsHMMER $distPredsRepeatMasker $minLenPred $minLenPredHMMER $minLenPredRepeatMasker $pattLtrs $includeLTRs);


####### used in "extractHmmerRM.pl":
our $varExtr = 30; # to allow small variation between respective extremities for cases of overlap between 2 predictions (default= 30).
our $minLenPred = 20; # minimum length for any prediction to be included in the analyses (default = 20).
our $maxPercOutOverlap = 0.20; # percentage used to calculate the maximum size allowed for the region outside the overlapping region, to consider the overlap as (possibly) acceptable (default = 0.20 == 20% - a lower value would increase the number of predictions to be considered in the following analyses).
our $maxRegionOutOverlap = 100; # (predefined) maximum length allowed for region(s) outside the overlapping region between 2 predictions (default = 100 - a lower value would increase the number of predictions to be considered in the following analyses).
#   $usingSubseqs = "yes"(default) - CHANGE IT IN "ParamsGeneral.pm".



####### used in "finalCandidsHmmerRM.pl":
our $maxPercCandidOutOverlap = 0.50; # maximum percentage of the smallest CANDIDATE allowed to be outside the overlapping region, to consider the overlap as allowed (default = 0.50 == 50% - a lower value would increase the number of CANDIDATES for a specific tool).

#   $filterHMMER     = 1.0e-03(default) - to filter HMMER's predictions based on e-values - CHANGE IT IN "ParamsGeneral.pm".
our $distPredsHMMER  = 300; # maximum distance between 2 predictions to consider them inside the SAME candidate (default = 300).
our $minLenPredHMMER = 20;  # minimum length for a HMMER's prediction to be included in the analyses (default = 20).

#   $filterRepeatMasker    = 225(default) - to filter RepeatMasker's predictions based on SWscores - CHANGE IT IN "ParamsGeneral.pm".
our $distPredsRepeatMasker = 300; # maximum distance between 2 predictions to consider them inside the SAME candidate (default = 300).
our $minLenPredRepeatMasker = 20; # minimum length for a RepeatMasker's prediction to be included in the analyses (default = 20).
our $pattLtrs    = ("-LTR|_LTR"); # the patterns used in Repbase to describe a sequence as a specific LTR one.
our $includeLTRs = "no";          # to consider ("yes") or not ("no") LTR predictions in the analyses (default = "no").


