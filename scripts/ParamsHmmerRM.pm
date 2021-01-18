# Developed by Carlos Fischer - 20 Oct 2020

# Parameters used in "extractHmmerRM.pl" and "finalCandidsHmmerRM.pl" scripts.

###########################################################################################

package ParamsHmmerRM;
# ParamsHmmerRM.pm

use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw($varExtr $minOverlPred $maxOverlCand $distPredsHMMER $distPredsRepeatMasker $minLenPred $minSwscore $pattLtrs $includeLTRs);




####### used in "extractHmmerRM.pl":
our $varExtr = 50; # to allow small variation between respective extremities for cases of overlapping between 2 predictions (default= 50).
our $minOverlPred = 100; # minimum value to consider that there is overlapping between two predictions (default = 100).
#   $usingSubseqs = "yes"(default) - CHANGE IT IN "ParamsGeneral.pm".




####### used in "finalCandidsHmmerRM.pl":
# our $varExtr = 50; # (see above).
# our $minOverlPred = 100; # (see above).
our $maxOverlCand   = 100; # maximum value of overlap between 2 candidates to consider them SEPARATE final candidates (default = 100).

#   $filterHMMER    = 1.0e-05(default) - to filter HMMER's predictions based on e-values - CHANGE IT IN "ParamsGeneral.pm".
our $distPredsHMMER = 300;  # maximum distance between 2 predictions to consider them inside the SAME candidate (default = 300).

#   $filterRepeatMasker    = 225(default) - to filter RepeatMasker's predictions based on SWscores - CHANGE IT IN "ParamsGeneral.pm".
our $distPredsRepeatMasker = 500; # maximum distance between 2 predictions to consider them inside the SAME candidate (default = 500).
our $minLenPred            = 100; # minimum length for a prediction to be included in the analyses (default = 100).
our $minSwscore            = 500; # minimum score to consider a "short" prediction of RepeatMasker in the analyses: if (length < $minLenPred), its SWscore must be >= $minSwscore (default = 500).

our $pattLtrs    = ("-LTR|_LTR"); # the patterns used in Repbase to describe a sequence as a specific LTR one.
our $includeLTRs = "no";          # to consider ("yes") or not ("no") LTR predictions in the analyses (default = "no").


