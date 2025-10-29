MEDOC_public
Public release of the MEDOC algorithm for the prediction of protein protonation states MEDOC
Please use latest version ! 1.4.4 as of 27/10/2025

Three version are provided depending on the precision used to store float : single (s),long (l) and longdouble (ld)
Presumably single is faster but higher precision allow to tackle longer sequences, although after the considerable speedups from V1.4.4+, this may no longer matter.
A log of the changes is available in a pdf on github.


To use MEDOC, open a terminal and type : python3 MEDOC_VX.py [command line arguments]
The following command line arguments can be specified to change the default behavior.


Major parameters:
Handle      Dtype   Default     Description
    -s      str     seq.fasta   Name of the sequence file which must be a single line of uninterrupted natural amino acids. Default
    -t      float   298.        Temperature in Kelvin
    -nn     int     2           Number of neighbors that influence the ionization free energy of a residue. Default is 2.

Cosmetic parameters:
Handle      Dtype   Default     Description
    -ds     int     0           Detailed suffix for file names. If 1, all file names will have additional ditails such as t, nn, pt, p.
    -dl     int     2           Detail level for figure and message printing: 0 is no figure, 1 adds the net charge vs pH plot and the charge density vs pH (if site specific is on) , 2 adds the mesostate plot and the plot of all sites protonation curve for all residue types at once (if site specific is on), 3 adds the separate site specific curves for each amino acid type (if site specific is on), 4 adds a figure for each protonation sites with an analysis of the cooperativity and transition asymmetry (if site-specific is on).
    -pr     float   0.01        Plots resolution in the pH space : Default is 0.01. Greatly affects figure printing speed, but has no effect on the speed of the algorithm itself.
    -pt     str     I           Prediction type : affects the free energy used for each charge context. Can be Implicit (I) or Unshifted (U). Due to the size and non-exhaustive nature of the explicit context database, it is not available to the public.
    -ss     int     0           Site specific : 0 is the global prediction (default), one is the site specific prediction
    -be     float   -709*R*T    This is the value for the reference energy. Base value is -709*R*T. This can be changed to get around the partition function "breaking" due to exceeding max floating point value. As of 1.4.4+ this should no longer be an issue.
    -p      str     d           Pruning. Can be set to d (during), very much recomended, of a (after). If after is used, the selection of states will be more accurate, but you will be limited to sequnces of less than 15 ionizable AAs or so.
    -md     float   0           Max_difference for the kept states. Only valid if -pruning is used. This max diff is the maximum difference between the overall states ensemble free energy and the discarded states ensemble free energy. If <=0, all states are Kept, increasingly large value results in less states being kept, with a minimum of 1. Remember to change md as well.
    -pH     float   0 14        Pair of values representing the upper and lower limit of the pH range to be plotted. Does not affect calculation speed, only figure printing.
If the procedure does not get to the end and has an error message, this may be due to a system dependent floating point issues. Try a machine with higher precision. This was tested on a linux machine with a 64-bit  python distribution.

You will get the error message : RuntimeWarning: divide by zero encountered in log This is normal.


