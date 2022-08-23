# Pango2022
Scripts for tracing the trade with ddRADseq SNP data.

Scripts start on JSalmona Account for initial demultiplexing

Scripts and workflow on IDumville account from samfiles
As follows;

01.alignmentstats_coverage.sh - getting coverage, extract statistics
  Calls coverage.sh loci.sh (x)cov.plot.all.R loci.plot.R (via intemediate loci.sh, but no longer availble?)
  
02.stacks_populations.sh - running gstacks (on seperate miseq/other datasets), populations (some seperation for lowreads / different localities, getting private alleles (PAs), unused script for error rate
  Calls run.gstacks.sh run.population.sh (change Rs, minmac 2-3, hap exports, output file) run.indv.popn.sh run.indvseedpopn.sh
  
02b.ploidy.sh - getting ploidy from bams, check for contamination
  calls nQuire_ploidy.sh run_vanquish.sh (which calls vanquish_ploidy_plot.R)
  
03.ANGSD.sh - estimate GL and getting beagle, running ngsADMIX, producing likelihood file, plotting + geoplotting
  Calls ANGSD.sh, ngsadmix.sh cuttingbeagle.sh ngsadmix_newbeagle.sh  ngsadmix_plots.R ngsadmix_geoplots.R pcangsd_pango.sh angsd.he.indv.sh IBD_plotting.R

03b.IBD_script.sh
  Calls cuttingbeagle.sh
