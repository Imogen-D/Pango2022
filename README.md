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

  Calls nQuire_ploidy.sh run_vanquish.sh (which calls vanquish_ploidy_plot.R)
  
03.ANGSD.sh - estimate GL and getting beagle, running ngsADMIX, producing likelihood file, plotting + geoplotting - IBD on this script obselete

  Calls ANGSD.sh, ngsadmix.sh cuttingbeagle.sh ngsadmix_newbeagle.sh  ngsadmix_plots.R ngsadmix_geoplots.R pcangsd_pango.sh angsd.he.indv.sh IBD_plotting.R

03b.IBD_script.sh - does IBD 

  Calls cuttingbeagle.sh; To run IBD, 1) remake list of samples to exclude 2) run beagle filtering 3) run pcagnsd with all 4) concat X + Auto 5) run IBD with all three datasets

04.locator.sh - runs locator, makes vcfs (indvidual and all lineages), plotting

  Calls locator.sh iterateseedlocator.sh ALvcfmaking.sh indvlocvcfs.sh plot_locator_sbatch.sh zipping.sh
  
05.iterativeangsd.sh - to see if iterating angsd over minimal gmin calls changes he (it doesn't)

  Calls angsd.he.indv.sh
  
06.privatealleleresampling.R - resampling PAs

08.assignmentmethods.sh - run BONE or rubias or assignPOP (tested but not in manuscript) plus make initial file from vcf

  Calls run.assignment.sh
