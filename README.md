# Pango2022

Scripts for tracing the trade with ddRADseq SNP data.

This GitHub is a work in progress (I apologise for the inconsistent numbering). All scripts excuted by Imogen Dumville, written by Imogen Dumville with aid of Jordi Salmona. 

**Scripts start on JSalmona Account for initial demultiplexing + alignment**

syntenyanalysis.sh - finds which chromosome in reference is the X chromosome

  *Calls mummer.sh*

script.2.demultiplex_data.sh - demultiplexing

  _Calls process_radtags.pg.sh_
  
script.3.align_rad_data.sh - trimmming then aligns all the data to autosomes, X, mito and checks files

  _Calls trimmomatic_low_PE.sh index_genome.sh bwa_miseq2021.sh_
  
script.3b.chromosomes.sh -indexes and blasts the mtDNA and Xchromosome as well as seperating the X from autosomes

  _Calls index_genome.sh blastmito.sh blastxchr.sh sepx_a_stats.sh_



**Scripts and workflow on IDumville account from samfiles**

As follows;

01.alignmentstats_coverage.sh - getting coverage, extract statistics

 _Calls coverage.sh loci.sh (x)cov.plot.all.R loci.plot.R (via intemediate loci.sh, but no longer availble - I had an internet issue in my project so may have been lost here)_
  
02.stacks_populations.sh - running gstacks (on seperate miseq/other datasets), populations (some seperation for lowreads / different localities, getting private alleles (PAs), unused script for error rate

 _Calls run.gstacks.sh run.population.sh (change Rs, minmac 2-3, hap exports, output file) run.indv.popn.sh run.indvseedpopn.sh_
  
02b.ploidy.sh - getting ploidy from bams, check for contamination

 _Calls nQuire_ploidy.sh run_vanquish.sh (which calls vanquish_ploidy_plot.R)_
  
03.ANGSD.sh - estimate GL and getting beagle, running ngsADMIX, producing likelihood file, plotting + geoplotting - IBD on this script obselete

  _Calls ANGSD.sh, ngsadmix.sh cuttingbeagle.sh ngsadmix_newbeagle.sh  ngsadmix_plots.R ngsadmix_geoplots.R pcangsd_pango.sh angsd.he.indv.sh IBD_plotting.R_

03b.IBD_script.sh - does IBD 

  _Calls cuttingbeagle.sh; To run IBD, 1) remake list of samples to exclude 2) run beagle filtering 3) run pcagnsd with all 4) concat X + Auto 5) run IBD with all three datasets_

04.locator.sh - runs locator, makes vcfs (indvidual and all lineages), plotting

 _Calls locator.sh iterateseedlocator.sh ALvcfmaking.sh indvlocvcfs.sh plot_locator_sbatch.sh zipping.sh_
  
05.iterativeangsd.sh - to see if iterating angsd over minimal gmin calls changes he (it doesn't)

_Calls angsd.he.indv.sh_
  
script.iterate_p_alleles.sh - resampling PAs

 _Calls run.iterativepop.sh_

08.assignmentmethods.sh - run BONE or rubias or assignPOP (tested but not in manuscript) plus make initial file from vcf

 _Calls run.assignment.sh_
