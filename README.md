# nutvar2
Classifier of the potential pathogenicity of human genomic truncations

# INSTALLATION

# You need to have installed the variable effect predictor version 81.
# You can obtain it from:

# http://jul2015.archive.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer
# Follow instructions to download and install it
# Install cache number 45 (Homo sapiens GRC37)

# NUTVAR2 INSTALLATION
  
Download the .zip files with  git clone https://github.com/manueltar/NutVar2.git

unzip  nutvar2-master.zip

cd nutvar2-master
  # Install NutVar2
  
  ./install.sh
  

# RUNNING


# There are two different options to run NutVar2: using SnpEff as the only predictor for variant outcome, using VEP as the only predictor for variant outcome.

**Note: VEP is a very comprehensive and informative predictor of variant outcome, but its time of execution is way larger than that of SnpEff. If VEP is going to be used in large input vcf files, cut files in smaller subfiles prior to running NutVar2.


cd nutvar2-master

# 1- SnpEff

nutvar2-master$ ./NutVar2_snpEff.sh path-to-nutvar2-master/ user.vcf data/final

Test: nutvar2-master$ ./NutVar2_snpEff.sh ~/Downloads/nutvar2-master example.vcf data/final

# 2- VEP

nutvar2-master$ ./NutVar2_VEP.sh path-to-nutvar2-master/ user.vcf data/final path-to-VEP81/

Test:nutvar2-master$ ./NutVar2_VEP.sh  ~/Downloads/nutvar2-master example.vcf data/final path-to-VEP81/


# Issues

  # Installation ISSUE!! The genomes are huge where are we going to allocate them online for the user to download? Right now I retrieve them from my local disk. There is a genome for snpeff and a genome for VEP. Consider it.
  
  # Installation ISSUE!! Don't think its really necessary to set the path variables because the API is being installed directly with VEP
  
  # Running Issue 2 !!-> Allow some MPI option specially for VEP and also intra bash script, run script 25 while running the rest of scripts.
  
  # Running Issue 3: ExAC (22 Gb input vcf) has generated in the run 53 Gb of intermediate and final files. Consider erasing files whenever they are no longer needed downstream in the pipeline.
  
  # Update the VEP parser
  # Optional parse from the snpEff and VEP if the splice variant entails splice abrogation
  # Convert to GRCh38
  
# nutvar2_NEW
