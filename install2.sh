
# First install and decompress files

# This line decompresses a big file needed to work

echo Decompressing CDS_genomic_coordinates_full_compresed.txt.tar.gz
tar -xzvf data/build_tables/CDS_genomic_coordinates_full_compresed.txt.tar.gz -C data/build_tables/
rm data/build_tables/CDS_genomic_coordinates_full_compresed.txt.tar.gz

echo Decompressing snpEff v3_6
mkdir SOFTWARE/snpEff
unzip SOFTWARE/snpEff_v3_6_core.zip -d ./SOFTWARE/
rm SOFTWARE/snpEff_v3_6_core.zip

# The genome files are obtained from ELSEWHERE (#ISSUE 1.1) ALL OF THESE FILES ARE OBTAINED FROM LOCAL THIS HAS TO CHANGE
# Cahce files from ENSEMBL are obtained from their website include snpeff ones in the package somehow
echo obtaining the human genome version GRCh37.75

# Obtain the genome from SnpEff database

wget http://downloads.sourceforge.net/project/snpeff/databases/v3_6/snpEff_v3_6_GRCh37.75.zip

mv  snpEff_v3_6_GRCh37.75.zip ./SOFTWARE/snpEff/
unzip SOFTWARE/snpEff/snpEff_v3_6_GRCh37.75.zip -d ./SOFTWARE/snpEff/
rm SOFTWARE/snpEff/snpEff_v3_6_GRCh37.75.zip


# VARIANT EFFECT PREDICTOR

# Decompress the VEP.zip file

#unzip SOFTWARE/ensembl-tools-release-81.zip -d ./SOFTWARE/
#echo INSTALLING VEP
#perl SOFTWARE/ensembl-tools-release-81/scripts/variant_effect_predictor/INSTALL.pl

# Answer Yes to install the API, and yes to install cache files the one needed is 
# 45 : homo_sapiens_vep_81_GRCh37.tar.gz


# Creating directories needed to run nutvar2

mkdir data/intermediate/

