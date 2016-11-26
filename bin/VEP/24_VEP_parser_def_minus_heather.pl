###############
###############

#!/usr/bin/perl

use strict;
use warnings;
use Time::localtime;


## Declare the array of accepted chromosomes. We won't accept variants in patches.

my @AcceptedChrArray=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");


## Create the hash for accepted chromosomes.

my %AcceptedChr=();
foreach my $AcceptedChrTmp(@AcceptedChrArray){$AcceptedChr{"$AcceptedChrTmp"}=1;}


## Create the hash with the effects we are going to accept. These are: missense_variant, stop_gained, frameshift_variant, stop_retained_variant, start_retained_variant,
## synonymous_variant,splice_acceptor_variant, splice_donor_variant, splice_region_variant, inframe_deletion, disruptive_inframe_deletion, exon_loss_variant,
## and start_lost.


my %AcceptedEffects=();
$AcceptedEffects{"splice_acceptor_variant"}=1;
$AcceptedEffects{"splice_donor_variant"}=1;
#$AcceptedEffects{"splice_region_variant"}=1;

#~ $AcceptedEffects{"exon_loss_variant"}=1;
#~ $AcceptedEffects{"transcript_ablation"}=1;

#~ $AcceptedEffects{"inframe_deletion"}=1;
#~ $AcceptedEffects{"disruptive_inframe_deletion"}=1;

$AcceptedEffects{"missense_variant"}=1; # missense_variant    (CURRENT_SVN) SO Accession:	SO:0001583 ANNOVAR:nonsynonymous SNV, VAAST:non_synonymous_codon, missense, missense codon, vep:NON_SYSNONYMOUS_CODING, SO:0001584, SO:0001783

$AcceptedEffects{"stop_gained"}=1; # término SO: stop_gained    (CURRENT_SVN) SO Accession:	SO:0001587       ANNOVAR:stopgain, nonsense, nonsense codon, vep:STOP_GAINED, stop gained, VAAST:stop_gained

$AcceptedEffects{"frameshift_variant"}=1; # término SO: frameshift_variant    (CURRENT_SVN) SO Accession:	SO:0001589  ANNOVAR:frameshift block substitution, frameshift variant, frameshift_, frameshift_coding, vep:FRAME_SHIFT, VAAST:frameshift_variant


#~ $AcceptedEffects{"stop_retained_variant"}=1; # término SO: stop_retained_variant  SO Accession:	SO:0001567  (CURRENT_SVN) vep:SYNONYMOUS_STOP, stop retained variant, VAAST:stop_retained
#~ $AcceptedEffects{"start_retained_variant"}=1; # término SO: start_retained_variant    (CURRENT_SVN) SO Accession:	SO:0002019 synonymous: vep??????????
$AcceptedEffects{"synonymous_variant"}=1; # término SO: synonymous_variant    SO Accesion: SO:0001819		ANNOVAR:synonymous SNV, silent mutation, silent substitution, silent_mutation, coding-synon, vep:SYNONYMOUS_CODING, synonymous codon, synonymous_coding, synonymous_codon, VAAST:synonymous_codon, SO:0001588, SO:0001588

#~ $AcceptedEffects{"start_lost"}=1;

## This is paralelized in 10 files so input and output files vary from 01 to 10, see the BASH SCRIPT to run it in a cluster.
## As the input file has been splited in 10 files and it is
## time consuming to copy the heather lines to each of them they have been copied in a separate file (heathers_VEP.txt) which is used to load the hashes.


my $list_2= $ARGV[0];
my $input2= $ARGV[1];
my $OUT_vep_parsed=$ARGV[2];
my $vep_NMD=$ARGV[3];

my %hash1=();
my %hash_BIAS=();

my $time='['. timestamp(). ']'."\n";
#~ print "Start reading INFILE:$time\n";

if(open (LIST_2, $list_2))
{
## Here we declare the hashes we are going to use to load the fields from the heather lines. 

my %HI0_vep=();
my %IH0_vep=();
my %HI1_vep=();
my %IH1_vep=();
my %HI2_vep=();
my @accepted_keys=();

	while (my $line = <LIST_2>)
	{
	chomp ($line);
		#~ print "------------------->$line\n";
		if($line=~/^#/)
		{
			if($line=~/^##INFO=<ID=CSQ.+ Format: (.+)">/)	
			{
				my @tmp_HI1 = split (/\|/,$1);
							
				for(my $i=0;$i<scalar(@tmp_HI1);$i++)
				{
					$HI1_vep{$tmp_HI1[$i]}=$i;
					$IH1_vep{$i}=$tmp_HI1[$i];
				}		
				my @dumkeys=keys (%HI1_vep);
				foreach my $dumkeys_tok(@dumkeys)
				{
					#~ print "H:$dumkeys_tok:H\n";
				}	
			}		
			elsif (index($line,'#CHROM')==0) 	
			{
				$line=~/^#(.+)/;
				my @tmp_HI0 = split (/\t/,$1);
									
				for(my $i=0;$i<scalar(@tmp_HI0);$i++)						
				{
					$HI0_vep{$tmp_HI0[$i]}=$i;
					$IH0_vep{$i}=$tmp_HI0[$i];
				}			
				my @dumkeys=keys (%HI0_vep);
				foreach my $dumkeys_tok(@dumkeys)
				{
				#~ print "H***:$dumkeys_tok:**H\n";
				}
			}	
		}
		elsif ($line !~ /^#/)
		{
							#~ print "Hello_worldI:$line\n";
							my @tmp = split (/\t/,$line);
							#~ print "El array es:@tmp\n";
							
							#~ print "***$HI0_vep{"CHROM"}\n";
							#~ print "***$tmp[$HI0_vep{"CHROM"}]\n";
					  					  
							my $CHROM="NaN";if (exists($HI0_vep{"CHROM"})){$CHROM=$tmp[$HI0_vep{"CHROM"}];}unless(defined($CHROM)){$CHROM="NaNein";}
							my $FILTER="NaN";if (exists($HI0_vep{"FILTER"})){$FILTER=$tmp[$HI0_vep{"FILTER"}];}unless(defined($FILTER)){$FILTER="NaNein";}
							my $INFO="NaN";if (exists($HI0_vep{"INFO"})){$INFO=$tmp[$HI0_vep{"INFO"}];}unless(defined($INFO)){$INFO="NaNein";}
						
							#~ print "****************$CHROM\t$FILTER\t$INFO\n";
							
							#~ exit;
						## First filter; not funny chromosomes, neither Mitochndrial. Not bad quality variants not variants without CSQ= field.
						
						if((exists($AcceptedChr{"$CHROM"})) and ($INFO =~ /CSQ\=/))
						
							{
								my $POS="NaN";if (exists($HI0_vep{"POS"})){$POS=$tmp[$HI0_vep{"POS"}];}unless(defined($POS)){$POS="NaNein";}
								my $ID="NaN";if (exists($HI0_vep{"ID"})){$ID=$tmp[$HI0_vep{"ID"}];}unless(defined($ID)){$ID="NaNein";}
								my $REF="NaN";if (exists($HI0_vep{"REF"})){$REF=$tmp[$HI0_vep{"REF"}];}unless(defined($REF)){$REF="NaNein";}
								my $ALT="NaN";if (exists($HI0_vep{"ALT"})){$ALT=$tmp[$HI0_vep{"ALT"}];}unless(defined($ALT)){$ALT="NaNein";}
								my $QUAL="NaN";if (exists($HI0_vep{"QUAL"})){$QUAL=$tmp[$HI0_vep{"QUAL"}];}unless(defined($QUAL)){$QUAL="NaNein";}
								my $FORMAT="NaN";if (exists($HI0_vep{"FORMAT"})){$FORMAT=$tmp[$HI0_vep{"FORMAT"}];}unless(defined($FORMAT)){$FORMAT="NaNein";}
								my $Sophia_genetics="NaN";if (exists($HI0_vep{"Sophia_genetics"})){$Sophia_genetics=$tmp[$HI0_vep{"Sophia_genetics"}];}unless(defined($Sophia_genetics)){$Sophia_genetics="NaNein";}
								
								#~ print "****************$CHROM\t$POS\t$ID\t$REF\t$QUAL\t$FILTER\n";
								#~ exit;
								
								my @INFO_derived = split ("\;",$INFO);
								
								## If possible we want to calculate the AFassessed.
								
								my $AC="NaN";
								my $AN="NaN";
								my $AF="NaN";
								my $AFassessed="NaN";
								
								foreach my $INFO_derived_tok (@INFO_derived)
								
									{											
										if($INFO_derived_tok=~/^AC\=/){$INFO_derived_tok=~s/^AC\=//;$AC=$INFO_derived_tok;}
										if($INFO_derived_tok=~/^AN\=/){$INFO_derived_tok=~s/^AN\=//;$AN=$INFO_derived_tok;}
										if($INFO_derived_tok=~/^AF\=/){$INFO_derived_tok=~s/^AF\=//;$AF=$INFO_derived_tok;}
									}
								
								if(($AC=~/[0-9]/)and($AN=~/[0-9]/)){$AFassessed=$AC/$AN;}
								
																
								## Here we start parsing the CSQ fields
								
								my @CSQfield = split ("CSQ\=",$INFO);
								
								
								my @CSQfield_derived =  split ("\,",$CSQfield[1]);
										
										foreach my $CSQfield_derived_tok (@CSQfield_derived)
										
										{	
											{
												my @CSQfield_derived_tok_tmp=split(/\|/,$CSQfield_derived_tok);
												
												## VEP has compound effects; i.e. stop_gained&NMD_variant&feature_truncation. We need to add these compound
												## effects to our hash of accepteed effects. Otherwise we will miss compound effects carrying the effects
												## we selected at the beggining.
												
												## First we select compound effects with our selected effects.
												
												my @compound_effects=grep (/missense_variant&.*/|/stop_gained&.*/|/frameshift_variant&.*/|/stop_retained_variant&.*/|/synonimous_variant&.*/|/start_retained_variant&.*/|/disruptive_inframe_deletion&.*/|/start_lost&.*/|/coding_sequence_variant&.*/|/splice_acceptor_variant&.*/|/splice_donor_variant&.*/|/exon_loss_variant&.*/|/transcript_ablation&.*/|/inframe_deletion&.*/,@CSQfield_derived_tok_tmp);

													## Then we add them to the hash.
													
													foreach my $compound_effects_tok(@compound_effects)
															{
																if (exists($HI2_vep{$compound_effects_tok})){$HI2_vep{$compound_effects_tok}++;} else{$HI2_vep{$compound_effects_tok}=1;}
																@accepted_keys=keys(%HI2_vep);
																foreach my $accepted_keys_tok (@accepted_keys)
																{
																	$AcceptedEffects{$accepted_keys_tok}=1;
																}														
															}
												
												
												## Define the Consequence field.
												
												my $Consequence="NaN";if (exists($HI1_vep{"Consequence"})){$Consequence=$CSQfield_derived_tok_tmp[$HI1_vep{"Consequence"}]} unless(defined($Consequence)){ print "ERROR in Consequence";}
												
												## If the consequence is among our accepted effects then we parse the fields it carries.
												
												if(exists($AcceptedEffects{$Consequence}))
												{	
													#print "HELLO_WORLD_II:$Consequence\n";
													#exit;
													
													my ($PUBMED,$Existing_variation,$CLIN_SIG,$AMR_MAF,$Amino_acids, $SYMBOL_SOURCE,$Protein_position,$Allele,$Codons,$CANONICAL, $DISTANCE,$ASN_MAF,$HIGH_INF_POS,$BIOTYPE,$CDS_position, $MOTIF_POS,$Feature_type,$AA_MAF,$Feature,$MOTIF_NAME, $CCDS,$STRAND,$EUR_MAF,$EXON,$EA_MAF, $PolyPhen,$MOTIF_SCORE_CHANGE,$cDNA_position,$ENSP,$Gene, $INTRON,$SYMBOL,$SIFT,$GMAF,$AFR_MAF, $DOMAINS)=("NaN")x36;

													if (exists($HI1_vep{"PUBMED"})){$PUBMED=$CSQfield_derived_tok_tmp[$HI1_vep{"PUBMED"}];unless (defined($PUBMED)){$PUBMED="NaN";}}else{ print "ERROR in PUBMED \n";}
													# NOT AVAILABLE W/O INTERNET CONNECTION my $HGVSc="NaN";if (exists($HI1_vep{"HGVSc"})){$HGVSc=$CSQfield_derived_tok_tmp[$HI1_vep{"HGVSc"}];if ($HGVSc!~/\w/){$HGVSc="NaN";}}else{ print "ERROR in HGVSc \n";}
													if (exists($HI1_vep{"Existing_variation"})){$Existing_variation=$CSQfield_derived_tok_tmp[$HI1_vep{"Existing_variation"}];unless (defined($Existing_variation)){$Existing_variation="NaN";}}else{ print "ERROR in Existing_variation \n";}
													if (exists($HI1_vep{"CLIN_SIG"})){$CLIN_SIG=$CSQfield_derived_tok_tmp[$HI1_vep{"CLIN_SIG"}];unless (defined($CLIN_SIG)){$CLIN_SIG="NaN";}} else{ print "ERROR in CLIN_SIG \n";}
													if (exists($HI1_vep{"AMR_MAF"})){$AMR_MAF=$CSQfield_derived_tok_tmp[$HI1_vep{"AMR_MAF"}];unless (defined($AMR_MAF)){$AMR_MAF="NaN";}}else{ print "ERROR in AMR_MAF\n";}
													if (exists($HI1_vep{"Amino_acids"})){$Amino_acids=$CSQfield_derived_tok_tmp[$HI1_vep{"Amino_acids"}];unless (defined($Amino_acids)){$Amino_acids="NaN";}} else{ print "ERROR in Amino_acids \n";}
													if (exists($HI1_vep{"SYMBOL_SOURCE"})){$SYMBOL_SOURCE=$CSQfield_derived_tok_tmp[$HI1_vep{"SYMBOL_SOURCE"}];unless (defined($SYMBOL_SOURCE)){$SYMBOL_SOURCE="NaN";}} else{ print "ERROR in SYMBOL_SOURCE\n";}
													if (exists($HI1_vep{"Protein_position"})){$Protein_position=$CSQfield_derived_tok_tmp[$HI1_vep{"Protein_position"}];unless (defined($Protein_position)){$Protein_position="NaN";}} else{ print "ERROR in Protein_position \n";}
													if (exists($HI1_vep{"Allele"})){$Allele=$CSQfield_derived_tok_tmp[$HI1_vep{"Allele"}];unless (defined($Allele)){$Allele="NaN";}}else{ print "ERROR in Allele \n";}	
													if (exists($HI1_vep{"Codons"})){$Codons=$CSQfield_derived_tok_tmp[$HI1_vep{"Codons"}];unless (defined($Codons)){$Codons="NaN";}}else{ print "ERROR in Codons \n";}
													if (exists($HI1_vep{"CANONICAL"})){$CANONICAL=$CSQfield_derived_tok_tmp[$HI1_vep{"CANONICAL"}];unless (defined($CANONICAL)){$CANONICAL="NaN";}}else{ print "ERROR in CANONICAL \n";}
													if (exists($HI1_vep{"DISTANCE"})){$DISTANCE=$CSQfield_derived_tok_tmp[$HI1_vep{"DISTANCE"}];unless (defined($DISTANCE)){$DISTANCE="NaN";}}else{ print "ERROR in DISTANCE \n";}
													if (exists($HI1_vep{"ASN_MAF"})){$ASN_MAF=$CSQfield_derived_tok_tmp[$HI1_vep{"ASN_MAF"}];unless (defined($ASN_MAF)){$ASN_MAF="NaN";}}else{ print "ERROR in ASN_MAF \n";}
													if (exists($HI1_vep{"HIGH_INF_POS"})){$HIGH_INF_POS=$CSQfield_derived_tok_tmp[$HI1_vep{"HIGH_INF_POS"}];unless (defined($HIGH_INF_POS)){$HIGH_INF_POS="NaN";}} else{ print "ERROR in HIGH_INF_POS \n";}
													if (exists($HI1_vep{"BIOTYPE"})){$BIOTYPE=$CSQfield_derived_tok_tmp[$HI1_vep{"BIOTYPE"}];unless (defined($BIOTYPE)){$BIOTYPE="NaN";}}else{ print "ERROR in BIOTYPE \n";}
													if (exists($HI1_vep{"CDS_position"})){$CDS_position=$CSQfield_derived_tok_tmp[$HI1_vep{"CDS_position"}];unless (defined($CDS_position)){$CDS_position="NaN";}}else{ print "ERROR in CDS_position \n";}
													if (exists($HI1_vep{"MOTIF_POS"})){$MOTIF_POS=$CSQfield_derived_tok_tmp[$HI1_vep{"MOTIF_POS"}];unless (defined($MOTIF_POS)){$MOTIF_POS="NaN";}}else{ print "ERROR in MOTIF_POS \n";}
													if (exists($HI1_vep{"Feature_type"})){$Feature_type=$CSQfield_derived_tok_tmp[$HI1_vep{"Feature_type"}];unless (defined($Feature_type)){$Feature_type="NaN";}}else{ print "ERROR in Feature_type \n";}
													if (exists($HI1_vep{"AA_MAF"})){$AA_MAF=$CSQfield_derived_tok_tmp[$HI1_vep{"AA_MAF"}];unless (defined($AA_MAF)){$AA_MAF="NaN";}}else{ print "ERROR in AA_MAF \n";}
													# NOT AVAILABLE W/O INTERNET CONNECTION my $HGVSp="NaN";if (exists($HI1_vep{"HGVSp"})){$HGVSp=$CSQfield_derived_tok_tmp[$HI1_vep{"HGVSp"}];if ($HGVSp!~/\w/){$HGVSp="NaN";}}else{ print "ERROR in HGVSp \n";}
													if (exists($HI1_vep{"Feature"})){$Feature=$CSQfield_derived_tok_tmp[$HI1_vep{"Feature"}];unless (defined($Feature)){$Feature="NaN";}}else{ print "ERROR in Feature \n";}	
													if (exists($HI1_vep{"MOTIF_NAME"})){$MOTIF_NAME=$CSQfield_derived_tok_tmp[$HI1_vep{"MOTIF_NAME"}];unless (defined($MOTIF_NAME)){$MOTIF_NAME="NaN";}}else{ print "ERROR in MOTIF_NAME \n";}
													if (exists($HI1_vep{"CCDS"})){$CCDS=$CSQfield_derived_tok_tmp[$HI1_vep{"CCDS"}];unless (defined($CCDS)){$CCDS="NaN";print "A:$CCDS:A";}}else{ print "ERROR in CCDS\n";}	
													if (exists($HI1_vep{"STRAND"})){$STRAND=$CSQfield_derived_tok_tmp[$HI1_vep{"STRAND"}];unless (defined($STRAND)){$STRAND="NaN";}}else{ print "ERROR in STRAND \n";}
													if (exists($HI1_vep{"EUR_MAF"})){$EUR_MAF=$CSQfield_derived_tok_tmp[$HI1_vep{"EUR_MAF"}];unless (defined($EUR_MAF)){$EUR_MAF="NaN";}}else{ print "ERROR in EUR_MAF \n";}
													if (exists($HI1_vep{"EXON"})){$EXON=$CSQfield_derived_tok_tmp[$HI1_vep{"EXON"}];unless (defined($EXON)){$EXON="NaN";}}else{ print "ERROR in EXON \n";}
													if (exists($HI1_vep{"EA_MAF"})){$EA_MAF=$CSQfield_derived_tok_tmp[$HI1_vep{"EA_MAF"}];unless (defined($EA_MAF)){$EA_MAF="NaN";}}else{ print "ERROR in EA_MAF \n";}
													if (exists($HI1_vep{"PolyPhen"})){$PolyPhen=$CSQfield_derived_tok_tmp[$HI1_vep{"PolyPhen"}];unless (defined($PolyPhen)){$PolyPhen="NaN";}}else{ print "ERROR in PolyPhen\n";}
													if (exists($HI1_vep{"MOTIF_SCORE_CHANGE"})){$MOTIF_SCORE_CHANGE=$CSQfield_derived_tok_tmp[$HI1_vep{"MOTIF_SCORE_CHANGE"}];unless (defined($MOTIF_SCORE_CHANGE)){$MOTIF_SCORE_CHANGE="NaN";}}else{ print "ERROR in MOTIF_SCORE_CHANGE \n";}
													if (exists($HI1_vep{"cDNA_position"})){$cDNA_position=$CSQfield_derived_tok_tmp[$HI1_vep{"cDNA_position"}];unless (defined($cDNA_position)){$cDNA_position="NaN";}}else{ print "ERROR in cDNA_position\n";}
													if (exists($HI1_vep{"ENSP"})){$ENSP=$CSQfield_derived_tok_tmp[$HI1_vep{"ENSP"}];unless (defined($ENSP)){$ENSP="NaN";}}else{ print "ERROR in ENSP \n";}
													if (exists($HI1_vep{"Gene"})){$Gene=$CSQfield_derived_tok_tmp[$HI1_vep{"Gene"}];unless (defined($Gene)){$Gene="NaN";}}else{ print "ERROR in Gene\n";}
													if (exists($HI1_vep{"INTRON"})){$INTRON=$CSQfield_derived_tok_tmp[$HI1_vep{"INTRON"}];unless (defined($INTRON)){$INTRON="NaN";}}else{ print "ERROR in INTRON \n";}
													if (exists($HI1_vep{"SYMBOL"})){$SYMBOL=$CSQfield_derived_tok_tmp[$HI1_vep{"SYMBOL"}];unless (defined($SYMBOL)){$SYMBOL="NaN";}}else{ print "ERROR in SYMBOL \n";}
													if (exists($HI1_vep{"SIFT"})){$SIFT=$CSQfield_derived_tok_tmp[$HI1_vep{"SIFT"}];unless (defined($SIFT)){$SIFT="NaN";}}else{ print "ERROR in SIFT\n";}
													if (exists($HI1_vep{"GMAF"})){$GMAF=$CSQfield_derived_tok_tmp[$HI1_vep{"GMAF"}];unless (defined($GMAF)){$GMAF="NaN";}}else{ print "ERROR in GMAF\n";}
													if (exists($HI1_vep{"AFR_MAF"})){$AFR_MAF=$CSQfield_derived_tok_tmp[$HI1_vep{"AFR_MAF"}];unless (defined($AFR_MAF)){$AFR_MAF="NaN";}}else{ print "ERROR in AFR_MAF\n";}
													if (exists($HI1_vep{"DOMAINS"})){$DOMAINS=$CSQfield_derived_tok_tmp[$HI1_vep{"DOMAINS"}];unless (defined($DOMAINS)){$DOMAINS="NaN";}}else{ print "ERROR in DOMAINS \n";}

													$Amino_acids=~s/\//$Protein_position/g;
													
													#~ print  "La línea de salida es: $CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$Consequence\t$Allele\t$Gene\t$Feature\t$Feature_type\t$cDNA_position\t$CDS_position\t$Protein_position\t$Amino_acids\t$Codons\t$Existing_variation\t$CANONICAL\t$PolyPhen\t$EXON\t$INTRON\t$MOTIF_NAME\t$MOTIF_POS\t$HIGH_INF_POS\t$MOTIF_SCORE_CHANGE\t$GMAF\t$SYMBOL\t$SYMBOL_SOURCE\t$BIOTYPE\t$DISTANCE\t$STRAND\t$ENSP\t$PUBMED\t$CCDS\t$AFassessed\t$AA_MAF\t$EA_MAF\t$AFR_MAF\t$AMR_MAF\t$ASN_MAF\t$EUR_MAF\t$DOMAINS\t$CLIN_SIG\t$SIFT\n";
													#~ print OUT_vep_parsed "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$Consequence\t$Allele\t$Gene\t$Feature\t$Feature_type\t$cDNA_position\t$CDS_position\t$Protein_position\t$Amino_acids\t$Codons\t$Existing_variation\t$CANONICAL\t$PolyPhen\t$EXON\t$INTRON\t$MOTIF_NAME\t$MOTIF_POS\t$HIGH_INF_POS\t$MOTIF_SCORE_CHANGE\t$GMAF\t$SYMBOL\t$SYMBOL_SOURCE\t$BIOTYPE\t$DISTANCE\t$STRAND\t$ENSP\t$PUBMED\t$CCDS\t$AFassessed\t$AA_MAF\t$EA_MAF\t$AFR_MAF\t$AMR_MAF\t$ASN_MAF\t$EUR_MAF\t$DOMAINS\t$CLIN_SIG\t$SIFT\n";
													#~ print "La línea de salida es:$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$Consequence\t$Allele\t$Gene\t$Feature\t$Feature_type\t$cDNA_position\t$CDS_position\t$Protein_position\t$Amino_acids\t$Codons\t$Existing_variation\t$CANONICAL\t$PolyPhen\t$EXON\t$INTRON\t$MOTIF_NAME\t$MOTIF_POS\t$HIGH_INF_POS\t$MOTIF_SCORE_CHANGE\t$GMAF\t$SYMBOL\t$SYMBOL_SOURCE\t$BIOTYPE\t$DISTANCE\t$STRAND\t$ENSP\t$PUBMED\t$CCDS\t$AFassessed\t$AA_MAF\t$EA_MAF\t$AFR_MAF\t$AMR_MAF\t$ASN_MAF\t$EUR_MAF\t$DOMAINS\t$CLIN_SIG\t$SIFT\n";
													my $string1=join(";",$Consequence,$SYMBOL,$Feature,$BIOTYPE,$Amino_acids,$AFassessed,$CLIN_SIG,$PolyPhen,$SIFT,$AFassessed,$GMAF,$AA_MAF,$EA_MAF,$AFR_MAF,$AMR_MAF,$ASN_MAF,$EUR_MAF,$DOMAINS,$MOTIF_NAME,$MOTIF_POS);
													my $string2=join("\t",$CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER);
													my $string3=join("\t",$string2,$string1);
													
													$hash_BIAS{$string2}{$SYMBOL}{$Feature}{$Consequence}{$string3}=1;
							
													#~ print "------------------------------------->hash_BIAS\t$string2\t$SYMBOL\t$Feature\t$Consequence\t$string3\n";
												
												}
										}
									}
						}
		}
	}# while

close LIST_2;
}else {print "No se pudo abrir el fichero $list_2 generado por VEP \n";}

#~ exit;

my %hashA=();

$time='['. timestamp(). ']'."\n";
print "Reading_gtf_INTRONS_Plus_NMD_threshold.txt:$time\n";

if(open(INPUT2,$input2))
{
	#~ gtf_INTRONS_Plus_NMD_threshold.txt
	
	#~ ##ENSMUST       ENSMUSG HGNC    CHROM   strand  Flag_monoexon   INTRON_COORDS(joined_by_;_and___)       NMD_Threshold
	#~ ENST00000000233 ENSG00000004059 ARF5    7       +       0       127228620__127229136;127229218__127229538;127229649__127230119;127230192__127231016;127231143__127231266        127231092
	#~ ENST00000000412 ENSG00000003056 M6PR    12      -       0       9094537__9095011;9095139__9096000;9096132__9096396;9096507__9098013;9098181__9098824;9099002__9102083   9095062

	
	while (my $line = <INPUT2>)
	{
		chomp ($line);
		#~ print "LINE:$line:DD\n";
		unless($line=~/^#/)
		{
			my @tmp=split(/\t/,$line);
		
			my $ENSMUST2=$tmp[0];
			my $ENSMUSG2=$tmp[1];
			my $HGNC2=$tmp[2];
			my $CHROM2=$tmp[3];
			my $strand2=$tmp[4];
			my $Flag_monoexon=$tmp[5];
			
			my $INTRON_coordinates=$tmp[6];
			my $NMD_Threshold=$tmp[7];
			
			#~ $hashA{$ENSMUST2}{$ENSMUSG2}{$HGNC2}{$CHROM2}{$strand2}{$Flag_monoexon}{$INTRON_coordinates}{$NMD_Threshold}=1;
			$hashA{$ENSMUST2}{$Flag_monoexon}{$INTRON_coordinates}{$NMD_Threshold}{$strand2}=1;
			#~ print "hashA\t$ENSMUST2\t$ENSMUSG2\t$HGNC2\t$CHROM2\t$strand2\t$Flag_monoexon\t$INTRON_coordinates\t$NMD_Threshold\n";
			#~ exit;
		
		}
		
	}
};

#~ exit;

if(open (OUTPUT,'>'.$OUT_vep_parsed) && open(OUTPUT2,'>'.$vep_NMD))
{
	# BIAS! Here we print only one effect per transcript!
	# 1 stop_gain, > 2 frameshift  > 3 splice_donor > 4 splice_acceptor > 5 missense_variant > 6 synonymous_variant
	
	foreach my $variant_key_tok(sort keys %hash_BIAS)
	{
		foreach my $Gene_name_tok(sort keys %{$hash_BIAS{$variant_key_tok}})
		{
		foreach my $transcript_ID_tok(sort keys %{$hash_BIAS{$variant_key_tok}{$Gene_name_tok}})
		{
			my @Effects_to_one_transcript_tmp=sort keys%{$hash_BIAS{$variant_key_tok}{$Gene_name_tok}{$transcript_ID_tok}};
			
			# HIERARCHY
			
			my $hierarchical_effect="NaN";
			my %hierarchical_hash;
			
			foreach my $tmp_tok(@Effects_to_one_transcript_tmp)
			{
				my @tmp2=split(/\&/,$tmp_tok);
				
				foreach my $effect_tok(@tmp2)
				{
					$hierarchical_hash{$effect_tok}=1;
				}
				
			}
			
			my @Effects_to_one_transcript_tmp2=sort keys%hierarchical_hash;
			
			#~ print "---------->$variant_key_tok\t$Gene_name_tok\t$transcript_ID_tok\t".join("<->",@Effects_to_one_transcript_tmp2)."\t$hierarchical_effect--->";
			
			
			
			# First stop_gained
			
			if(grep ( $_ eq 'stop_gained', @Effects_to_one_transcript_tmp2))
			{
				$hierarchical_effect="stop_gained";
			}
			elsif(grep ( $_ eq 'frameshift_variant', @Effects_to_one_transcript_tmp2))
			{
				$hierarchical_effect="frameshift_variant";
			}
			elsif(grep ( $_ eq 'splice_donor_variant', @Effects_to_one_transcript_tmp2))
			{
				$hierarchical_effect="splice_donor_variant";
			}
			elsif(grep ( $_ eq 'splice_acceptor_variant', @Effects_to_one_transcript_tmp2))
			{
				$hierarchical_effect="splice_acceptor_variant";
			}
			elsif(grep ( $_ eq 'missense_variant', @Effects_to_one_transcript_tmp2))
			{
				$hierarchical_effect="missense_variant";
			}
			elsif(grep ( $_ eq 'synonymous_variant', @Effects_to_one_transcript_tmp2))
			{
				$hierarchical_effect="synonymous_variant";
			}
			#~ print "$hierarchical_effect\n";
			my $Flag_WARNING="NaN";
			
			# Checking steps
			
			if(exists($hashA{$transcript_ID_tok}))
			{
				foreach my $Flag_monoexon_tok(sort keys%{$hashA{$transcript_ID_tok}})
				{
					foreach my $INTRON_coords_tok(sort keys%{$hashA{$transcript_ID_tok}{$Flag_monoexon_tok}})
					{
						foreach my $NMD_Threshold_tok(sort keys%{$hashA{$transcript_ID_tok}{$Flag_monoexon_tok}{$INTRON_coords_tok}})
						{
						foreach my $strand_tok(sort keys%{$hashA{$transcript_ID_tok}{$Flag_monoexon_tok}{$INTRON_coords_tok}{$NMD_Threshold_tok}})
						{
							
							if($Flag_monoexon_tok == 1 && $hierarchical_effect eq "splice_donor_variant")
							{
								print "WARNING\tVARIANT\t$variant_key_tok\tin_transcript\t$transcript_ID_tok\tcauses\t$hierarchical_effect\tbut\ttranscrip_is_monoexon";
								$Flag_WARNING=1;
								
							}
							elsif($Flag_monoexon_tok == 1 && $hierarchical_effect eq "splice_acceptor_variant")
							{
								print "WARNING\tVARIANT\t$variant_key_tok\tin_transcript\t$transcript_ID_tok\tcauses\t$hierarchical_effect\tbut\ttranscrip_is_monoexon";
								$Flag_WARNING=1;
								
							}
							else
							{
								$Flag_WARNING=0;
								
								if($hierarchical_effect eq "stop_gained" || $hierarchical_effect eq "splice_donor_variant" || $hierarchical_effect eq "splice_acceptor_variant")
								{
									# Calculate NMD for stop_gained splice_donor_variant splice_acceptor_varian
									
									my $NMD="NaN";
									
									# First extacte POS of the variant
									
									my @tmp2=split(/\t/,$variant_key_tok);
									
									my $CHROM=$tmp2[0];
									my $POS=$tmp2[1];
									my $REF=$tmp2[3];
									my $ALT=$tmp2[4];
									
									if($NMD_Threshold_tok ne 'NaN') # NMD threahold calculated for multiexon transcripts that can suffer NMD
									{
										if($strand_tok eq '+')
										{
											if($POS < $NMD_Threshold_tok)
											{
												$NMD="NMD_positive";
											}
											else
											{
												$NMD="NMD_negative";
											}
											
										}
										elsif($strand_tok eq '-')
										{
											if($POS > $NMD_Threshold_tok)
											{
												$NMD="NMD_positive";
											}
											else
											{
												$NMD="NMD_negative";
											}
										}
									}
									else
									{
										$NMD="NMD_negative";
									}
									
									print OUTPUT2 "$CHROM\t$Gene_name_tok\t$POS\t$REF\t$ALT\t$hierarchical_effect\t$transcript_ID_tok\t$NMD\n";
									#~ print "---------------------------->$CHROM\t$Gene_name_tok\t$POS\t$REF\t$ALT\t$hierarchical_effect\t$transcript_ID_tok\t$NMD\n";
									#~ exit;
								}
							}
							
						}
						}
					}	
				}
				
				
			}
			else
			{
				print "$transcript_ID_tok not present in ENSEMBL Homo_sapiens.GRCh37.75.gtf\n";
				$Flag_WARNING=1;
			}
			
			
			
			
			# Printing key
			
			if($Flag_WARNING == 0)
			{
				foreach my $string_effect(sort keys %{$hash_BIAS{$variant_key_tok}{$Gene_name_tok}{$transcript_ID_tok}})
				{
				foreach my $string_tok(sort keys %{$hash_BIAS{$variant_key_tok}{$Gene_name_tok}{$transcript_ID_tok}{$string_effect}})
				{
					my @tmp3=split(/\t/,$string_tok);
					
					for(my $i=0;$i<scalar(@tmp3);$i++)
					{
						if($i != scalar(@tmp3)-1)
						{
							print OUTPUT "$tmp3[$i]\t";
							#~ print  "$tmp3[$i]\t";
						}
						elsif($i == scalar(@tmp3)-1)
						{
							my @tmp5=split(/\;/,$tmp3[scalar(@tmp3)-1]);
							
							my $candidate=shift(@tmp5);
							
							unshift(@tmp5,$hierarchical_effect);
							
							print OUTPUT join(";",@tmp5)."\n";
							#~ print  join(";",@tmp5)."\n";
						}
					}
					
					
				}
				}
			}
		}
		}
	}
}


sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
