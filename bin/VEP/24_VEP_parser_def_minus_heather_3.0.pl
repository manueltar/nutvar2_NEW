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


my $input_vcf= $ARGV[0];
my $list_2= $ARGV[1];
my $input2= $ARGV[2];
my $OUT_VEP_parsed=$ARGV[3];
my $VEP_NMD=$ARGV[4];



my %hash1=();
my %hash0=();
my %hash_deletion=();
my %hash_BIAS=();

my $time='['. timestamp(). ']'."\n";
#~ print "Start reading INPUT_vcf:$time\n";

if(open (INPUT_vcf, $input_vcf))
{
## Here we declare the hashes we are going to use to load the fields from the heather lines. 

	while (my $line = <INPUT_vcf>)
	{
		chomp ($line);
		#~ print "-------------->$line\n";
		if ($line !~ /^#/)
		{
			my @tmp=split(/\t/,$line);
			
			my $CHROM=$tmp[0];
			my $POS=$tmp[1];
			my $ID=$tmp[2];
			my $REF=$tmp[3];
			my $ALT=$tmp[4];
			
			#~ print "****$CHROM\t$ID\t$POS\t$REF\t$ALT\n";
			if(length($REF) <= length($ALT))
			{
				$hash0{$CHROM}{$POS}{$ALT}{$REF}{$ID}=1;
				#~ print"hash0\t$CHROM\t$POS\t$ALT\t$REF\t$ID\n";
				#~ exit;
			}
			elsif(length($REF) > length($ALT))
			{
				$hash_deletion{$CHROM}{$POS}{$ALT}{$REF}{$ID}=1;
				#~ print"hash_deletion\t$CHROM\t$POS\t$ALT\t$REF\t$ID\n";
				#~ exit;
			}
			
		}
	}
}


$time='['. timestamp(). ']'."\n";
#~ print "Start reading variant_effect_output.txt:$time\n";

if(open (LIST_2, $list_2))
{

#~ #Uploaded_variation     Location        Allele  Gene    Feature Feature_type    Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids     Codons  Existing_variation      Extra
#~ DummyID 1:231487315-231487316   AGGT    ENSG00000010072 ENST00000391858 Transcript      frameshift_variant      1996-1997       716-717 239     K/KGX   aaa/aaAGGTa     -       IMPACT=HIGH;STRAND=1;VARIANT_CLASS=insertion;SYMBOL=SPRTN;SYMBOL_SOURCE=HGNC;HGNC_ID=25356;BIOTYPE=protein_coding;CCDS=CCDS31054.1;ENSP=ENSP00000375731;SWISSPROT=SPRTN_HUMAN;UNIPARC=UPI000004C614;EXON=4/4;HGVSc=ENST00000391858.4:c.718_721dupGGTA;HGVSp=ENSP00000375731.4:p.Thr241ArgfsTer18;HGVS_OFFSET=5
#~ DummyID 1:231487315-231487316   AGGT    ENSG00000010072 ENST00000295050 Transcript      frameshift_variant,splice_region_variant        1052-1053       716-717 239     K/KGX   aaa/aaAGGTa     -       IMPACT=HIGH;STRAND=1;VARIANT_CLASS=insertion;SYMBOL=SPRTN;SYMBOL_SOURCE=HGNC;HGNC_ID=25356;BIOTYPE=protein_coding;CANONICAL=YES;CCDS=CCDS1594.1;ENSP=ENSP00000295050;SWISSPROT=SPRTN_HUMAN;TREMBL=L8E708_HUMAN;UNIPARC=UPI000006D601;EXON=4/5;DOMAINS=hmmpanther:PTHR21220;HGVSc=ENST00000295050.7:c.718_718+3dupGGTA;HGVS_OFFSET=5

	while (my $line = <LIST_2>)
	{
	chomp ($line);
		#~ print "$line\n";
		
		if ($line !~ /^#/)
		{
			#~ print "Hello_worldI:$line\n";
			my @tmp = split (/\t/,$line);
			#~ print "El array es:\t".join("<->",@tmp)."\n";
			
			my $ID=$tmp[0];
			my $Location=$tmp[1];
			
			my @tmp_location=split(/\:/,$Location);
			
			my $CHROM=$tmp_location[0];
			my $POS_string=$tmp_location[1];
			my $POS=$POS_string;
			$POS =~ s/\-\d+//g;
			my $ALT=$tmp[2];
			my $ENSMUST=$tmp[4];
			my $Consequence_string=$tmp[6];
			my $INFO=$tmp[13];
			
			# Consequences array # Prioritization 1 csq per transcript
			
			my @tmp_csq=split(/\,/,$Consequence_string);
			
			my $csq_def="NaN";
			
			if(grep($_ eq 'stop_gained',@tmp_csq))
			{
				$csq_def="stop_gained";
			}
			elsif(grep($_ eq 'frameshift_variant',@tmp_csq))
			{
				$csq_def="frameshift_variant";
			}
			elsif(grep($_ eq 'splice_donor_variant',@tmp_csq))
			{
				$csq_def="splice_donor_variant";
			}
			elsif(grep($_ eq 'splice_acceptor_variant',@tmp_csq))
			{
				$csq_def="splice_acceptor_variant";
			}
			elsif(grep($_ eq 'missense_variant',@tmp_csq))
			{
				$csq_def="missense_variant";
			}
			elsif(grep($_ eq 'synonymous_variant',@tmp_csq))
			{
				$csq_def="synonymous_variant";
			}
			
			# SYMBOL array
			
			my @tmp_SYMBOL=split(/\;SYMBOL=/,$INFO);
			
			my $SYMBOL=$tmp_SYMBOL[1];
			$SYMBOL=~ s/\;.+//g;
			
			# BIOTYPE array
			my @tmp_BIOTYPE=split(/\;BIOTYPE=/,$INFO);
			my $Transcript_BioType=$tmp_BIOTYPE[1];
			$Transcript_BioType=~ s/\;.+//g;
			
			# AA change
			
			my $Amino_Acid_Change=$tmp[10];
			
			####
			
			my $string1=join(";",$csq_def,$SYMBOL,$ENSMUST,$Transcript_BioType,$Amino_Acid_Change);
			
			#~ print "----1--->$CHROM\t$POS_string\t$POS\t$ALT\t$ENSMUST\t$csq_def\t".join("<->",@tmp_csq)."\t$SYMBOL\t$Transcript_BioType\t$Amino_Acid_Change\n";
			#~ print "----10--->$ENSMUST\t$Consequence_string\t$INFO\n";
			#~ exit;
			
			if(exists($AcceptedChr{$CHROM}) && exists($AcceptedEffects{$csq_def}))
			{
				#~ print "----2--->$CHROM\t$POS\t->$ALT<-\n";
				
				if($ALT eq '-')
				{
					# This is an deletion
					
					#~ print "----deletion--->$CHROM\t$POS\t$ALT\t$SYMBOL\t$ENSMUST\t$csq_def\n";
					
					my $POS_adjusted=$POS-1;
					
					foreach my $ALT_tok(sort keys %{$hash_deletion{$CHROM}{$POS_adjusted}})
					{
					foreach my $REF_tok(sort keys %{$hash_deletion{$CHROM}{$POS_adjusted}{$ALT_tok}})
					{
						#~ print "--------------------------------------------------deletion--->$CHROM\t$POS_adjusted\t$REF_tok\t$ALT_tok\t$SYMBOL\t$ENSMUST\t$csq_def\n";
						
						my $QUAL=100;
						my $FILTER="PASS";
						
						my $string2=join("\t",$CHROM,$POS_adjusted,$ID,$REF_tok,$ALT_tok,$QUAL,$FILTER);
						my $string3=join("\t",$string2,$string1);
						
						my @tmp_DEL=split(/\-/,$POS_string);
						my $POS_END=$tmp_DEL[1];
						my $offset=$POS_END - $POS_adjusted +1;
						
						#~ print "--->$POS_adjusted\n";
						#~ print "--->$POS_END\n";
						#~ print "--->$REF_tok\n";
						
						
						if($offset eq length($REF_tok))
						{
							$hash_BIAS{$string2}{$SYMBOL}{$ENSMUST}{$csq_def}{$string3}=1;
							#~ print "hash_BIAS\t$string2\t$SYMBOL\t$ENSMUST\t$csq_def\t$string3\n";
						}
						
					}
					}
				}
				elsif(exists($hash0{$CHROM}{$POS}{$ALT}))
				{
					# This is a missense/ synonimous
					foreach my $REF_tok(sort keys %{$hash0{$CHROM}{$POS}{$ALT}})
					{
						#~ print "--------------------------------------------------missense/syn--->$CHROM\t$POS\t$REF_tok\t$ALT\t$SYMBOL\t$ENSMUST\t$csq_def\n";
						
						my $QUAL=100;
						my $FILTER="PASS";
						my $string2=join("\t",$CHROM,$POS,$ID,$REF_tok,$ALT,$QUAL,$FILTER);
						my $string3=join("\t",$string2,$string1);
						$hash_BIAS{$string2}{$SYMBOL}{$ENSMUST}{$csq_def}{$string3}=1;
						
						#~ print "hash_BIAS\t$string2\t$SYMBOL\t$ENSMUST\t$csq_def\t$string3\n";
					}
				}
				elsif(exists($hash0{$CHROM}{$POS}))
				{
					#~ print "----4--->$CHROM\t$POS\t$ALT\n";
					foreach my $ALT_tok(sort keys %{$hash0{$CHROM}{$POS}})
					{
					foreach my $REF_tok(sort keys %{$hash0{$CHROM}{$POS}{$ALT_tok}})
					{
						my @insertion_point=split(/$ALT/,$ALT_tok);
						
						my $REF=$insertion_point[0];
						#~ print "SCALAR_insertion_poin\t".scalar(@insertion_point)."\t$REF\n";
						
						if($REF_tok eq $REF)
						{
							#~ print "--------------------------------------------------insertion--->$CHROM\t$POS\t$REF_tok\t$ALT_tok\t$SYMBOL\t$ENSMUST\t$csq_def\n";
							
							# This is an insertion
							my $QUAL=100;
							my $FILTER="PASS";
							my $string2=join("\t",$CHROM,$POS,$ID,$REF_tok,$ALT_tok,$QUAL,$FILTER);
							my $string3=join("\t",$string2,$string1);
							$hash_BIAS{$string2}{$SYMBOL}{$ENSMUST}{$csq_def}{$string3}=1;
							
							#~ print "hash_BIAS\t$string2\t$SYMBOL\t$ENSMUST\t$csq_def\t$string3\n";
						}
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
#~ print "Reading_gtf_INTRONS_Plus_NMD_threshold.txt:$time\n";

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
}
if(open (OUTPUT,'>'.$OUT_VEP_parsed) && open(OUTPUT2,'>'.$VEP_NMD))
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
			
			# First stop_gained
			
			if(grep ( $_ eq 'stop_gained', @Effects_to_one_transcript_tmp))
			{
				$hierarchical_effect="stop_gained";
			}
			elsif(grep ( $_ eq 'frameshift_variant', @Effects_to_one_transcript_tmp))
			{
				$hierarchical_effect="frameshift_variant";
			}
			elsif(grep ( $_ eq 'splice_donor_variant', @Effects_to_one_transcript_tmp))
			{
				$hierarchical_effect="splice_donor_variant";
			}
			elsif(grep ( $_ eq 'splice_acceptor_variant', @Effects_to_one_transcript_tmp))
			{
				$hierarchical_effect="splice_acceptor_variant";
			}
			elsif(grep ( $_ eq 'missense_variant', @Effects_to_one_transcript_tmp))
			{
				$hierarchical_effect="missense_variant";
			}
			elsif(grep ( $_ eq 'synonymous_variant', @Effects_to_one_transcript_tmp))
			{
				$hierarchical_effect="synonymous_variant";
			}
			
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
				foreach my $string_tok(sort keys %{$hash_BIAS{$variant_key_tok}{$Gene_name_tok}{$transcript_ID_tok}{$hierarchical_effect}})
				{
					print OUTPUT "$string_tok\n";
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
