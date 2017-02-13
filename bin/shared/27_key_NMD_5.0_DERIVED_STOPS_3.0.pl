##	Script to calculate if the frameshift derived first stop-codon entails NMD. 2014.
##	Manuel Tardáguila Sancho. Master Student (Supervisor A. Rausell Telenti Lab, CHUV-SIB, Lausanne) 

use strict;
use warnings;
use Time::localtime;

my %hash2=();
my %hash3=();
my %hash_splice=();
my %hash_derived_PTC=();
my %hash_EXON=();
my %RESULT_hash=();
my %RESULT_hash_derived_PTC=();

my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $input3=$ARGV[2];
my $output=$ARGV[3];

my $time='['. timestamp(). ']'."\n";
#~ print "Reading_gtf_INTRONS_Plus_NMD_threshold.txt:$time\n";

my %hash0=();

if(open (INPUT1,$input1))
{
	#~ gtf_INTRONS_Plus_NMD_threshold.txt
	
	#~ ##ENSMUST       ENSMUSG HGNC    CHROM   strand  Flag_monoexon   INTRON_COORDS(joined_by_;_and___)       NMD_Threshold
	#~ ENST00000000233 ENSG00000004059 ARF5    7       +       0       127228620__127229136;127229218__127229538;127229649__127230119;127230192__127231016;127231143__127231266        127231092
	#~ ENST00000000412 ENSG00000003056 M6PR    12      -       0       9094537__9095011;9095139__9096000;9096132__9096396;9096507__9098013;9098181__9098824;9099002__9102083   9095062

	
	while (my $line = <INPUT1>)
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
			
			$hash0{$ENSMUST2}{$NMD_Threshold}{$strand2}=1;
			#~ print "hash0\t$ENSMUST2\t$NMD_Threshold\t$strand2\n";
			#~ exit;
		
		}
	}
}

$time='['. timestamp(). ']'."\n";
#~ print "Start charging hash1:$time\n";

if (open (INPUT2, $input2))
{
	#INPUT2= 1GK_frameshift_derived_PTCs.txt
	
	#  X       GRIA3   122336600       T       TG      frameshift_variant&feature_elongation   ENST00000371264 ->DERIVED_PTC:122338313
#  X       GRIA3   122336600       T       TG      frameshift_variant&feature_elongation   ENST00000371266 ->DERIVED_PTC:122338313
#  X       IL9R    155231142       GT      G       frameshift_variant&feature_truncation   ENST00000369423 ->DERIVED_PTC:155232595
#  X       IL9R    155231142       GT      G       frameshift_variant&feature_truncation   ENST00000540897 ->DERIVED_PTC:155232595
#  X       LCA10   153149707       C       CG      frameshift_variant&feature_elongation   ENST00000357566 ->DERIVED_PTC:153149717
#  X       LCA10   153149708       G       GC      frameshift_variant&feature_elongation   ENST00000357566 ->DERIVED_PTC:153149717
#  X       LCA10   153151280       G       GCC     frameshift_variant&feature_elongation   ENST00000357566 ->DERIVED_PTC:153152479
#  X       LCA10   153151284       GT      G       frameshift_variant&feature_truncation   ENST00000357566 ->DERIVED_PTC:153152479

	
while(my $line=<INPUT2>)
	{
		chomp $line;
		#print "$line\n";
		#~ print "La línea es:$line\n";
		if($line=~/>DERIVED_PTC:/)
		{
			#~ print "La línea es:$line\n";
			
			my @tmp=split("\t->DERIVED_PTC:",$line);
			my $parental=$tmp[0];
			my $derived=$tmp[1];
			my ($CHROM,$SYMBOL,$POS,$REF,$ALT,$Effect,$ENST,$POS_derived)="NaN"x8;
			if($parental=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
			{
				$CHROM=$1;
				$SYMBOL=$2;
				$POS=$3;
				$REF=$4;
				$ALT=$5;
				$Effect=$6;
				$ENST=$7;
				
			}
			if($derived=~/(.+)/)
			{
				$POS_derived=$1;
			}
			
			my $NMD="NaN";
									
			# First extacte POS of the variant
			
			if(exists($hash0{$ENST}))
			{
				foreach my $NMD_Threshold_tok(sort keys%{$hash0{$ENST}})
				{
					foreach my $strand_tok(sort keys%{$hash0{$ENST}{$NMD_Threshold_tok}})
					{
						if($NMD_Threshold_tok ne 'NaN') # NMD threahold calculated for multiexon transcripts that can suffer NMD
						{
							if($strand_tok eq '+')
							{
								if($POS_derived < $NMD_Threshold_tok)
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
								if($POS_derived > $NMD_Threshold_tok)
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
					}
				}
			}
			else
			{
				print "WARNING\t$ENST\tnot present in ENSEMBL Homo_sapiens.GRCh37.75.gtf\n"; #checked 0
			}
			$hash_derived_PTC{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}{$POS_derived}{$NMD}=1;
			
			#~ print "hash_derived_PTC_1\t$CHROM\t$SYMBOL\t$ENST\t$POS\t$REF\t$ALT\t$Effect\t$POS_derived\t$NMD\n";
			#~ exit;
		}
		elsif($line=~/->STOP_LOST/)
		{
			
			my @tmp=split("\t-->STOP_LOST",$line);
			my $parental=$tmp[0];
			my $derived=$tmp[1];
			my ($CHROM,$SYMBOL,$POS,$REF,$ALT,$Effect,$ENST,$POS_derived)="NaN"x8;
			if($parental=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
			{
				$CHROM=$1;
				$SYMBOL=$2;
				$POS=$3;
				$REF=$4;
				$ALT=$5;
				$Effect=$6;
				$ENST=$7;
				
			}
			
			$POS_derived=$POS;
			
			#~ # All frameshift resulting in an stop outside the CDS are assumed to be NMD_negative
			#~ 
			my $NMD="NaN";
			#~ 
			$hash_derived_PTC{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}{$POS_derived}{$NMD}=1;
			#~ print "hash_derived_PTC_2\t$CHROM\t$SYMBOL\t$ENST\t$POS\t$REF\t$ALT\t$Effect\t$POS_derived\t$NMD\n";
			#~ exit;
		}	
	}
}else{print "Unable to open INPUT2";}

#~ exit;

$time='['. timestamp(). ']'."\n";
#~ print "Start charging hash2:$time\n";

if(open (INPUT3, $input3) && open (OUTPUT,'>'.$output))
{
	#INPUT1=XXX_out_vep_parsed.vcf
	
	#	13      32906565        .       CA      C       33.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000380152;protein_coding;;0.5;;;;0.5;;;;NaN;NaN;NaN;NaN;PIRSF_domain:PIRSF002397;;
	#	13      32906565        .       CA      C       33.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000530893;protein_coding;;0.5;;;;0.5;;;;NaN;NaN;NaN;NaN;NaN;;
	#	13      32906565        .       CA      C       33.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000544455;protein_coding;;0.5;;;;0.5;;;;NaN;NaN;NaN;NaN;PIRSF_domain:PIRSF002397;;
	#	13      32906576        .       CA      C       36.73   .       frameshift_variant&feature_truncation;BRCA2;ENST00000380152;protein_co
	
	#INPUT1b=XXX-eff_parsed.vcf
	
	#13      32906565        .       CA      C       33.73   .       frameshift_variant;BRCA2;ENST00000380152;protein_coding;T317;10;;CODING;aca/;3418;HIGH;1;0.5;NaN;NaN    LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        
	#13      32906565        .       CA      C       33.73   .       frameshift_variant;BRCA2;ENST00000530893;protein_coding;T194;10;;CODING;aca/;480;HIGH;1;0.5;WARNING_TRANSCRIPT_INCOMPLETE;NaN   LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        
	#13      32906565        .       CA      C       33.73   .       frameshift_variant;BRCA2;ENST00000544455;protein_coding;T317;10;;CODING;aca/;3418;HIGH;1;0.5;NaN;NaN    LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        
	#13      32906576        .       CA      C       36.73   .       frameshift_variant;BRCA2;ENST00000380152;protein_coding;Q321;10;;CODING;caa/;3418;HIGH;1;0.5;NaN;NaN    LOF_positive_snpEff;BRCA2;ENSG00000139618;6;0.50        

while(my $line=<INPUT3>)
	{
		chomp $line;
		#print "$line\n";
		if($line=~/^([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t(.+)/)
		{
			my $CHROM=$1;
			my $POS=$2;
			my $REF=$3;
			my $ALT=$4;
			my $fields=$5;
			#~ print "hello_world:$CHROM\t$POS\t$fields\t$REF\t$ALT\n";
			#~ exit;
			my @fields_tmp=split(";",$fields);
			my $Effect=$fields_tmp[0];
			my $SYMBOL=$fields_tmp[1];
			my $ENST=$fields_tmp[2];
			if(exists($hash_derived_PTC{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}))
			{
				foreach my $POS_derived_tok(sort{$a<=>$b} keys%{$hash_derived_PTC{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}})
				{
					foreach my $NMD_tok(sort keys%{$hash_derived_PTC{$CHROM}{$SYMBOL}{$ENST}{$POS}{$REF}{$ALT}{$Effect}{$POS_derived_tok}})
					{
						print OUTPUT "$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$ENST\t->DERIVED:$POS_derived_tok\t$NMD_tok\n";
					}
				}
				
			}
			else
			{
				my $POS_derived_tok=$POS;
				my $NMD_tok="NaN";
				
				print OUTPUT "$CHROM\t$SYMBOL\t$POS\t$REF\t$ALT\t$Effect\t$ENST\t->DERIVED:$POS_derived_tok\t$NMD_tok\n";
			}
			
		}
	}
}else{print "Unable to open INPUT1\n";}



sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
