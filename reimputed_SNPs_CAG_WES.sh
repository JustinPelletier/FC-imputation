#!/bin/bash
#SBATCH --account=ctb-hussinju
#SBATCH --time=20:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output=reimputed_SNPs_CAG_WES.out


module load vcftools

path_new=/lustre04/scratch/justinp/topmed_new/CAG_new
path_old=/lustre04/scratch/justinp/topmed_new/CAG_old
path_imputation=/lustre06/project/6065672/shared/Cartagene/Exome/results/SNP_calling_Step2/GRCh38


module load plink
module load vcftools
module load bcftools

source /project/6005588/shared/virtualenv_python_3.7.0/bin/activate


#more info https://github.com/hmgu-itg/VCF-liftover
#loop through the directory of each chip
for chip in  'WES' #'omni' 'gsa_4224' 'gsa_5300' 'gsa_760' 'gsa_17K' 
do
	echo $chip
	cd $path_imputation
	
	#correct IDs
        #echo "correct sample ID"
        #bcftools query -l $path_imputation/102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz > $path_imputation/102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz.ID

        #head -102 $path_imputation/102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz.ID > $path_imputation/102.tmp
        #tail -96 $path_imputation/102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz.ID > $path_imputation/96.tmp

        #cut -d "_" -f 1 $path_imputation/102.tmp > $path_imputation/102.tmp.tmp
        #paste $path_imputation/102.tmp.tmp $path_imputation/102.tmp.tmp | sed 's/\t/_/g' > $path_imputation/102.final

        #rev $path_imputation/96.tmp | sed 's/\./_/g'| cut -d "_" -f 1 | rev  > $path_imputation/96.tmp.tmp
        #paste $path_imputation/96.tmp.tmp $path_imputation/96.tmp.tmp | sed 's/\t/_/g' > $path_imputation/96.final

        #cat $path_imputation/102.final $path_imputation/96.final  > $path_imputation/198.final
        #paste $path_imputation/102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz.ID $path_imputation/198.final | sed 's/\t/\//g' > $path_imputation/to_sed.sh.tmp
        #cat $path_imputation/line_CHROM.header $path_imputation/102.final $path_imputation/96.final  > $path_imputation/198.final
        #awk '$1="sed -i \"s/"$1"/g\" $path_imputation/cp_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf"' $path_imputation/to_sed.sh.tmp > $path_imputation/to_sed.sh

	#LAUNCH $path_imputation/to_sed.sh to CHANGE FOR GOOD SAMPLE ID



	#lift over to grch38
	#echo "lift-over"
	#python3 ~/projects/def-hussinju/shared/bin/Crossmap/CrossMap.py vcf /lustre07/scratch/justinp/topmed_new/imputation_pipeline/hg19ToHg38.over.chain.gz $path_imputation/102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz /lustre07/scratch/justinp/topmed_new/imputation_pipeline/hg38.fa $path_imputation/grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf  #--compress
	#bgzip -f $path_imputation/grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf
	
	#cp $path_imputation/102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz $path_imputation/grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz
	#REmove indels
	echo "remove INDELS"
	#this one has the samples ID manually changed
	vcftools --gzvcf $path_imputation/sample_id_changed_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz --remove-indels --recode --out $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5
	#vcftools --gzvcf $path_imputation/102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz --remove-indels --recode --out $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5
	#bgzip -f $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf


	#split lines that have SNP and INDELS on the same position
	echo "split double allelic/indels in two lines"
	for f in $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf ; do f=${f/\.recode\.vcf/} ; bcftools norm -m -any $f.recode.vcf > $f.splitMulti.vcf ; done
	#bgzip -f $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.splitMulti.vcf
	bgzip -f $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf

	
	#remove the indels
	echo "remove indels"
	grep "^#" $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.splitMulti.vcf > $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.splitMulti.recode.vcf
	grep -v "^#" $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.splitMulti.vcf | awk '{if($5!="*")print}' >> $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.splitMulti.recode.vcf
	
	bgzip -f $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.splitMulti.recode.vcf
	bgzip -f $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.splitMulti.vcf

	#change the IDs in the vcf file
	echo "change ID"
	#NO indels version
	zgrep "^#" $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.splitMulti.recode.vcf.gz > $path_imputation/good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf
	zgrep -v "^#" $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.splitMulti.recode.vcf.gz | awk '$3=$1":"$2":"$4":"$5' | sed 's/ /\t/g' >> $path_imputation/good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf
	
	#Invert alleles
	#echo "invert alleles"
	#zgrep -v "^#" $path_imputation/no_indels_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.splitMulti.recode.vcf.gz | awk '$3=$1":"$2":"$5":"$4' | sed 's/ /\t/g' >> $path_imputation/good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf
	
	
	

	#NORMAL
	#zgrep "^#" $path_imputation/grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz > $path_imputation/good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf
	#zgrep -v "^#" $path_imputation/grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz | awk '$3=$1":"$2":"$4":"$5' | sed 's/ /\t/g' >> $path_imputation/good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf
	bgzip -f $path_imputation/good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf

	
	#Get a list of indels ID
	#zgrep -v "^#" $path_imputation/with_indels_good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz |cut -f3 > $path_imputation/with_indels.ID
	
	#zgrep -v "^#" $path_imputation/good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz |cut -f3 > $path_imputation/no_indels.ID

	#zgrep -vf $path_imputation/no_indels.ID $path_imputation/with_indels.ID > $path_imputation/only_indels.ID




done

echo "end"


