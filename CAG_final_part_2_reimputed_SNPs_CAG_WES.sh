#!/bin/bash
#SBATCH --account=ctb-hussinju
#SBATCH --time=120:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=1
#SBATCH --output=CAG_final_part_2_reimputed_SNPs_CAG_WES.out


path_new=/lustre07/scratch/justinp/topmed_new/CAG_new/
path_CAG=/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_from_CaG_TopMed/match_with_Exome_samples
path_old=/lustre04/scratch/justinp/topmed_new/CAG_old
path_imputation=/lustre06/project/6065672/shared/Cartagene/Exome/results/SNP_calling_Step2/GRCh38


module load plink
module load vcftools
module load bcftools

#source /project/6005588/shared/virtualenv_python_3.7.0/bin/activate


#more info https://github.com/hmgu-itg/VCF-liftover
#loop through the directory of each chip
for chip in  'CAG' #'gsa_760' 'gsa_4224' 'gsa_5300' #'gsa_17K' #'axiom' 'omni' 'gsa_17K' 
do
	echo $chip

		rm $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP.R2
	#loop through each chromosome
	 	
		#split by chr
		#uncomment for remake first line to split chr
		#zgrep -w "chr" $path_imputation/with_indels_good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz |sed 's/ /\t/g' | cut -f 3 > $path_new/WES/CAG_MERGE/WES_filtered_out_chr.ID
		#FOR NO INDEL VERSION
		echo "NO INDELS"
		zgrep -v "^#" $path_imputation/good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz |sed 's/ /\t/g' | cut -f 3 > $path_new/WES/CAG_MERGE/${chip}_no_indels_WES_filtered_out_chr.ID
			
			
		#zgrep -w "chr" $path_imputation/good_ID_grch38_102_96.CaG.Exomes.snp.recalibrated.PASS.miss0.5.recode.vcf.gz |sed 's/ /\t/g' | cut -f 1,2 > $path_new/$chip/CAG_filtered_out_chr.bed
		#mv  $path_new/$chip/CAG_filtered_out_chr.bed  $path_new/$chip/CAG_filtered_out_chr.bed_tmp
		#mv  $path_new/$chip/CAG_filtered_out_chr.bed  $path_new/$chip/CAG_filtered_out_chr.bed_tmp
		#sed -i 's/ /\t/g' $path_new/$chip/CAG_filtered_out_chr.bed_tmp | cut -f 1,2 > $path_new/$chip/CAG_filtered_out_chr.bed
		#rm $path_new/$chip/CAG_filtered_out_chr.bed_tmp
		
		
		#find the re-imputed_sites
		echo "find reimputed sites in $chip"
		#zgrep -wf $path_new/WES/CAG_MERGE/WES_filtered_out_chr.ID $path_new/$chip/chr.dose.vcf.gz  > $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP.vcf
		#no indels version
		zgrep "^#"  $path_CAG/sample_id_changed_Cartagene.CaG_allCHR.allSamples.matchExomesSamples.vcf.gz > $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP.vcf
		zgrep -wf $path_new/WES/CAG_MERGE/${chip}_no_indels_WES_filtered_out_chr.ID $path_CAG/sample_id_changed_Cartagene.CaG_allCHR.allSamples.matchExomesSamples.vcf.gz >> $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP.vcf
		#zgrep -wf $path_new/WES/CAG_MERGE/WES_filtered_out_chr.ID $path_new/$chip/chr.dose.vcf.gz | grep "IMPUTED" > $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP.vcf
		#zgrep -wf $path_new/WES/CAG_MERGE/WES_filtered_out_chr.ID $path_new/$chip/chr.dose.vcf.gz | grep "IMPUTED" > $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP.vcf
		#vcftools --gzvcf $path_new/$chip/chr.dose.vcf.gz --positions $path_new/$chip/CAG_filtered_out_chr.bed --recode --out $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr	
		#vcftools --gzvcf $path_new/$chip/chr.dose.vcf.gz --positions $path_new/$chip/CAG_filtered_out_chr.bed --recode --out $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr	
		#mv $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr.recode.vcf $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP.vcf
		
		echo "creates info files for common variants"
		cut -f 1,2,3,8  $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP.vcf > $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr.info
		cut -f 1,2,3 $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr.info > $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr.bed
		cut -d ";" -f 3 $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr.info > $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr.tmp
		paste $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr.bed $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr.tmp >  $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr.R2
		
		sed -i 's/R2=//g' $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP_chr.R2
		echo "bgzip final file"
		bgzip -f $path_new/WES/CAG_MERGE/${chip}_re_imputed_SNP.vcf
done

#echo "sbatch gzip"
#sbatch gzip_CAG_new.sh
echo "end"
