#!/bin/bash
#SBATCH --account=ctb-hussinju
#SBATCH --time=120:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=1
#SBATCH --output=final_index_POPRES_tmp_QC_maf_hwe_ld_geno_common_CAG_CPTP_POPRES.out

path_common=/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES
path_merged_JH=/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_all_arrays_09122021/maf0.01.hwe0.000001.miss5perc/Imputed/merging/noArchi
path_merged_CAG=/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_from_CaG_TopMed/merged/noArchi



module load tabix/0.2.6
module load plink/1.9b_6.21-x86_64

rm $path_common/gp_09_files_to_merge_plink.txt
#rm $path_common/files_to_merge_plink.txt
#THIS PART WAS SPLIT IN multiple BATCH FOR FASTER ACCESS.
i='JH_imputation' 
echo $i
#plink --vcf $path_merged_JH/JH_imputation_CAG_CPTP_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.0_JH_merge.vcf --double-id --keep-allele-order --vcf-min-gp 0.9 --make-bed --out $path_common/gp_09_CAG_CPTP_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$i


i='CAG_imputation'
echo $i
#plink --vcf $path_merged_CAG/CAG_imputation_CAG_CPTP_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.0_CAG_merge.vcf --double-id --keep-allele-order --vcf-min-gp 0.9 --make-bed --out $path_common/gp_09_CAG_CPTP_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$i	

	



echo "dataset for JH imputation"

	plink --bfile $path_common/gp_09_CAG_CPTP_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_JH_imputation --merge-list $path_common/JH_imputation_gp_09_files_to_merge_plink.txt --double-id --keep-allele-order --make-bed --out $path_common/JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL

	plink --bfile $path_common/JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --double-id --keep-allele-order --maf 0.01 --hwe 0.00001 --geno 0.01 --make-bed --out $path_common/filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL 

	plink --bfile $path_common/filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL  --double-id --indep-pairwise 1000 50 0.05 --out $path_common/tmp_filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL
	
	plink --bfile $path_common/filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --double-id --extract $path_common/tmp_filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL.prune.in --make-bed --out $path_common/pruned_filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL 
	#plink --bfile $path_common/filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --double-id --extract $path_common/tmp_filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL.prune.in --make-bed --out $path_common/pruned_filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL 




echo "dataset for CAG imputation"

	#plink --bfile $path_common/gp_09_CAG_CPTP_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_CAG_imputation --merge-list $path_common/CAG_imputation_gp_09_files_to_merge_plink.txt --double-id --keep-allele-order --make-bed --out $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL

        #plink --bfile $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --double-id --keep-allele-order --maf 0.01 --hwe 0.00001 --geno 0.01 --make-bed --out $path_common/filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL

        #plink --bfile $path_common/filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL  --double-id --indep-pairwise 1000 50 0.05 --out $path_common/tmp_filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL

        #plink --bfile $path_common/filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --double-id --extract $path_common/tmp_filtered_JH_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL.prune.in --make-bed --out $path_common/pruned_filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL
        #plink --bfile $path_common/filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --double-id --extract $path_common/tmp_filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL.prune.in --make-bed --out $path_common/pruned_filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL



rm $path_common/*.nosex




cohort='1000G'
echo "extract all ind from $cohort"

plink --bfile $path_common/pruned_filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --keep $path_common/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort.sample.ID --make-bed --out $path_common/JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort
#plink --bfile $path_common/pruned_filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --keep $path_common/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort.sample.ID --make-bed --out $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort




#do every cohort togheter
cohort='QC_REF'
echo "extract all ind from $cohort"

plink --bfile $path_common/pruned_filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --keep $path_common/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort.sample.ID --make-bed --out $path_common/JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort
#plink --bfile $path_common/pruned_filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --keep $path_common/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort.sample.ID --make-bed --out $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort




#do every cohort togheter
cohort='CAG_all'
echo "extract all ind from $cohort"

cut -f1,2 $path_common/gp_09_CAG_CPTP_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_JH_imputation.fam > $path_common/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_JH_imputation.sample.ID
cut -f1,2 $path_common/gp_09_CAG_CPTP_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_CAG_imputation.fam > $path_common/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_CAG_imputation.sample.ID



plink --bfile $path_common/pruned_filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --keep $path_common/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_JH_imputation.sample.ID --make-bed --out $path_common/JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort
plink --bfile $path_common/pruned_filtered_JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --keep $path_common/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_CAG_imputation.sample.ID --make-bed --out $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort
#plink --bfile $path_common/pruned_filtered_CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_ALL --keep $path_common/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_CAG_imputation.sample.ID --make-bed --out $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$cohort




echo "sbatch flash PCA"
sbatch /lustre07/scratch/justinp/topmed_new/imputation_pipeline/1000G_PCA/final_flash_PCA_CAG_CPTP_QCREF_POPRES_1000G.sh


echo "end"







