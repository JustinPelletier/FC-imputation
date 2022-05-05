#!/bin/bash
#SBATCH --account=ctb-hussinju
#SBATCH --time=20:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --output=final_flash_PCA_CAG_QCREF_genotyped.out

module load vcftools
module load tabix
module load bcftools

path=/lustre07/scratch/justinp/topmed_new/QC_REF
path_common=/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/genotyped
#path_common=/lustre07/scratch/justinp/topmed_new/Common_positions_imputed_CAG_CPTP_QCREF_POPRES


#DO a PCA on 1000G WGS sequences

cd $path_common


echo "QC REF"
/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/final_genotyped_CAG_QC_REF_common_hg19_QC_REF --suffix _QC_REF.txt --outload $path_common/PCA/loadings_QC_REF.txt --outmeansd $path_common/PCA/meansd_QC_REF.txt


mv *.txt $path_common/PCA





chip='CAG_euro'
#chip='CAG'
#do a projection on POPRES


echo "CAG projection QC REF"
/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/final_genotyped_CAG_QC_REF_common_hg19_$chip --project --inmeansd $path_common/PCA/meansd_QC_REF.txt --outproj $path_common/PCA/CAG_projections_QC_REF_$chip.txt --inload $path_common/PCA/loadings_QC_REF.txt -v











echo "end"
