#!/bin/bash
#SBATCH --account=ctb-hussinju
#SBATCH --time=20:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --output=final_flash_PCA_CAG_CPTP_QCREF_POPRES_1000G.out

module load vcftools
module load tabix
module load bcftools

path=/lustre07/scratch/justinp/topmed_new/QC_REF
path_common=/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES
#path_common=/lustre07/scratch/justinp/topmed_new/Common_positions_imputed_CAG_CPTP_QCREF_POPRES


#path flashpca
#~/projects/def-hussinju/shared/bin/flashpca_x86-64
#doccumentation :   https://github.com/gabraham/flashpca#flashpcaR


#DO a PCA on 1000G WGS sequences
echo "1000G"
/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_1000G --suffix _JH_1000G.txt --outload $path_common/PCA/loadings_JH_1000G.txt --outmeansd $path_common/PCA/meansd_JH_1000G.txt
#/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_000G --suffix _CAG_1000G.txt --outload $path_common/PCA/loadings_CAG_1000G.txt --outmeansd $path_common/PCA/meansd_CAG_1000G.txt




echo "QC REF"
/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_QC_REF --suffix _JH_QC_REF.txt --outload $path_common/PCA/loadings_JH_QC_REF.txt --outmeansd $path_common/PCA/meansd_JH_QC_REF.txt
#/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_QC_REF --suffix _CAG_QC_REF.txt --outload $path_common/PCA/loadings_CAG_QC_REF.txt --outmeansd $path_common/PCA/meansd_CAG_QC_REF.txt



#mv eigenvalues_POPRES.txt pve_POPRES.txt pcs_POPRES.txt eigenvectors_POPRES.txt /lustre07/scratch/justinp/topmed_new/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA
mv *.txt $path_common/PCA





chip='CAG_all'
#chip='CAG'
#do a projection on POPRES
echo "CAG projection 1000G"
/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$chip --project --inmeansd $path_common/PCA/meansd_JH_1000G.txt --outproj $path_common/PCA/JH_projections_JH_1000G_$chip.txt --inload $path_common/PCA/loadings_JH_1000G.txt -v
/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$chip --project --inmeansd $path_common/PCA/meansd_JH_1000G.txt --outproj $path_common/PCA/CAG_projections_CAG_1000G_$chip.txt --inload $path_common/PCA/loadings_JH_1000G.txt -v
#/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$chip --project --inmeansd $path_common/PCA/meansd_CAG_1000G.txt --outproj $path_common/PCA/CAG_projections_CAG_1000G_$chip.txt --inload $path_common/PCA/loadings_CAG_1000G.txt -v




echo "CAG projection QC REF"
/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/JH_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$chip --project --inmeansd $path_common/PCA/meansd_JH_QC_REF.txt --outproj $path_common/PCA/JH_projections_JH_QC_REF_$chip.txt --inload $path_common/PCA/loadings_JH_QC_REF.txt -v
/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$chip --project --inmeansd $path_common/PCA/meansd_JH_QC_REF.txt --outproj $path_common/PCA/CAG_projections_CAG_QC_REF_$chip.txt --inload $path_common/PCA/loadings_JH_QC_REF.txt -v
#/lustre06/project/6005588/shared/bin/flashpca_x86-64 --bfile $path_common/CAG_final_filtered_common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_$chip --project --inmeansd $path_common/PCA/meansd_CAG_QC_REF.txt --outproj $path_common/PCA/CAG_projections_CAG_QC_REF_$chip.txt --inload $path_common/PCA/loadings_CAG_QC_REF.txt -v











echo "end"
