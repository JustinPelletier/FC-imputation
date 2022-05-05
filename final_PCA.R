library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(data.table)
library(RColorBrewer)


#CAG_G1000<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/pcs_CAG_1000G.txt", header=T)
JH_G1000<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/pcs_JH_1000G.txt", header=T)

#CAG_QC_REF<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/pcs_CAG_QC_REF.txt", header=T)
JH_QC_REF<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/pcs_JH_QC_REF.txt", header=T)


projection_CAG_1000G<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CAG_projections_CAG_1000G_CAG_all.txt", header=T)
projection_CAG_QC_REF<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CAG_projections_CAG_QC_REF_CAG_all.txt", header=T)

projection_JH_1000G<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/JH_projections_JH_1000G_CAG_all.txt", header=T)
projection_JH_QC_REF<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/JH_projections_JH_QC_REF_CAG_all.txt", header=T)


#GENOTYPED DATA
geno_QC_REF<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/genotyped/PCA/pcs_QC_REF.txt", header=T)

geno_projection_CAG_QC_REF<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/genotyped/PCA/CAG_projections_QC_REF_CAG_euro.txt", header=T)




#-------------------------------LABELS for the cohort/chip----------------------------------------



id_1000G<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_1000G.sample.ID", header=F)
names(id_1000G)<-c("FID", "IID")
id_1000G$cohort<-"1000G"

id_QC_REF<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_QC_REF.sample.ID", header=F)
names(id_QC_REF)<-c("FID", "IID")
id_QC_REF$cohort<-"QC_REF"



#id_axiom<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_axiom.sample.ID", header=F)
#names(id_axiom)<-c("FID", "IID")
#id_axiom$cohort<-"JH_axiom"

#id_omni<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_omni.sample.ID", header=F)
#names(id_omni)<-c("FID", "IID")
#id_omni$cohort<-"JH_omni"

id_gsa_760<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_gsa_760.sample.ID", header=F)
names(id_gsa_760)<-c("FID", "IID")
id_gsa_760$cohort<-"JH_gsa_760"

id_gsa_4224<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_gsa_4224.sample.ID", header=F)
names(id_gsa_4224)<-c("FID", "IID")
id_gsa_4224$cohort<-"JH_gsa_4224"

id_gsa_5300<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_gsa_5300.sample.ID", header=F)
names(id_gsa_5300)<-c("FID", "IID")
id_gsa_5300$cohort<-"JH_gsa_5300"

id_gsa_17K<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/common_1000G_WGS_grcH38_filtered_all_chr_R2_0.8_gsa_17K.sample.ID", header=F)
names(id_gsa_17K)<-c("FID", "IID")
id_gsa_17K$cohort<-"JH_gsa_17K"




CAG_id_gsa_760<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/CAG_imputation_gsa_760.sample.ID", header=F)
names(CAG_id_gsa_760)<-c("FID", "IID")
CAG_id_gsa_760$cohort<-"CAG_gsa_760"

CAG_id_gsa_4224<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/CAG_imputation_gsa_4224.sample.ID", header=F)
names(CAG_id_gsa_4224)<-c("FID", "IID")
CAG_id_gsa_4224$cohort<-"CAG_gsa_4224"

CAG_id_gsa_5300<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/CAG_imputation_gsa_5300.sample.ID", header=F)
names(CAG_id_gsa_5300)<-c("FID", "IID")
CAG_id_gsa_5300$cohort<-"CAG_gsa_5300"

CAG_id_gsa_17k<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/CAG_imputation_gsa_17k.sample.ID", header=F)
names(CAG_id_gsa_17k)<-c("FID", "IID")
CAG_id_gsa_17k$cohort<-"CAG_gsa_17K"



pop_id<-rbind( id_1000G, id_QC_REF, id_gsa_760, id_gsa_4224, id_gsa_5300, id_gsa_17K, CAG_id_gsa_760, CAG_id_gsa_4224, CAG_id_gsa_5300, CAG_id_gsa_17k)


#------------------------------add label---------------------------------------


JH_G1000$population<- pop_id$cohort[match(JH_G1000$IID, pop_id$IID)]

JH_QC_REF$population<- pop_id$cohort[match(JH_QC_REF$IID, pop_id$IID)]

geno_QC_REF$population<-pop_id$cohort[match(geno_QC_REF$IID, pop_id$IID)]

projection_CAG_1000G$population<- pop_id$cohort[match(projection_CAG_1000G$IID, pop_id$IID)]
projection_JH_1000G$population<- pop_id$cohort[match(projection_JH_1000G$IID, pop_id$IID)]

projection_CAG_QC_REF$population<- pop_id$cohort[match(projection_CAG_QC_REF$IID, pop_id$IID)]
projection_JH_QC_REF$population<- pop_id$cohort[match(projection_JH_QC_REF$IID, pop_id$IID)]



geno_projection_CAG_QC_REF$population<-pop_id$cohort[match(geno_projection_CAG_QC_REF$IID, pop_id$IID)]





JH_G1000$cohort<-"1000G"
JH_QC_REF$cohort<-"QC_REF"
geno_QC_REF$cohort<-"QC_REF"

projection_CAG_1000G$cohort<-"CAG"
projection_CAG_QC_REF$cohort<-"CAG"
projection_JH_1000G$cohort<-"CAG"
projection_JH_QC_REF$cohort<-"CAG"
geno_projection_CAG_QC_REF$cohort<-"CAG"

projection_CAG_1000G$Population<-projection_CAG_1000G$population
projection_CAG_QC_REF$Population<-projection_CAG_QC_REF$population
projection_JH_1000G$Population<-projection_JH_1000G$population
projection_JH_QC_REF$Population<-projection_JH_QC_REF$population
geno_projection_CAG_QC_REF$Population<-geno_projection_CAG_QC_REF$population




#-------------------------------POP label-------------------------------------

#incomplete
#pop_POPRES<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/v2_POPRES_ethnicity.info", header=F)
#names(pop_POPRES)<-c("FID", "IID", "Europe", "Ethnicity", "origin", "country" )
#pop_POPRES<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/v3_POPRES_ethnicity.info", header=T)

#missing_info<-POPRES[which(!(POPRES$IID %in% pop_POPRES$IID)),]$IID
#POPRES$population<- pop_POPRES$Ethnicity[match(POPRES$IID, pop_POPRES$IID)]
#POPRES[which(is.na(POPRES$population)),]$population="unknown"
#POPRES$country<- pop_POPRES$country[match(POPRES$IID, pop_POPRES$IID)]
#POPRES[which(is.na(POPRES$country)),]$country="unknown"
#POPRES$origin<- pop_POPRES$origin[match(POPRES$IID, pop_POPRES$IID)]
#POPRES[which(is.na(POPRES$origin)),]$origin="unknown"
#POPRES$Europe<- pop_POPRES$Europe[match(POPRES$IID, pop_POPRES$IID)]
#POPRES[which(is.na(POPRES$Europe)),]$Europe="unknown"

pop_1000G<-fread("/lustre06/project/6005588/shared/References/1000G/1000G_PCA_country.info", sep = "\t", header = TRUE)

JH_G1000$Population<-pop_1000G$Population[match(JH_G1000$IID, pop_1000G$Sample)]



pop_QC_REF<-fread("/lustre06/project/6065672/shared/Quebec_NeuroX_RefDataset/Merged_AllPanels/metadata_for_umap/modified_QC_ref_population.info", sep = " ", header=F)
names(pop_QC_REF)<-c("IID", "Population")

JH_QC_REF$Population<-pop_QC_REF$Population[match(JH_QC_REF$IID, pop_QC_REF$IID)]

geno_QC_REF$Population<-pop_QC_REF$Population[match(geno_QC_REF$IID, pop_QC_REF$IID)]



#------------------------------PCA Ref --------------------------------------------#


data_final_1000G<-JH_G1000

pop_type<-unique(data_final_1000G$Population)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

data_final_1000G$Population<-factor(data_final_1000G$Population,levels=pop_type)



plot_1000G<-ggplot(data_final_1000G, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/1000G_PC1_PC2.png", plot_1000G)
p0<-ggplot(data_final_1000G, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25)) + facet_wrap(~Population) #+xlim(c(-0.1,0.55)) +ylim(c(-0.1,0.5))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/1000G_split_pop_PC1_PC2.png", p0)



for (i in 1:length(pop_type)){

        print(pop_type[i])
        #get the data for one pop
        data_final<- data_final_1000G[which(data_final_1000G$Population==pop_type[i]),]

        #plot the data
        plot<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population)) + geom_point(alpha=0.5)+scale_color_manual(values=col_vector[i])  +xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25))
        ggsave(paste0("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/1000G_", pop_type[i], "_PC1_PC2.png"), plot)
}





#keep only the europeans with 1000G------------------------------------------------#

#define european in 1000G
list_european<-c("GBR", "FIN", "CEU", "IBS","TSI")
data_1000G_euro<-data_final_1000G[which(data_final_1000G$Population %in% list_european),]
PC1_min_euro<-range(data_1000G_euro$PC1)[1]
PC1_max_euro<-range(data_1000G_euro$PC1)[2]

PC2_min_euro<-range(data_1000G_euro$PC2)[1]
PC2_max_euro<-range(data_1000G_euro$PC2)[2]


corresponding_ID<-read.table("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/corresponding_ID_CAG_JH.sample.ID", header=F)
names(corresponding_ID)<-c("CAG_ID", "JH_ID")

projection_CAG_1000G$JH_IID<-corresponding_ID$JH_ID[match(projection_CAG_1000G$IID, corresponding_ID$CAG_ID)]
projection_CAG_QC_REF$JH_IID<-corresponding_ID$JH_ID[match(projection_CAG_QC_REF$IID, corresponding_ID$CAG_ID)]




euro_JH_ID<-projection_JH_1000G[which(projection_JH_1000G$PC1<=PC1_max_euro & projection_JH_1000G$PC1>=PC1_min_euro & projection_JH_1000G$PC2>=PC2_min_euro & projection_JH_1000G$PC2<=PC2_max_euro ),]$IID

euro_CAG_ID<-projection_CAG_1000G[which(projection_CAG_1000G$PC1<=PC1_max_euro & projection_CAG_1000G$PC1>=PC1_min_euro & projection_CAG_1000G$PC2>=PC2_min_euro & projection_CAG_1000G$PC2<=PC2_max_euro ),]$JH_IID

list_strict_european<-intersect(euro_JH_ID, euro_CAG_ID)
df_CAG_ID_euro<-data.frame(list_strict_european, list_strict_european)
#write.table(df_CAG_ID_euro, "/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/genotyped/CAG_euro.sample.ID", row.names=F, col.names=F)


euro_projection_JH_QC_REF<-projection_JH_QC_REF[which(projection_JH_QC_REF$IID %in% list_strict_european),]

euro_projection_CAG_QC_REF<-projection_CAG_QC_REF[which(projection_CAG_QC_REF$JH_IID %in% list_strict_european),]





#----------------------------comparison for the same PC the diff between CAG and JH-------------------- 


merge_CAG_JH_QC_REF<-euro_projection_JH_QC_REF
merge_CAG_JH_QC_REF$PC1_CAG<-euro_projection_CAG_QC_REF$PC1[match(euro_projection_CAG_QC_REF$JH_IID, euro_projection_JH_QC_REF$IID)]
merge_CAG_JH_QC_REF$PC2_CAG<-euro_projection_CAG_QC_REF$PC2[match(euro_projection_CAG_QC_REF$JH_IID, euro_projection_JH_QC_REF$IID)]
merge_CAG_JH_QC_REF$PC3_CAG<-euro_projection_CAG_QC_REF$PC3[match(euro_projection_CAG_QC_REF$JH_IID, euro_projection_JH_QC_REF$IID)]
merge_CAG_JH_QC_REF$PC4_CAG<-euro_projection_CAG_QC_REF$PC4[match(euro_projection_CAG_QC_REF$JH_IID, euro_projection_JH_QC_REF$IID)]
merge_CAG_JH_QC_REF$PC5_CAG<-euro_projection_CAG_QC_REF$PC5[match(euro_projection_CAG_QC_REF$JH_IID, euro_projection_JH_QC_REF$IID)]
merge_CAG_JH_QC_REF$PC6_CAG<-euro_projection_CAG_QC_REF$PC6[match(euro_projection_CAG_QC_REF$JH_IID, euro_projection_JH_QC_REF$IID)]
merge_CAG_JH_QC_REF$PC7_CAG<-euro_projection_CAG_QC_REF$PC7[match(euro_projection_CAG_QC_REF$JH_IID, euro_projection_JH_QC_REF$IID)]
merge_CAG_JH_QC_REF$PC8_CAG<-euro_projection_CAG_QC_REF$PC8[match(euro_projection_CAG_QC_REF$JH_IID, euro_projection_JH_QC_REF$IID)]
merge_CAG_JH_QC_REF$PC9_CAG<-euro_projection_CAG_QC_REF$PC9[match(euro_projection_CAG_QC_REF$JH_IID, euro_projection_JH_QC_REF$IID)]
merge_CAG_JH_QC_REF$PC10_CAG<-euro_projection_CAG_QC_REF$PC10[match(euro_projection_CAG_QC_REF$JH_IID, euro_projection_JH_QC_REF$IID)]



list_PC<-c("PC1", "PC2" , "PC3", "PC4",  "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

for (PC in 1:length(list_PC)){
        print(list_PC[PC])
	
	
	png(paste0("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/comparison_",list_PC[PC],"CAG_vs_JH.png"), width=2300, height=1200)
	correlation_PC<-ggplot(merge_CAG_JH_QC_REF , aes(x=merge_CAG_JH_QC_REF[,PC+2], y=merge_CAG_JH_QC_REF[,PC+15])) +geom_point(alpha=0.3) + geom_abline(col="red") + xlab(paste0(list_PC[PC]," Impute-Merge Method"))+ ylab(paste0(list_PC[PC], " Merge-Impute Method"))+ theme(text = element_text(size = 20))
	print(correlation_PC)
	dev.off()
}




#-----------------correlation between genotyped PC and imputed PC----------------------------------
library(corrplot)


geno_QC_REF$PC1_imputed<-JH_QC_REF$PC1[match(JH_QC_REF$IID, geno_QC_REF$IID)]
geno_QC_REF$PC2_imputed<-JH_QC_REF$PC2[match(JH_QC_REF$IID, geno_QC_REF$IID)]
geno_QC_REF$PC3_imputed<-JH_QC_REF$PC3[match(JH_QC_REF$IID, geno_QC_REF$IID)]
geno_QC_REF$PC4_imputed<-JH_QC_REF$PC4[match(JH_QC_REF$IID, geno_QC_REF$IID)]
geno_QC_REF$PC5_imputed<-JH_QC_REF$PC5[match(JH_QC_REF$IID, geno_QC_REF$IID)]
geno_QC_REF$PC6_imputed<-JH_QC_REF$PC6[match(JH_QC_REF$IID, geno_QC_REF$IID)]
geno_QC_REF$PC7_imputed<-JH_QC_REF$PC7[match(JH_QC_REF$IID, geno_QC_REF$IID)]
geno_QC_REF$PC8_imputed<-JH_QC_REF$PC8[match(JH_QC_REF$IID, geno_QC_REF$IID)]
geno_QC_REF$PC9_imputed<-JH_QC_REF$PC9[match(JH_QC_REF$IID, geno_QC_REF$IID)]
geno_QC_REF$PC10_imputed<-JH_QC_REF$PC10[match(JH_QC_REF$IID, geno_QC_REF$IID)]



data_correlation<-data.frame( geno_QC_REF[,3:12], geno_QC_REF[,16:25])
#data_correlation<-data.frame(geno_QC_REF$IID, geno_QC_REF[,3:12], geno_QC_REF[,16:25])

matrix_cor<- cor(data_correlation)

png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/genotyped/PCA/correlation_PCA_imputed_vs_genotyped.png",  width=2300, height=1200)
p_cor<-corrplot(matrix_cor, type="upper", order="hclust", tl.col="black", tl.srt=45)

print(p_cor)
dev.off()



symnum(matrix_cor, abbr.colnames=FALSE)

#remove PC4, PC5, PC7 because they are problematic
data_without<-matrix_cor[-c(4,5,7, 14,15,17),-c(4,5,7, 14,15,17)]


#
data_cor_df<-data.frame(abs(matrix_cor))

#perfect_cor<-data_cor_df[which(data_cor_df$PC1 ==1 | data_cor_df$PC3 ==1 | data_cor_df$PC4 ==1| data_cor_df$PC5 ==1| data_cor_df$PC6 ==1| data_cor_df$PC7 ==1| data_cor_df$PC8 ==1| data_cor_df$PC9 ==1| data_cor_df$PC10 ==1| data_cor_df$PC1_imputed ==1| data_cor_df$PC2_imputed ==1 | data_cor_df$PC3_imputed ==1| data_cor_df$PC4_imputed ==1| data_cor_df$PC5_imputed ==1| data_cor_df$PC6_imputed ==1| data_cor_df$PC7_imputed ==1| data_cor_df$PC8_imputed ==1| data_cor_df$PC9_imputed ==1| data_cor_df$PC10_imputed ==1),]



#get the second most correlated
PC1<-data_cor_df[which(data_cor_df$PC1==max( data_cor_df$PC1[data_cor_df$PC1!=max(data_cor_df$PC1)] )),]
PC2<-data_cor_df[which(data_cor_df$PC2==max( data_cor_df$PC2[data_cor_df$PC2!=max(data_cor_df$PC2)] )),]
PC3<-data_cor_df[which(data_cor_df$PC3==max( data_cor_df$PC3[data_cor_df$PC3!=max(data_cor_df$PC3)] )),]
PC4<-data_cor_df[which(data_cor_df$PC4==max( data_cor_df$PC4[data_cor_df$PC4!=max(data_cor_df$PC4)] )),]
PC5<-data_cor_df[which(data_cor_df$PC5==max( data_cor_df$PC5[data_cor_df$PC5!=max(data_cor_df$PC5)] )),]
PC6<-data_cor_df[which(data_cor_df$PC6==max( data_cor_df$PC6[data_cor_df$PC6!=max(data_cor_df$PC6)] )),]
PC7<-data_cor_df[which(data_cor_df$PC7==max( data_cor_df$PC7[data_cor_df$PC7!=max(data_cor_df$PC7)] )),]
PC8<-data_cor_df[which(data_cor_df$PC8==max( data_cor_df$PC8[data_cor_df$PC8!=max(data_cor_df$PC8)] )),]
PC9<-data_cor_df[which(data_cor_df$PC9==max( data_cor_df$PC9[data_cor_df$PC9!=max(data_cor_df$PC9)] )),]
PC10<-data_cor_df[which(data_cor_df$PC10==max( data_cor_df$PC10[data_cor_df$PC10!=max(data_cor_df$PC10)] )),]


















#-----------------------ks test between arrays-------------------------------------------------------------------



list_chip_JH<-c("JH_gsa_760", "JH_gsa_4224", "JH_gsa_5300", "JH_gsa_17K")
list_chip_geno<-list_chip_JH

list_chip_CAG<-c("CAG_gsa_760", "CAG_gsa_4224", "CAG_gsa_5300", "CAG_gsa_17K")

#list_PC<-c("PC1", "PC2" , "PC3")
list_PC<-c("PC1", "PC2" , "PC3", "PC4",  "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")



#Genotyped


ks_geno_D<-c()
ks_geno_pval<-c()

chip1_geno<-c()
chip2_geno<-c()
list_geno_PC<-c()

for (PC in 1:length(list_PC)){
        print(list_PC[PC])
        countPC<-0

        for (chip1 in 1:length(list_chip_geno)){
                print(list_chip_geno[chip1])

                data_chip1<-geno_projection_CAG_QC_REF[which(geno_projection_CAG_QC_REF$Population==list_chip_geno[chip1]),PC+2]

                for (chip2 in 1:length(list_chip_geno)){
                        print(list_chip_geno[chip2])

                        data_chip2<-geno_projection_CAG_QC_REF[which(geno_projection_CAG_QC_REF$Population==(list_chip_geno[chip2])),PC+2]
			chip1_geno<-c(chip1_geno, list_chip_geno[chip1])
                        chip2_geno<-c(chip2_geno, list_chip_geno[chip2])

			#data_chip2<-geno_QC_REF[,PC+2]
                        #chip1_geno<-c(chip1_geno, list_chip_geno[chip1])
                        #chip2_geno<-c(chip2_geno, "QC REF")


                        head(data_chip2)

                        #ks test
                        test<-ks.test(data_chip1, data_chip2)


                        ks_geno_D<-c(ks_geno_D, as.numeric(test[1]))
                        ks_geno_pval<-c(ks_geno_pval, as.numeric(test[2]))
                        #list_geno_PC<-c(list_geno_PC, list_PC[PC])
                        countPC<-countPC+1
                }
        }

        list_geno_PC<-c(list_geno_PC, rep(c(list_PC[PC]),each=countPC))
}



data_ks_geno<-data.frame(list_geno_PC, chip1_geno, chip2_geno, ks_geno_D, ks_geno_pval)
data_ks_geno$significant<-"pval>0.05"
data_ks_geno[which(data_ks_geno$ks_geno_pval<0.05),]$significant<-"pval<=0.05"


data_ks_geno[which(data_ks_geno$chip1_geno=="JH_gsa_760"),]$chip1_geno<-"gsa_760"
data_ks_geno[which(data_ks_geno$chip1_geno=="JH_gsa_4224"),]$chip1_geno<-"gsa_4224"
data_ks_geno[which(data_ks_geno$chip1_geno=="JH_gsa_5300"),]$chip1_geno<-"gsa_5300"
data_ks_geno[which(data_ks_geno$chip1_geno=="JH_gsa_17K"),]$chip1_geno<-"gsa_17K"
data_ks_geno[which(data_ks_geno$chip2_geno=="JH_gsa_760"),]$chip2_geno<-"gsa_760"
data_ks_geno[which(data_ks_geno$chip2_geno=="JH_gsa_4224"),]$chip2_geno<-"gsa_4224"
data_ks_geno[which(data_ks_geno$chip2_geno=="JH_gsa_5300"),]$chip2_geno<-"gsa_5300"
data_ks_geno[which(data_ks_geno$chip2_geno=="JH_gsa_17K"),]$chip2_geno<-"gsa_17K"




data_ks_geno$list_geno_PC<-factor(data_ks_geno$list_geno_PC, levels=c("PC1", "PC2" , "PC3", "PC4",  "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))



png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/genotyped/PCA/QC_REF_geno_chip_ks_test_D.png", width=2300, height=1200)
#png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/genotyped/PCA/geno_chip_ks_test_D.png", width=2300, height=1200)
geno_plot<-ggplot(data_ks_geno , aes(x=chip1_geno, y=chip2_geno, fill=ks_geno_D, color=significant)) +geom_tile() +facet_wrap(~list_geno_PC) +xlab("Genotyping array")+ylab("Genotyping array")+ theme(text = element_text(size = 20))
print(geno_plot)
dev.off()
#ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/geno_chip_ks_test_D.png", geno_plot)



png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/genotyped/PCA/geno_chip_ks_test_pval.png", width=2300, height=1200)
geno_plot<-ggplot(data_ks_geno , aes(x=chip1_geno, y=chip2_geno, fill=ks_geno_pval)) +geom_tile() +facet_wrap(~list_geno_PC) +xlab("Genotyping array")+ylab("Genotyping array")+ theme(text = element_text(size = 20))
print(geno_plot)
dev.off()







#Imputed

ks_JH_D<-c()
ks_JH_pval<-c()

chip1_JH<-c()
chip2_JH<-c()
list_JH_PC<-c()

for (PC in 1:length(list_PC)){
	print(list_PC[PC])
	countPC<-0
	
	index_chip<-1
	for (chip1 in 1:length(list_chip_JH)){
		print(list_chip_JH[chip1])

                data_chip1<-euro_projection_JH_QC_REF[which(euro_projection_JH_QC_REF$Population==list_chip_JH[chip1]),PC+2]

		for (chip2 in index_chip:length(list_chip_JH)){
			print(list_chip_JH[chip2])
	
                        data_chip2<-euro_projection_JH_QC_REF[which(euro_projection_JH_QC_REF$Population==(list_chip_JH[chip2])),PC+2]
			chip1_JH<-c(chip1_JH, list_chip_JH[chip1])
			chip2_JH<-c(chip2_JH, list_chip_JH[chip2])
			
			#data_chip2<-JH_QC_REF[,PC+2]
                        #chip1_JH<-c(chip1_JH, list_chip_JH[chip1])
                        #chip2_JH<-c(chip2_JH, "QC REF")

			#head(data_chip2)

			#ks test 
			test<-ks.test(data_chip1, data_chip2)
			

			ks_JH_D<-c(ks_JH_D, as.numeric(test[1])) 
			ks_JH_pval<-c(ks_JH_pval, as.numeric(test[2])) 
			#list_JH_PC<-c(list_JH_PC, list_PC[PC])	
			countPC<-countPC+1
		}
		
		index_chip<-index_chip+1
	}
	
	list_JH_PC<-c(list_JH_PC, rep(c(list_PC[PC]),each=countPC))
}


data_ks_JH<-data.frame(list_JH_PC, chip1_JH, chip2_JH, ks_JH_D, ks_JH_pval)
data_ks_JH$significant<-"pval>0.05"
data_ks_JH[which(data_ks_JH$ks_JH_pval<0.05),]$significant<-"pval<=0.05"


data_ks_JH[which(data_ks_JH$chip1_JH=="JH_gsa_760"),]$chip1_JH<-"gsa_760"
data_ks_JH[which(data_ks_JH$chip1_JH=="JH_gsa_4224"),]$chip1_JH<-"gsa_4224"
data_ks_JH[which(data_ks_JH$chip1_JH=="JH_gsa_5300"),]$chip1_JH<-"gsa_5300"
data_ks_JH[which(data_ks_JH$chip1_JH=="JH_gsa_17K"),]$chip1_JH<-"gsa_17K"
data_ks_JH[which(data_ks_JH$chip2_JH=="JH_gsa_760"),]$chip2_JH<-"gsa_760"
data_ks_JH[which(data_ks_JH$chip2_JH=="JH_gsa_4224"),]$chip2_JH<-"gsa_4224"
data_ks_JH[which(data_ks_JH$chip2_JH=="JH_gsa_5300"),]$chip2_JH<-"gsa_5300"
data_ks_JH[which(data_ks_JH$chip2_JH=="JH_gsa_17K"),]$chip2_JH<-"gsa_17K"




data_ks_JH$list_JH_PC<-factor(data_ks_JH$list_JH_PC, levels=rev(c("PC1", "PC2" , "PC3", "PC4",  "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")))



png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/QC_REF_JH_chip_ks_test_D.png", width=2300, height=1200)
#png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/JH_chip_ks_test_D.png", width=2300, height=1200)
JH_plot<-ggplot(data_ks_JH , aes(x=chip1_JH, y=chip2_JH, fill=ks_JH_D, color=significant)) +geom_tile() +facet_wrap(~list_JH_PC) +xlab("Genotyping array")+ylab("Genotyping array")+ theme(text = element_text(size = 20))
print(JH_plot)
dev.off()
#ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/JH_chip_ks_test_D.png", JH_plot)



png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/JH_chip_ks_test_pval.png", width=2300, height=1200)
JH_plot<-ggplot(data_ks_JH , aes(x=chip1_JH, y=chip2_JH, fill=ks_JH_pval)) +geom_tile() +facet_wrap(~list_JH_PC) +xlab("Genotyping array")+ylab("Genotyping array")+ theme(text = element_text(size = 20))
print(JH_plot)
dev.off()







ks_CAG_D<-c()
ks_CAG_pval<-c()

chip1_CAG<-c()
chip2_CAG<-c()
list_CAG_PC<-c()


for (PC in 1:length(list_PC)){
        print(list_PC[PC])
	countPC<-0

        index_chip<-1

        for (chip1 in 1:length(list_chip_CAG)){
                print(list_chip_CAG[chip1])

                data_chip1<-euro_projection_CAG_QC_REF[which(euro_projection_CAG_QC_REF$Population==list_chip_CAG[chip1]),PC+2]

                for (chip2 in index_chip:length(list_chip_CAG)){
                        print(list_chip_CAG[chip2])

                        data_chip2<-euro_projection_CAG_QC_REF[which(euro_projection_CAG_QC_REF$Population==(list_chip_CAG[chip2])),PC+2]
                        chip1_CAG<-c(chip1_CAG, list_chip_CAG[chip1])
                        chip2_CAG<-c(chip2_CAG, list_chip_CAG[chip2])

			#data_chip2<-JH_QC_REF[,PC+2]
                        #chip1_CAG<-c(chip1_CAG, list_chip_CAG[chip1])
                        #chip2_CAG<-c(chip2_CAG, "QC REF")

                        #ks test
                        test<-ks.test(data_chip1, data_chip2)

                        ks_CAG_D<-c(ks_CAG_D, as.numeric(test[1]))
                        ks_CAG_pval<-c(ks_CAG_pval, as.numeric(test[2]))
                        #list_CAG_PC<-c(list_CAG_PC, list_PC[PC])
	                countPC<-countPC+1
		}
                index_chip<-index_chip+1

        }

	list_CAG_PC<-c(list_CAG_PC, rep(c(list_PC[PC]),each=countPC))

}

data_ks_CAG<-data.frame(list_CAG_PC, chip1_CAG, chip2_CAG, ks_CAG_D, ks_CAG_pval)
data_ks_CAG$significant<-"pval>0.05"
data_ks_CAG[which(data_ks_CAG$ks_CAG_pval<0.05),]$significant<-"pval<=0.05"

data_ks_CAG[which(data_ks_CAG$chip1_CAG=="CAG_gsa_760"),]$chip1_CAG<-"gsa_760"
data_ks_CAG[which(data_ks_CAG$chip1_CAG=="CAG_gsa_4224"),]$chip1_CAG<-"gsa_4224"
data_ks_CAG[which(data_ks_CAG$chip1_CAG=="CAG_gsa_5300"),]$chip1_CAG<-"gsa_5300"
data_ks_CAG[which(data_ks_CAG$chip1_CAG=="CAG_gsa_17K"),]$chip1_CAG<-"gsa_17K"
data_ks_CAG[which(data_ks_CAG$chip2_CAG=="CAG_gsa_760"),]$chip2_CAG<-"gsa_760"
data_ks_CAG[which(data_ks_CAG$chip2_CAG=="CAG_gsa_4224"),]$chip2_CAG<-"gsa_4224"
data_ks_CAG[which(data_ks_CAG$chip2_CAG=="CAG_gsa_5300"),]$chip2_CAG<-"gsa_5300"
data_ks_CAG[which(data_ks_CAG$chip2_CAG=="CAG_gsa_17K"),]$chip2_CAG<-"gsa_17K"




data_ks_CAG$list_CAG_PC<-factor(data_ks_CAG$list_CAG_PC, levels=rev(c("PC1", "PC2" , "PC3", "PC4",  "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")))



png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/QC_REF_CAG_chip_ks_test_D.png", width=2300, height=1200)
#png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/CAG_chip_ks_test_D.png", width=2300, height=1200)
CAG_plot<-ggplot(data_ks_CAG , aes(x=chip1_CAG, y=chip2_CAG, fill=ks_CAG_D, color=significant)) +geom_tile() +facet_wrap(~list_CAG_PC) +xlab("Genotyping array")+ylab("Genotyping array")+ theme(text = element_text(size = 20))   
print(CAG_plot)
dev.off()
#ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/CAG_chip_ks_test_D.png", CAG_plot)




png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/CAG_chip_ks_test_pval.png", width=2300, height=1200)
CAG_plot<-ggplot(data_ks_CAG , aes(x=chip1_CAG, y=chip2_CAG, fill=ks_CAG_pval)) +geom_tile() +facet_wrap(~list_CAG_PC) +xlab("Genotyping array")+ylab("Genotyping array")+ theme(text = element_text(size = 20))
print(CAG_plot)
dev.off()







#----------------------both in one plot-------------------------------

data_ks_JH$Method<-"Impute-Merge"
data_ks_CAG$Method<-"Merge-Impute"
data_ks_geno$Method<-"Genotyped"
data_ks_JH$D_CAG<-data_ks_CAG$ks_CAG_D
data_ks_JH$D_geno<-data_ks_geno$ks_geno_D



for (PC in 1:length(list_PC)){
        print(list_PC[PC])


	PC_data_ks_JH<-data_ks_JH[which(data_ks_JH$list_JH_PC==list_PC[PC]),]
	PC_data_ks_CAG<-data_ks_CAG[which(data_ks_CAG$list_CAG_PC==list_PC[PC]),]
        names(PC_data_ks_CAG)<-names(data_ks_JH)

	data_plot<-rbind(PC_data_ks_JH, PC_data_ks_CAG)		



	png(paste0("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/", list_PC[PC],"_chip_ks_test_D.png"), width=2300, height=1200)
	JH_plot<-ggplot(data_plot , aes(x=chip1_JH, y=chip2_JH, fill=ks_JH_D, color=significant)) +geom_tile() +facet_wrap(~Method) +xlab("Genotyping array")+ylab("Genotyping array")+ theme(text = element_text(size = 20))
	print(JH_plot)
	dev.off()
#ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/JH_chip_ks_test_D.png", JH_plot)

}


for (PC in 1:length(list_PC)){
        print(list_PC[PC])


        PC_data_ks_JH<-data_ks_JH[which(data_ks_JH$list_JH_PC==list_PC[PC]),]
PC_data_ks_JH$D_diff<-PC_data_ks_JH$D_CAG-PC_data_ks_JH$ks_JH_D

        png(paste0("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/", list_PC[PC],"_chip_ks_difference_D.png"), width=2300, height=1200)
        JH_plot<-ggplot(PC_data_ks_JH , aes(x=chip1_JH, y=chip2_JH, fill=D_diff)) +geom_tile() +xlab("Genotyping array")+ylab("Genotyping array")+ theme(text = element_text(size = 20))
        print(JH_plot)
        dev.off()



}



#-----compare D geno and JH + CAG----------------


tmp_CAG<-data_ks_JH
tmp_JH<-data_ks_JH

tmp_CAG$D<-tmp_CAG$D_CAG
tmp_JH$D<-tmp_JH$ks_JH_D

tmp_CAG$Method<-"Merge-Impute"
tmp_JH$Method<-"Impute-Merge"

geno_vs_imputed<-rbind(tmp_CAG, tmp_JH)
geno_vs_imputed$PC<-geno_vs_imputed$list_JH_PC





png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/genotyped/PCA/D_comparison_geno_JH_CAG_ks_test_D.png", width=2300, height=1200)
JH_plot<-ggplot(geno_vs_imputed , aes(x=D, y=D_geno)) +geom_point(aes(shape=PC, colour = Method), size=8) + geom_abline(col="red") + xlab("D in Imputed CaG")+ylab("D in Genotyped CaG")+ theme(text = element_text(size = 40))+ guides(fill=guide_legend(title="Comparison"))#+geom_smooth(method = "lm", se = FALSE)
#JH_plot<-ggplot(data_JH_D , aes(x=ks_JH_D, y=D_CAG)) +geom_point(aes(colour = factor(chip_comapre))) + geom_abline(col="red") + xlab("D in Impute-Merge Method")+ylab("D in Merge-Impute Method")+ theme(text = element_text(size = 20))+ guides(fill=guide_legend(title="Comparison")))
print(JH_plot)
dev.off()





png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/D_comparison_JH_CAG_ks_test_D.png", width=2300, height=1200)
JH_plot<-ggplot(geno_vs_imputed , aes(x=ks_JH_D, y=D_CAG)) +geom_point(aes(colour = PC), size=8) + geom_abline(col="red") + xlab("D in Impute-Merge")+ylab("D in Merge-Impute")+ theme(text = element_text(size = 40))+ theme(legend.title=element_blank())#+geom_smooth(method = "lm", se = FALSE)
#JH_plot<-ggplot(data_JH_D , aes(x=ks_JH_D, y=D_CAG)) +geom_point(aes(colour = factor(chip_comapre))) + geom_abline(col="red") + xlab("D in Impute-Merge Method")+ylab("D in Merge-Impute Method")+ theme(text = element_text(size = 20))+ guides(fill=guide_legend(title="Comparison")))
print(JH_plot)
dev.off()








cp_data_ks_CAG<-data_ks_CAG
names(cp_data_ks_CAG)<-names(data_ks_JH)[1:7]
tmp_data_final<-rbind(data_ks_JH[,1:7], cp_data_ks_CAG)
data_final_double<-tmp_data_final[which(tmp_data_final$chip1_JH!= tmp_data_final$chip2_JH),]
data_final_double$comparison<-paste0(data_final_double$chip1_JH, " - ", data_final_double$chip2_JH, " - " , data_final_double$list_JH_PC)
data_final_double$chip_comapre<-paste0(data_final_double$chip1_JH, " - ", data_final_double$chip2_JH)


data_JH_D<-data_final_double[which(data_final_double$Method=="Impute-Merge"),]
data_CAG_D<-data_final_double[which(data_final_double$Method=="Merge-Impute"),]

data_JH_D$D_CAG<-data_CAG_D$ks_JH_D[match(data_JH_D$comparison, data_CAG_D$comparison)]




data_JH_D$PC<-data_JH_D$list_JH_PC
data_JH_D$Comparison<-data_JH_D$chip_comapre



#mean point in distribution(0.04358611,  0.0442329) pente: 1.014839360521047



sizes <- c("PC1"=19,"PC2"=17,"PC3"=15,"PC4"=13,"PC5"=11,"PC6"=9,"PC7"=7,"PC8"=5,"PC9"=3,"PC10"=2)

png("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA_final/D_comparison_JH_CAG_ks_test_D.png", width=2300, height=1200)
#JH_plot<-ggplot(data_JH_D , aes(x=ks_JH_D, y=D_CAG)) +geom_point(aes(colour = factor(list_JH_PC))) + geom_abline(col="red") +geom_smooth(method = "lm", se = FALSE)+ xlab("D in Impute-Merge Method")+ylab("D in Merge-Impute Method")+ theme(text = element_text(size = 3330))+ guides(fill=guide_legend(title="Comparison"))
JH_plot<-ggplot(data_JH_D , aes(x=ks_JH_D, y=D_CAG)) +geom_point(aes(colour = Comparison, size=PC)) + geom_abline(col="red") + xlab("D in Impute-Merge Method")+ylab("D in Merge-Impute Method")+ theme(text = element_text(size = 40))+ guides(fill=guide_legend(title="Comparison"))+ scale_size_manual(values=sizes) #+scale_size(range = c(2,20)) #+ guides(colour = guide_legend(override.aes = list(size=10)))+ guides(color= guide_legend(),  size=guide_legend(override.aes = list(size = legend_size)))
print(JH_plot)
dev.off()





























#---------------------------1000G----------------------------------------#


pca_1000G<-G1000

#choose the data to plot
data_final_1000G<-pca_1000G[which(pca_1000G$cohort=="1000G"),]




pop_type<-unique(data_final_1000G$Population)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

data_final_1000G$Population<-factor(data_final_1000G$Population,levels=pop_type)



plot_1000G<-ggplot(data_final_1000G, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/1000G_PC1_PC2.png", plot_1000G)
p0<-ggplot(data_final_1000G, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25)) + facet_wrap(~Population) #+xlim(c(-0.1,0.55)) +ylim(c(-0.1,0.5))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/1000G_split_pop_PC1_PC2.png", p0)



for (i in 1:length(pop_type)){
	
	print(pop_type[i])
	#get the data for one pop
	data_final<- data_final_1000G[which(data_final_1000G$Population==pop_type[i]),]
	
	#plot the data
	plot<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population)) + geom_point(alpha=0.5)+scale_color_manual(values=col_vector[i])  +xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25))
	ggsave(paste0("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/1000G_", pop_type[i], "_PC1_PC2.png"), plot)

}






#-------------------------PROJECTIONS 1000G-------------------------------------------------------

#ADD CAG
data_final<-rbind(data_final_1000G, projection_CAG_1000G)


pop_type<-unique(data_final$Population)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

data_final$Population<-factor(data_final$Population,levels=pop_type)



plot_1000G<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CAG_1000G_PC1_PC2.png", plot_1000G)
p0<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25)) + facet_wrap(~Population) #+xlim(c(-0.1,0.55)) +ylim(c(-0.1,0.5))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CAG_1000G_split_pop_PC1_PC2.png", p0)


pop_type_CAG<-unique(data_final[which(data_final$cohort=="CAG"),]$Population)

for (i in 1:length(pop_type_CAG)){
	

        print(pop_type_CAG[i])
        #get the data for one pop
        data_final_tmp<-rbind(data_final_1000G, projection_CAG_1000G[which(projection_CAG_1000G$Population==pop_type_CAG[i]),])
	#data_final_tmp<-data_final[which(data_final$Population==pop_type_CAG[i] | data_final$cohort=="1000G"),]
	#data_final_tmp$Population<-factor(data_final$Population,levels=pop_type)
	col_vector_tmp<-col_vector
	col_vector_tmp[27]<-"black"

        #plot the data
        plot<-ggplot(data_final_tmp, aes(x=PC1, y=PC2, color=Population)) + geom_point(alpha=0.5)+scale_color_manual(values=col_vector_tmp)  +xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25))
        ggsave(paste0("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CAG_1000G_", pop_type_CAG[i], "_PC1_PC2.png"), plot)

}


#extract europeans from CAG

data_CAG_EUR<-data_final[which((data_final$PC1>=-0.12 & data_final$PC1<=-0.02) & (data_final$PC2>=0.15 & data_final$PC2<=0.25) & (data_final$cohort=="CAG")),]





#ADD CPTP
data_final<-rbind(data_final_1000G, projection_CPTP_1000G)

pop_type<-unique(data_final$Population)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

data_final$Population<-factor(data_final$Population,levels=pop_type)



plot_1000G<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CPTP_1000G_PC1_PC2.png", plot_1000G)
p0<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25)) + facet_wrap(~Population) #+xlim(c(-0.1,0.55)) +ylim(c(-0.1,0.5))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CPTP_1000G_split_pop_PC1_PC2.png", p0)


pop_type_CPTP<-unique(data_final[which(data_final$cohort=="CanPath"),]$Population)

for (i in 1:length(pop_type_CPTP)){


        print(pop_type_CPTP[i])
        #get the data for one pop
        data_final_tmp<-rbind(data_final_1000G, projection_CPTP_1000G[which(projection_CPTP_1000G$Population==pop_type_CPTP[i]),])
        #data_final_tmp<-data_final[which(data_final$Population==pop_type_CPTP[i] | data_final$cohort=="1000G"),]
        #data_final_tmp$Population<-factor(data_final$Population,levels=pop_type)
        col_vector_tmp<-col_vector
        col_vector_tmp[27]<-"black"

        #plot the data
        plot<-ggplot(data_final_tmp, aes(x=PC1, y=PC2, color=Population)) + geom_point(alpha=0.5)+scale_color_manual(values=col_vector_tmp)  +xlim(c(-0.205,0.35)) +ylim(c(-0.25,.25))
        ggsave(paste0("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CPTP_1000G_", pop_type_CPTP[i], "_PC1_PC2.png"), plot)

}






#---------------------------QC_REF----------------------------------------#


pca_QC_REF<-QC_REF


data_final_QC_REF<-pca_QC_REF[which(pca_QC_REF$cohort=="QC_REF"),]



pop_type<-unique(data_final_QC_REF$Population)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

data_final_QC_REF$Population<-factor(data_final_QC_REF$Population,levels=pop_type)


plot_QC_REF<-ggplot(data_final_QC_REF, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.7)+xlim(c(-0.2,0.2)) +ylim(c(-0.2,0.2))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/QC_REF_PC1_PC2.png", plot_QC_REF)
p0<-ggplot(data_final_QC_REF, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.7)+xlim(c(-0.2,0.2)) +ylim(c(-0.2,0.2)) + facet_wrap(~Population) #+xlim(c(-0.1,0.55)) +ylim(c(-0.1,0.5))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/QC_REF_split_pop_PC1_PC2.png", p0)



for (i in 1:length(pop_type)){

        print(pop_type[i])
        #get the data for one pop
        data_final<- data_final_QC_REF[which(data_final_QC_REF$Population==pop_type[i]),]

        #plot the data
        plot<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population)) + geom_point(alpha=0.7)+scale_color_manual(values=col_vector[i])  +xlim(c(-0.2,0.2)) +ylim(c(-0.2,0.2))
        ggsave(paste0("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/QC_REF_", pop_type[i], "_PC1_PC2.png"), plot)

}

#--------------------------Projections---------------------------------------------#

#CAG
data_final<-rbind(data_final_QC_REF, projection_CAG_QC_REF)


pop_type<-unique(data_final$Population)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

data_final$Population<-factor(data_final$Population,levels=pop_type)



plot_QC_REF<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.2,0.2)) +ylim(c(-0.2,0.2))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CAG_QC_REF_PC1_PC2.png", plot_QC_REF)
p0<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.2,0.2)) +ylim(c(-0.2,0.2)) + facet_wrap(~Population) #+xlim(c(-0.1,0.55)) +ylim(c(-0.1,0.5))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CAG_QC_REF_split_pop_PC1_PC2.png", p0)


pop_type_CAG<-unique(data_final[which(data_final$cohort=="CAG"),]$Population)

for (i in 1:length(pop_type_CAG)){


        print(pop_type_CAG[i])
        #get the data for one pop
        data_final_tmp<-rbind(data_final_QC_REF, projection_CAG_QC_REF[which(projection_CAG_QC_REF$Population==pop_type_CAG[i]),])
        #data_final_tmp<-data_final[which(data_final$Population==pop_type_CAG[i] | data_final$cohort=="QC_REF"),]
        #data_final_tmp$Population<-factor(data_final$Population,levels=pop_type)
        col_vector_tmp<-col_vector
        col_vector_tmp[12]<-"black"

        #plot the data
        plot<-ggplot(data_final_tmp, aes(x=PC1, y=PC2, color=Population)) + geom_point(alpha=0.5)+scale_color_manual(values=col_vector_tmp)  +xlim(c(-0.2,0.2)) +ylim(c(-0.2,0.2))
        ggsave(paste0("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CAG_QC_REF_", pop_type_CAG[i], "_PC1_PC2.png"), plot)

}




#CPTP
data_final<-rbind(data_final_QC_REF, projection_CPTP_QC_REF)


pop_type<-unique(data_final$Population)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

data_final$Population<-factor(data_final$Population,levels=pop_type)



plot_QC_REF<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.2,0.2)) +ylim(c(-0.2,0.2))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CPTP_QC_REF_PC1_PC2.png", plot_QC_REF)
p0<-ggplot(data_final, aes(x=PC1, y=PC2, color=Population))+scale_color_manual(values=col_vector) + geom_point(alpha=0.5)+xlim(c(-0.2,0.2)) +ylim(c(-0.2,0.2)) + facet_wrap(~Population) #+xlim(c(-0.1,0.55)) +ylim(c(-0.1,0.5))
ggsave("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CPTP_QC_REF_split_pop_PC1_PC2.png", p0)


pop_type_CPTP<-unique(data_final[which(data_final$cohort=="CanPath"),]$Population)

for (i in 1:length(pop_type_CPTP)){


        print(pop_type_CPTP[i])
        #get the data for one pop
        data_final_tmp<-rbind(data_final_QC_REF, projection_CPTP_QC_REF[which(projection_CPTP_QC_REF$Population==pop_type_CPTP[i]),])
        #data_final_tmp<-data_final[which(data_final$Population==pop_type_CPTP[i] | data_final$cohort=="QC_REF"),]
        #data_final_tmp$Population<-factor(data_final$Population,levels=pop_type)
        col_vector_tmp<-col_vector
        col_vector_tmp[12]<-"black"

        #plot the data
        plot<-ggplot(data_final_tmp, aes(x=PC1, y=PC2, color=Population)) + geom_point(alpha=0.5)+scale_color_manual(values=col_vector_tmp)  +xlim(c(-0.2,0.2)) +ylim(c(-0.2,0.2))
        ggsave(paste0("/lustre06/project/6065672/justinp/Common_positions_imputed_CAG_CPTP_QCREF_POPRES/PCA/CPTP_QC_REF_", pop_type_CPTP[i], "_PC1_PC2.png"), plot)

}



























