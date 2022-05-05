library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)


data<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/no_archi_modified_merged_raw_reimputed_WES.diff", header=F)
names(data)<-c("SNP", "FID", "IID", "NEW", "OLD")
#WITH ARCHI
#data<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/modified_merged_raw_reimputed_WES.diff", header=T)


well_reimputed<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/no_archi_well_reimputed_100percent_WES.ID", header=F)
names(well_reimputed)<-c("SNP")

#add ref and alt alleles
data$tmp_SNP<-data$SNP
data_sep<-separate(data, tmp_SNP, c(NA, NA, "ref", "alt"), ":")
data<-data_sep

#MAF and R2 from recalculated
probs_metric<-read.table("/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_from_CaG_TopMed/match_with_Exome_samples/noArchi/probsMetrics/all_chr_CAG_WES.probsMetrics" ,header=F)
names(probs_metric)<-c("SNP", "alt", "MAF", "a", "R2", "c", "d", "e")


#file containing SNP ID, AF , MAF ,R2 
data_final_info<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/common_reimputed_SNPs_CAG.info", header=F)
names(data_final_info)<-c("SNP", "AF_imputed", "MAF_imputed", "R2")

data_final_info$R2<-probs_metric$R2[match(data_final_info$SNP, probs_metric$SNP)]



freq_omni<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/genotyped_raw_CAG.frq", header=F)
names(freq_omni)<-c("SNP", "minor_allele", "maf_frq")

freq_192<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/reimputed_grch38_CAG_difference.frq", header=T)
names(freq_192)<-c("chr", "SNP", "alt_allele", "ref_allele", "maf_192", "nb_chrom")


freq_topmed<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/freq_reimputed_SNP_CAG.BRAVO_TOPMed_Freeze_8.frq", header=F)
names(freq_topmed)<-c("chr", "pos", "SNP", "ref", "alt", "topmed_AAF")

#change AF to MAF
freq_topmed_tmp <- transform(freq_topmed, MAF = pmin(topmed_AAF, (1-topmed_AAF)))
freq_topmed<-freq_topmed_tmp





#------------------------------------------------AAF to MAF and MAR comparaison from allelel--------------------------------------------------#


#freq_test<-freq_topmed

#freq_test$maf_frq_omni <- freq_omni$maf_frq[match(freq_test$SNP, freq_omni$SNP)]
#freq_test$maf_allele_omni <- freq_omni$minor_allele[match(freq_test$SNP, freq_omni$SNP)]


#see if the alternative allele is the same
#setequal(freq_test$alt, freq_test$maf_allele_omni)
#nrow(freq_test)
#nrow(freq_test[which(freq_test$maf_allele_omni==freq_test$alt),])

#There are some discordance with the alternative and the maf of the chip. Mainly those switch are due to the alternative alle frequency beeing higher in CAg than in TOPMED



#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#get the overall percentage of well reimputed sites


#Set the well imputed sites to 0 Freq
well_reimputed$Freq<-0
well_reimputed$perc<-0

#count the number of occurence for every SNPs
occurences_SNP<-as.data.frame(table(data$SNP))
occurences_SNP$perc<-((occurences_SNP$Freq)/90)*100
names(occurences_SNP)<-c("SNP", "Freq", "perc")

#get the baddes SNPs
#very_bad<-occurences_SNP[which(occurences_SNP$perc >=90),]

#merge two df to get full coverage
data_final<-rbind(occurences_SNP, well_reimputed)

#add freq values

#MAF and R2 from recalculated
data_final$maf <- freq_omni$maf_frq[match(data_final$SNP, freq_omni$SNP)]
data_final$maf_192 <- probs_metric$MAF[match(data_final$SNP, probs_metric$SNP)]


#add freq topmed values
#data_final$maf_topmed<- freq_topmed$MAF[match(data_final$SNP, freq_topmed$SNP)]
#same but for AAF
data_final$aaf_topmed<- freq_topmed$topmed_AAF[match(data_final$SNP, freq_topmed$SNP)]
data_final$maf_topmed<- freq_topmed$MAF[match(data_final$SNP, freq_topmed$SNP)]



#add maf_imputed and R2 values
data_final$maf_imputed<-data_final_info$MAF_imputed[match(data_final$SNP, data_final_info$SNP)]
data_final$R2<-data_final_info$R2[match(data_final$SNP, data_final_info$SNP)]



#-----------------------------------ICI-----------------------------#



SNPs_common<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/imputed_COMMON_Hussin_CAG_re_imputed_SNP.ID", header=F)
names(SNPs_common)<-c("SNP")


list_common_SNPs<-SNPs_common$SNP

data_common_imputed<-data_final[which(data_final$SNP %in% list_common_SNPs),]

data_final<-data_common_imputed







#get the big MAF difference
data_final$MAF_diff<-abs(data_final$maf-data_final$maf_topmed)

#get the big MAF difference
#data_final$MAF_diff<-abs(data_final$maf-data_final$aaf_topmed)
#test_diff<-data_final[which(data_final$MAF_diff>0.5),]
#data_tmp<-data_final[which(data_final$MAF_diff<=0.5),]
#data_final<-data_tmp
#write.table(test_diff , "omni_MAF_diff_over_0.5.Rdata", row.name=F)



#get the opposite of the bad percentage
data_final$goodPerc<-100-data_final$perc

#split percentage by grouo
data_final$groupPerc<-0
data_final[which(data_final$goodPerc>80&data_final$goodPerc<=90),]$groupPerc<-"80%<Accuracy<=90%"
data_final[which(data_final$goodPerc>90&data_final$goodPerc<=95),]$groupPerc<-"90%<Accuracy<=95%"
data_final[which(data_final$goodPerc>95&data_final$goodPerc<=99),]$groupPerc<-"95%<Accuracy<=99%"
data_final[which(data_final$goodPerc>99&data_final$goodPerc<100),]$groupPerc<-"99%<Accuracy<100%"
data_final[which(data_final$goodPerc==100),]$groupPerc<-"98%<Accuracy<=100%"



#order groups
data_final$groupPerc<-factor(data_final$groupPerc,levels=c("80%<Accuracy<=90%","90%<Accuracy<=95%","95%<Accuracy<=99%","99%<Accuracy<100%","98%<Accuracy<=100%"))



occurence_groupPerc<-as.data.frame(table(data_final$groupPerc))
names(occurence_groupPerc)<-c("Accuracy", "NumberOfSNPs")






#get the opposite of the bad percentage 
data_final$goodPerc<-100-data_final$perc

#split percentage by grouo
data_final$groupPerc<-0
#data_final[which(data_final$goodPerc<=1),]$groupPerc<-"Accuracy=0%"
#data_final[which(data_final$goodPerc>0&data_final$goodPerc<=1),]$groupPerc<-"0%<Accuracy<=1%"
#data_final[which(data_final$goodPerc>1&data_final$goodPerc<=5),]$groupPerc<-"1%<Accuracy<=5%"
#data_final[which(data_final$goodPerc>5&data_final$goodPerc<=10),]$groupPerc<-"5%<Accuracy<=10%"
#data_final[which(data_final$goodPerc>10&data_final$goodPerc<=20),]$groupPerc<-"10%<Accuracy<=20%"
#data_final[which(data_final$goodPerc>20&data_final$goodPerc<=30),]$groupPerc<-"20%<Accuracy<=30%"
#data_final[which(data_final$goodPerc>30&data_final$goodPerc<=40),]$groupPerc<-"30%<Accuracy<=40%"
#data_final[which(data_final$goodPerc>40&data_final$goodPerc<=50),]$groupPerc<-"40%<Accuracy<=50%"
#data_final[which(data_final$goodPerc>50&data_final$goodPerc<=60),]$groupPerc<-"50%<Accuracy<=60%"
#data_final[which(data_final$goodPerc>60&data_final$goodPerc<=70),]$groupPerc<-"60%<Accuracy<=70%"
#data_final[which(data_final$goodPerc>70&data_final$goodPerc<=80),]$groupPerc<-"70%<Accuracy<=80%"
#data_final[which(data_final$goodPerc>80&data_final$goodPerc<=90),]$groupPerc<-"80%<Accuracy<=90%"
#data_final[which(data_final$goodPerc>90&data_final$goodPerc<=95),]$groupPerc<-"90%<Accuracy<=95%"
#data_final[which(data_final$goodPerc>95&data_final$goodPerc<=99),]$groupPerc<-"95%<Accuracy<=99%"
#data_final[which(data_final$goodPerc>99&data_final$goodPerc<100),]$groupPerc<-"99%<Accuracy<100%"
#data_final[which(data_final$goodPerc==100),]$groupPerc<-"98%<Accuracy<=100%"
data_final[which(data_final$goodPerc<=80),]$groupPerc<-"Accuracy<=80%"
data_final[which(data_final$goodPerc>80&data_final$goodPerc<=95),]$groupPerc<-"80%<Accuracy<=95%"
data_final[which(data_final$goodPerc>95&data_final$goodPerc<=98),]$groupPerc<-"95%<Accuracy<=98%"
data_final[which(data_final$goodPerc>98&data_final$goodPerc<=100),]$groupPerc<-"98%<Accuracy<=100%"



#order groups
#data_final$groupPerc<-factor(data_final$groupPerc,levels=c("Accuracy=0%","0%<Accuracy<=1%","1%<Accuracy<=5%","5%<Accuracy<=10%","10%<Accuracy<=20%","20%<Accuracy<=30%","30%<Accuracy<=40%","40%<Accuracy<=50%","50%<Accuracy<=60%","60%<Accuracy<=70%","70%<Accuracy<=80%","80%<Accuracy<=90%","90%<Accuracy<=95%","95%<Accuracy<=99%","99%<Accuracy<100%","98%<Accuracy<=100%"))
data_final$groupPerc<-factor(data_final$groupPerc,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=98%", "98%<Accuracy<=100%"))




occurence_groupPerc<-as.data.frame(table(data_final$groupPerc))
names(occurence_groupPerc)<-c("Accuracy", "NumberOfSNPs")





#---------------------------------------PLOTS---------------------------------------------




#get the bad accuracy SNP lst
#wrong_SNP<-data_final[which(data_final$goodPerc<=30),]
#write.table(wrong_SNP, file="SNP_reimputed_accuracy_under_30perc.df", col.names=T, row.names=F)

#MAF distribution
#p0<-ggplot(data_final, aes(x=maf)) + geom_histogram(fill="blue")
p0<-ggplot(data_final, aes(x=maf_topmed)) + geom_histogram(fill="blue")
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnitest_reimputed_frequency_maf.png",p0)
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnitest_reimputed_frequency_topmed_maf.png",p0)

#plot the results
#p1<-ggplot(data_final, aes(x=SNP, y=perc)) + geom_point(alpha=0.3) + xlab("SNPs ID") + ylab("percentage of bad re-imputation (%)")
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnitest_reimputed_SNP_bad_accuracy.png",p1)


#density plot of the maf in the dataset
p1<-ggplot(data_final, aes(x=R2)) +  geom_density() + xlab("R2") + ylab("density") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_dist_R2.png",p1)




#by group of goodPerc
#png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_reimputed_SNP_accuracy.png", width = 4000, height = 2000)
#by group of goodPerc
png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_reimputed_SNP_accuracy.png", width = 2000, height = 2000)
p2<-ggplot(occurence_groupPerc, aes(x=Accuracy, y=NumberOfSNPs, fill=as.factor(Accuracy))) +  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) + geom_text(aes(label=NumberOfSNPs), position=position_dodge(width=0.9), size=15)+   theme(legend.position = "none", text = element_text(size = 50))  
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnireimputed_SNP_accuracy.png",p2)
print(p2)
dev.off()



#R2 from accuracy group
#p7<-ggplot(data_final, aes(x=groupPerc, y=R2, fill=as.factor(groupPerc))) + geom_violin() + xlab("Percentage of accuracy") + ylab("R2") +  theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnireimputed_SNP_accuracy_group_R2.png", width = 20, height = 10 ,p7)

#distribution of R2 in each group
p8<-ggplot(data_final, aes(x=R2)) + geom_density() +ylim(c(0,45)) + xlab("R2") + ylab("density") + facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnireimputed_SNP_distribution_R2.png", width = 20, height = 10 ,p8)

p13<-ggplot(data_final, aes(x=R2)) + geom_histogram()  +ylim(c(0,155000))+ xlab("R2") + ylab("Number of variants") + facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnireimputed_SNP_distribution_R2_hist.png", width = 20, height = 10 ,p13)
p13<-ggplot(data_final, aes(x=R2)) + geom_histogram()+ xlab("R2") + ylab("Number of variants") + facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/test_CAG_vs_WES_omnireimputed_SNP_distribution_R2_hist.png", width = 20, height = 10 ,p13)


#--------------------------------------------------compare omni maf vs TOPMED maf----------------------------------------------------



#MAF
#p9<-ggplot(data_final, aes(x=maf, y=maf_topmed, color=as.factor(colors))) + geom_point(alpha=0.3) + xlab("MAF omni") + ylab("MAF TOPMed")+facet_wrap(~data_final$groupPerc) + theme(legend.position = "none")
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnireimputed_SNP_maf_CAG_vs_TOPMed.png", width = 20, height = 10 ,p9)



#for only weird
#p9<-ggplot(data_test, aes(x=maf, y=maf_topmed, color=as.factor(colors))) + geom_point(alpha=0.3) + xlab("MAF omni") + ylab("MAF TOPMed") + theme(legend.position = "none")
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnicolored_reimputed_SNP_maf_CAG_vs_TOPMed.png", width = 20, height = 10 ,p9)





#----------------------------------------------------correlation between accuracy and R2-------------------------------------------------


#write data to file so we can change the SNP ID with awk
write.table(data_final, file="/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/data_final.df", col.names=F, row.names=F)

#bash commands;
#cut -d" " -f1 data_final.df | cut -d ":" -f 1,2 | sed 's/:/./g' | sed 's/"//g' |sed 's/chr//g' > col1_data_final.df
#paste col1_data_final.df data_final.df > modified_data_final.df


#read new modified file
modified_data_final<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/modified_data_final.df", header=F)
names(modified_data_final)<-c("Position", "SNP", "Freq", "perc", "maf", "maf_192", "aaf_topmed", "maf_topmed", "maf_imputed", "R2", "MAF_diff", "goodPerc", "groupPerc")




#groups for R2
modified_data_final$groupR2<-0
modified_data_final[which(modified_data_final$R2<0.3),]$groupR2<-"R2<0.3"
modified_data_final[which(modified_data_final$R2>=0.3 & modified_data_final$R2<0.8),]$groupR2<-"0.3<=R2<0.8"
modified_data_final[which(modified_data_final$R2>=0.8),]$groupR2<-"R2>=0.8"

modified_data_final$groupR2<-factor(modified_data_final$groupR2,levels=c("R2<0.3", "0.3<=R2<0.8", "R2>=0.8"))



#plot the accuracy and color by R2 values
p4<-ggplot(modified_data_final, aes(x=Position, y=goodPerc)) + geom_bin2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis") + xlab("Chromosome and positon") +ylab("Accuracy of imputation") + facet_wrap(~groupR2)
#p4<-ggplot(modified_data_final, aes(x=Position, y=goodPerc)) + geom_bin2d()+  scale_fill_gradient(low="yellow",high="darkblue",trans="log10") + xlab("Chromosome and positon") +ylab("Accuracy of Re-imputation") +  theme(legend.position = "none",) + facet_wrap(~groupR2)

#p4<-ggplot(modified_data_final, aes(x=Position, y=goodPerc)) + geom_point(alpha=0.5, aes(colour = factor(groupR2))) + xlab("Chromosome and positon") +ylab("Accuracy of Re-imputation") +  theme(legend.position = "none",) + facet_wrap(~groupR2)
#p4<-ggplot(modified_data_final, aes(x=Position, y=goodPerc)) + geom_point(alpha=0.5, aes(colour = factor(groupR2))) + xlab("Chromosome and positon") +ylab("Accuracy of Re-imputation") +  theme(legend.position = "none",) + facet_wrap(~groupR2)
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnireimputed_SNP_accuracy_and_R2.png",p4)




#AAF
#p9<-ggplot(modified_data_final, aes(x=maf, y=maf_topmed, color=as.factor(groupR2))) + geom_point(alpha=0.3) + xlab("MAF CAG") + ylab("MAF TOPMed") + facet_wrap(~modified_data_final$groupPerc) + guides(fill=guide_legend(title="Pre-imputation filter"))#+ scale_fill_discrete(name = "Pre-imputation filter")#+ theme(legend.position = "none")
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnireimputed_SNP_aaf_CAG_vs_TOPMed.png", width = 20, height = 10 ,p9)



modified_data_final_1<-modified_data_final[which(modified_data_final$groupR2=="R2>=0.8"),]
modified_data_final_2<-modified_data_final[which(modified_data_final$groupR2=="0.3<=R2<0.8"),]
modified_data_final_3<-modified_data_final[which(modified_data_final$groupR2=="R2<0.3"),]




p9_1<-ggplot(modified_data_final_1, aes(x=maf, y=maf_topmed, color=groupR2)) + geom_point(alpha=0.3) + scale_color_manual(values=c("#00BCF4")) + xlab("MAF CAG") + ylab("MAF TOPMed") + facet_wrap(~modified_data_final_1$groupPerc) + guides(fill=guide_legend(title="Pre-imputation filter"))#+ scale_fill_discrete(name = "Pre-imputation filter")#+ theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnireimputed_SNP_maf_CAG_vs_TOPMed_R2_08.png", width = 20, height = 10 ,p9_1)



p9_2<-ggplot(modified_data_final_2, aes(x=maf, y=maf_topmed, color=groupR2)) + geom_point(alpha=0.3) + scale_color_manual(values=c("#7CAE00")) + xlab("MAF CAG") + ylab("MAF TOPMed") + facet_wrap(~modified_data_final_2$groupPerc) + guides(fill=guide_legend(title="Pre-imputation filter"))#+ scale_fill_discrete(name = "Pre-imputation filter")#+ theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnireimputed_SNP_maf_CAG_vs_TOPMed_R2_03_08.png", width = 20, height = 10 ,p9_2)




p9_3<-ggplot(modified_data_final_3, aes(x=maf, y=maf_topmed, color=groupR2)) + geom_point(alpha=0.3) + scale_color_manual(values=c("#F8766D")) + xlab("MAF CAG") + ylab("MAF TOPMed") + facet_wrap(~modified_data_final_3$groupPerc) + guides(fill=guide_legend(title="Pre-imputation filter"))#+ scale_fill_discrete(name = "Pre-imputation filter")#+ theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_omnireimputed_SNP_maf_CAG_vs_TOPMed_R2_03.png", width = 20, height = 10 ,p9_3)






#-----------------------------------------------Comparaison with Hussin method--------------------------------------------------





#read new modified file for hussin method
modified_data_final_hussin<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/modified_data_final.df", header=F)
names(modified_data_final_hussin)<-c("Position", "SNP", "Freq", "perc", "maf","maf_192", "aaf_topmed", "maf_topmed", "maf_imputed", "R2", "MAF_diff", "goodPerc", "groupPerc")


#groups for R2
modified_data_final_hussin$groupR2<-0
modified_data_final_hussin[which(modified_data_final_hussin$R2<0.3),]$groupR2<-"R2<0.3"
modified_data_final_hussin[which(modified_data_final_hussin$R2>=0.3 & modified_data_final_hussin$R2<0.8),]$groupR2<-"0.3<=R2<0.8"
modified_data_final_hussin[which(modified_data_final_hussin$R2>=0.8),]$groupR2<-"R2>=0.8"

modified_data_final_hussin$groupR2<-factor(modified_data_final_hussin$groupR2,levels=c("R2<0.3", "0.3<=R2<0.8", "R2>=0.8"))

#change the name of CAG data
modified_data_final_CAG<-modified_data_final


#----Chi square---


modified_data_final_hussin$method="Impute-Merge"
modified_data_final_CAG$method="Merge-Impute"

modified_data_final_hussin$perfect_accuracy=NA
modified_data_final_CAG$perfect_accuracy=NA


modified_data_final_hussin[which(modified_data_final_hussin$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="oui"
modified_data_final_hussin[which(modified_data_final_hussin$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="non"

modified_data_final_CAG[which(modified_data_final_CAG$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="oui"
modified_data_final_CAG[which(modified_data_final_CAG$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="non"

merged_data<-rbind(modified_data_final_hussin, modified_data_final_CAG)

table(merged_data$method , merged_data$perfect_accuracy)

test<-chisq.test(merged_data$method , merged_data$perfect_accuracy, correct=FALSE)


#-------------------------Number of SNP in each R2 category Hussin--------------------


excluded_GWAS<-nrow(modified_data_final_hussin[which(modified_data_final_hussin$groupR2=="R2<0.3"),])/nrow(modified_data_final_hussin)
excluded_GWAS<-nrow(modified_data_final_hussin[which(modified_data_final_hussin$groupR2=="R2<0.3" & modified_data_final_hussin$goodPerc==100),])/nrow((modified_data_final_hussin[which(modified_data_final_hussin$groupR2=="R2<0.3"),]))


wrong_R2_good_ac<-modified_data_final_hussin[which(modified_data_final_hussin$groupR2=="R2<0.3" & modified_data_final_hussin$groupPerc=="98%<Accuracy<=100%"),]
perc_wrong_R2_good_ac<-(nrow(wrong_R2_good_ac)/nrow(modified_data_final_hussin))*100

wrong_R2_bad_ac<-modified_data_final_hussin[which((modified_data_final_hussin$groupR2=="0.3<=R2<0.8" | modified_data_final_hussin$groupR2=="R2>=0.8" ) & modified_data_final_hussin$groupPerc!="98%<Accuracy<=100%"),]
perc_wrong_R2_bad_ac<-(nrow(wrong_R2_bad_ac)/nrow(modified_data_final_hussin))*100



 excluded_GWAS<-nrow(modified_data_final_hussin[which(modified_data_final_hussin$groupR2=="R2<0.3" & modified_data_final_hussin$groupPerc=="98%<Accuracy<=100%"),])/nrow((modified_data_final_hussin[which(modified_data_final_hussin$goodPerc<98),]))


excluded_GWAS<-nrow(modified_data_final_hussin[which(modified_data_final_hussin$R2>=0.3),])/nrow(modified_data_final_hussin)

excluded_GWAS<-nrow(modified_data_final_hussin[which(modified_data_final_hussin$R2>=0.3 & modified_data_final_hussin$goodPerc<98),])/nrow((modified_data_final_hussin[which(modified_data_final_hussin$R2>=0.3),]))


#------------------------------------Frequency spectrum exons --------------------------------------------

nrow(modified_data_final_hussin[(modified_data_final_hussin$maf<0.01),])/nrow(modified_data_final_hussin)

SFS<-ggplot(modified_data_final_hussin, aes(x=maf)) + geom_histogram() +theme_classic() +xlab("MAF")+ylab("Number of variants")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/SFS_hussin_WES.png", SFS)




all_cag<-read.table("/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_all_arrays_09122021/maf0.01.hwe0.000001.miss5perc/Imputed/merging/R2_computation/noArchi/probsMetrics_files/filtered_all_chr_R2_0.0_Hussin.vcf.maf_r2", header=F)
names(all_cag)<-c("SNP", "MAF" ,"R2")

nrow(all_cag[which(all_cag$MAF<0.01),])/nrow(all_cag)

SFS<-ggplot(all_cag, aes(x=MAF)) + geom_histogram() +theme_classic() +xlab("MAF")+ylab("Number of variants")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/SFS_CAG_all_genome.png", SFS)








#-----------------------------check region badly imputed in JH--------------------------------------------

data_chr1_JH<-modified_data_final_hussin[which(as.numeric(modified_data_final_hussin$Position)<2),]
data_chr1_CAG<-modified_data_final_CAG[which(as.numeric(modified_data_final_CAG$Position)<2),]

modified_data_final_hussin$chr<-floor(modified_data_final_hussin$Position)
modified_data_final_CAG$chr<-floor(modified_data_final_CAG$Position)


p1<-ggplot(modified_data_final_hussin, aes(x=Position, y=R2)) + geom_smooth() +facet_wrap(~chr)
p2<-ggplot(modified_data_final_CAG, aes(x=Position, y=R2)) + geom_smooth() + facet_wrap(~chr)
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/JH_chr_1_R2_from_position.png", p1)
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/CAG_chr_1_R2_from_position.png", p2)


#----------------comparaison MAF CAG vs Hussin----------------------


#GET the WES real freq for the individuals
FREQ_WES<-read.table("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/renamed_SNPs_raw_genotyped_grch38_WES_difference.frq", header=T)
names(FREQ_WES)<-c("chr", "SNP", "alt_allele", "ref_allele", "maf_WES", "nb_chrom")


modified_data_final_CAG$maf_WES<-FREQ_WES$maf_WES[match(modified_data_final_CAG$SNP, FREQ_WES$SNP)]
modified_data_final_hussin$maf_WES<-FREQ_WES$maf_WES[match(modified_data_final_hussin$SNP, FREQ_WES$SNP)]




modified_data_final_CAG$maf_CAG<-modified_data_final_CAG$maf
modified_data_final_hussin$maf_Hussin<-modified_data_final_hussin$maf

modified_data_final_hussin$maf_CAG <- modified_data_final_CAG$maf_CAG[match(modified_data_final_hussin$SNP, modified_data_final_CAG$SNP)]



pMAF<-ggplot(modified_data_final_hussin, aes(x=maf_Hussin, y=maf_CAG)) + geom_point(alpha=0.3) + scale_color_manual(values=c("#F8766D")) + xlab("MAF Hussin") + ylab("MAF CAG") + facet_wrap(~modified_data_final_hussin$groupPerc) + guides(fill=guide_legend(title="Pre-imputation filter"))#+ scale_fill_discrete(name = "Pre-imputation filter")#+ theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_Hussin_maf_comparaison.png", width = 20, height = 10 ,pMAF)



pMAF_topmed_CAG<-ggplot(modified_data_final_CAG, aes(x=Position, y=MAF_diff)) + geom_point(alpha=0.3) + scale_color_manual(values=c("#F8766D")) + xlab("Position") + ylab("MAF_CAG - MAF TOPMED (abs)") + facet_wrap(~modified_data_final_CAG$groupPerc) + guides(fill=guide_legend(title="Pre-imputation filter"))
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_mafTOPMED_comparaison.png", width = 20, height = 10 ,pMAF_topmed_CAG)

pMAF_topmed_hussin<-ggplot(modified_data_final_hussin, aes(x=Position, y=MAF_diff)) + geom_point(alpha=0.3) + scale_color_manual(values=c("#F8766D")) + xlab("Position") + ylab("MAF_CAG - MAF TOPMED (abs)") + facet_wrap(~modified_data_final_hussin$groupPerc) + guides(fill=guide_legend(title="Pre-imputation filter"))
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_vs_WES_mafTOPMED_comparaison.png", width = 20, height = 10 ,pMAF_topmed_hussin)




#-----------------------------------------------ACCURACY comparaison HJ vs CAG
modified_data_final_CAG$maf_Hussin<-modified_data_final_hussin$maf_Hussin[match(modified_data_final_CAG$SNP, modified_data_final_hussin$SNP)]
modified_data_final_CAG$goodPerc_hussin<-modified_data_final_hussin$goodPerc[match(modified_data_final_CAG$SNP, modified_data_final_hussin$SNP)]



 modified_data_final_CAG$accuracy_com<-(as.numeric(modified_data_final_CAG$goodPerc)-as.numeric(modified_data_final_CAG$goodPerc_hussin))
 mean(modified_data_final_CAG$accuracy_com)

comp_accuracy<-ggplot(modified_data_final_CAG , aes(x=as.numeric(goodPerc), y=as.numeric(goodPerc_hussin))) + geom_bin_2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis")+geom_abline(col="red") +ylab("Accuracy Impute-Merge") + xlab("Accuracy Merge-Impute") 
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_Hussin_accuracy_comparaison.png",comp_accuracy)



#--------------------------------DIST R2---------------------------


p1<-ggplot(modified_data_final_CAG, aes(x=R2)) +  geom_density() + ylim(c(0,30)) + xlab("R2") + ylab("density") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_vs_WES_dist_R2.png",p1)


p1<-ggplot(modified_data_final_hussin, aes(x=R2)) +  geom_density() + ylim(c(0,30))+ xlab("R2") + ylab("density") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_vs_WES_dist_R2.png",p1)







#-------------------JH relation between R2 and Accuracy -------------------


p1<-ggplot(modified_data_final_hussin, aes(x=R2, y=goodPerc)) + geom_bin2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis") + xlab("R2") +ylab("Accuracy of imputation") 
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_R2_vs_accuracy.png",p1)



p1<-ggplot(modified_data_final_CAG, aes(x=R2, y=goodPerc)) + geom_bin2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis") + xlab("R2") +ylab("Accuracy of imputation") 
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/CAG_R2_vs_accuracy.png",p1)















#------------------- FINALLLL     JH and CAG accuracy groups ------------------
modified_data_final_hussin$groupPerc_CAG<-modified_data_final_CAG$groupPerc[match(modified_data_final_hussin$SNP, modified_data_final_CAG$SNP)]
modified_data_final_hussin$R2_CAG<-modified_data_final_CAG$R2[match(modified_data_final_hussin$SNP, modified_data_final_CAG$SNP)]




occurence_groupPerc_JH<-as.data.frame(table(modified_data_final_hussin$groupPerc))
names(occurence_groupPerc_JH)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_JH$method<-"Impute-Merge"


occurence_groupPerc_CAG<-as.data.frame(table(modified_data_final_hussin$groupPerc_CAG))
names(occurence_groupPerc_CAG)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_CAG$method<-"Merge-Impute"

occurence_groupPerc_JH_CAG<-rbind(occurence_groupPerc_JH, occurence_groupPerc_CAG)

occurence_groupPerc_JH_CAG$Accuracy<-factor(occurence_groupPerc_JH_CAG$Accuracy,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=98%", "98%<Accuracy<=100%"))
occurence_groupPerc_JH_CAG$method<-factor(occurence_groupPerc_JH_CAG$method, levels=c( "Impute-Merge", "Merge-Impute"))


png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/Hussin_and_CAG_reimputed_SNP_accuracy.png", width = 2000, height = 2000)
p2<-ggplot(occurence_groupPerc_JH_CAG, aes(x=Accuracy, y=NumberOfSNPs, fill=as.factor(method))) +  geom_bar(stat="identity",position = "dodge")+ylab("Number of imputed SNPs") + theme_classic()+ theme(axis.text.x = element_text(angle = 35, vjust = 0.9, hjust=1)) + theme(text = element_text(size = 70)) +theme(legend.position="none") #+ geom_text(aes(label=NumberOfSNPs), position=position_dodge(width=0.9), size=15)+ theme(text = element_text(size = 70)) + theme(legend.title=element_blank())#+ guides(fill=guide_legend(title="Method"))
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_vs_WES_omnireimputed_SNP_accuracy.png",p2)
print(p2)
dev.off()






tmp_JH<-data.frame(modified_data_final_hussin$SNP, modified_data_final_hussin$groupPerc, modified_data_final_hussin$R2)
tmp_JH$method<-"Impute-Merge"
names(tmp_JH)<-c("SNP", "groupPerc", "R2", "Data_set")

tmp_CAG<-data.frame(modified_data_final_hussin$SNP, modified_data_final_hussin$groupPerc_CAG, modified_data_final_hussin$R2_CAG)
tmp_CAG$method<-"Merge-Impute"
names(tmp_CAG)<-c("SNP", "groupPerc", "R2", "Data_set")

modified_data_final_JH_CAG<-rbind(tmp_JH, tmp_CAG)


modified_data_final_JH_CAG$groupPerc<-factor(modified_data_final_JH_CAG$groupPerc,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=98%", "98%<Accuracy<=100%"))


p13<-ggplot(modified_data_final_JH_CAG, aes(x=R2, fill=as.factor(Data_set))) + geom_histogram(position = "dodge")+ xlab("R2") + ylab("Number of variants") + facet_wrap(~modified_data_final_JH_CAG$groupPerc, scales = "free_y") + guides(fill=guide_legend(title="Data set")) 

ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/test_Hussin_and_CAG_SNP_distribution_R2_hist.png", width = 20, height = 10 ,p13)






#------------------------Accuracy and R2 histogram in JH method--------------------------




data_r2<-data.frame(modified_data_final_hussin$SNP, modified_data_final_hussin$maf, modified_data_final_hussin$maf_WES, modified_data_final_hussin$R2 )
names(data_r2)<-c("SNP", "MAF", "MAF_WES", "score" )
data_r2$quality_metric<-"R2"

data_accuracy<-data.frame(modified_data_final_hussin$SNP, modified_data_final_hussin$maf, modified_data_final_hussin$maf_WES, (modified_data_final_hussin$goodPerc/100))
names(data_accuracy)<-c("SNP", "MAF", "MAF_WES",  "score")
data_accuracy$quality_metric<-"Accuracy"

#final dataset
data_accuracy_r2<-rbind(data_r2, data_accuracy)


p1<-ggplot(data_accuracy_r2, aes(x=score, fill=quality_metric)) + geom_histogram(position="dodge") + ylab("Number of variants") + xlab("Quality of imputation") + guides(fill=guide_legend(title="Quality Metric")) 
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/R2_and_accuracy_dist_JH.png", p1)




#split by maf groups
data_accuracy_r2$group_maf<-0


data_accuracy_r2[which(data_accuracy_r2$MAF_WES<0.01),]$group_maf<-"MAF < 0.01"
data_maf<-data_accuracy_r2[which(data_accuracy_r2$MAF_WES<0.01),]
p1<-ggplot(data_maf, aes(x=score, fill=quality_metric)) + geom_histogram(position="dodge") + ylab("Number of variants") + xlab("Quality of imputation") + guides(fill=guide_legend(title="Quality Metric"))
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/maf_001_R2_and_accuracy_dist_JH.png", p1)


data_accuracy_r2[which(data_accuracy_r2$MAF_WES>=0.01 & data_accuracy_r2$MAF_WES<0.05),]$group_maf<-"0.01 <= MAF < 0.05"
data_maf<-data_accuracy_r2[which(data_accuracy_r2$MAF_WES>=0.01 & data_accuracy_r2$MAF_WES<0.05),]
p1<-ggplot(data_maf, aes(x=score, fill=quality_metric)) + geom_histogram(position="dodge") + ylab("Number of variants") + xlab("Quality of imputation") + guides(fill=guide_legend(title="Quality Metric"))
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/maf_001_005_R2_and_accuracy_dist_JH.png", p1)


data_accuracy_r2[which(data_accuracy_r2$MAF_WES>=0.05 & data_accuracy_r2$MAF_WES<0.1),]$group_maf<-"0.05 <= MAF < 0.1"
data_maf<-data_accuracy_r2[which(data_accuracy_r2$MAF_WES>=0.05 & data_accuracy_r2$MAF_WES<0.1),]
p1<-ggplot(data_maf, aes(x=score, fill=quality_metric)) + geom_histogram(position="dodge") + ylab("Number of variants") + xlab("Quality of imputation") + guides(fill=guide_legend(title="Quality Metric"))
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/maf_005_01_R2_and_accuracy_dist_JH.png", p1)


data_accuracy_r2[which(data_accuracy_r2$MAF_WES>=0.1 & data_accuracy_r2$MAF_WES<0.25),]$group_maf<-"0.1 <= MAF < 0.25"
data_maf<-data_accuracy_r2[which(data_accuracy_r2$MAF_WES>=0.1 & data_accuracy_r2$MAF_WES<0.25),]
p1<-ggplot(data_maf, aes(x=score, fill=quality_metric)) + geom_histogram(position="dodge") + ylab("Number of variants") + xlab("Quality of imputation") + guides(fill=guide_legend(title="Quality Metric"))
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/maf_01_025_R2_and_accuracy_dist_JH.png", p1)


data_accuracy_r2[which(data_accuracy_r2$MAF_WES>=0.25 & data_accuracy_r2$MAF_WES<=0.5),]$group_maf<-"0.25 <= MAF <= 0.5"
data_maf<-data_accuracy_r2[which(data_accuracy_r2$MAF_WES>=0.25 & data_accuracy_r2$MAF_WES<=0.5),]
p1<-ggplot(data_maf, aes(x=score, fill=quality_metric)) + geom_histogram(position="dodge") + ylab("Number of variants") + xlab("Quality of imputation") + guides(fill=guide_legend(title="Quality Metric"))
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/maf_025_05_R2_and_accuracy_dist_JH.png", p1)


cp_data_accuracy_r2<-data_accuracy_r2
cp_data_accuracy_r2$group_maf<-"All"

final_data_accuracy_r2<-rbind(cp_data_accuracy_r2, data_accuracy_r2)

final_data_accuracy_r2$group_maf<-factor(final_data_accuracy_r2$group_maf,levels=c("All","MAF < 0.01", "0.01 <= MAF < 0.05", "0.05 <= MAF < 0.1", "0.1 <= MAF < 0.25", "0.25 <= MAF <= 0.5"))


#facet_wrap
png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/maf__R2_and_accuracy_dist_JH.png", width = 1000, height = 500)
p1<-ggplot(final_data_accuracy_r2, aes(x=score, fill=quality_metric)) + geom_histogram(position="dodge") +theme_classic() + ylab("Number of variants") + xlab("Quality of imputation") + facet_wrap(~final_data_accuracy_r2$group_maf, scales = "free_y") +guides(fill=guide_legend(title="Quality Metric")) + theme(text = element_text(size = 20))      
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/maf__R2_and_accuracy_dist_JH.png", p1)
print(p1)
dev.off()

min(data_accuracy$score)
mean(data_accuracy$score)
mean(data_r2$score)
mean(data_accuracy$MAF_WES)
length(unique(data_accuracy_r2[which(data_accuracy_r2$group_maf=="MAF < 0.01"),]$SNP))



#-----------------------Percentage of SNP with bad R2 and high accuracy JH-------------------------------



ratio_bad_r2_perfect_acc<-(nrow(modified_data_final_hussin[which(modified_data_final_hussin$R2<0.3 & modified_data_final_hussin$goodPerc == 100),])/nrow(modified_data_final_hussin))

ration_bad_acc_good_r2<-(nrow(modified_data_final_hussin[which(modified_data_final_hussin$R2>=0.3 & modified_data_final_hussin$goodPerc<100),])/nrow(modified_data_final_hussin))

























#---------------------list position that have a maf over 0.1--------------------------------------------

maf_over_0_1_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_WES>0.1),]
maf_over_0_1_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_WES>0.1),]
maf_over_0_1_hussin$dataset<-"Impute-Merge"
maf_over_0_1_cag$dataset<-"Merge-Impute"

tmp_maf_over_0_1_hussin<-data.frame(maf_over_0_1_hussin$SNP, maf_over_0_1_hussin$R2, maf_over_0_1_hussin$dataset, maf_over_0_1_hussin$groupPerc)
names(tmp_maf_over_0_1_hussin)<-c("SNP","R2", "dataset", "groupPerc")
tmp_maf_over_0_1_cag<-data.frame(maf_over_0_1_cag$SNP, maf_over_0_1_cag$R2, maf_over_0_1_cag$dataset, maf_over_0_1_cag$groupPerc)
names(tmp_maf_over_0_1_cag)<-c("SNP","R2", "dataset", "groupPerc")

maf_over_0_1_merged<-rbind(tmp_maf_over_0_1_hussin, tmp_maf_over_0_1_cag)




p1<-ggplot(maf_over_0_1_merged, aes(x=R2, color=dataset)) +  geom_density() + xlab("R2") + ylab("density") + ggtitle("MAF over 0.1 in WES") #+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/MAF_over_01_inWES_CAG_and_JH_dist_R2.png",p1)

p1<-ggplot(maf_over_0_1_merged, aes(x=R2, color=dataset, fill=dataset)) +  geom_histogram(position=position_dodge(), alpha=0.5) + xlab("R2") + ggtitle("MAF over 0.1 in WES") #+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/MAF_over_01_inWES_CAG_and_JH_dist_R2_hist.png",p1)








#------density of R2 per ,method per MAF group--------------------------------------------------

#-------MAF<0.01
list_under_0_01_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_192<0.01),]$SNP
list_under_0_01_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_192<0.01),]$SNP
common<-intersect(list_under_0_01_hussin, list_under_0_01_cag)


maf_under_0_01_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$SNP %in% common),]
maf_under_0_01_cag<-modified_data_final_CAG[which(modified_data_final_CAG$SNP %in% common),]
maf_under_0_01_hussin$dataset<-"Impute-Merge"
maf_under_0_01_cag$dataset<-"Merge-Impute"

tmp_maf_under_0_01_hussin<-data.frame(maf_under_0_01_hussin$SNP, maf_under_0_01_hussin$R2, maf_under_0_01_hussin$dataset, maf_under_0_01_hussin$groupPerc)
names(tmp_maf_under_0_01_hussin)<-c("SNP","R2", "dataset", "groupPerc")
tmp_maf_under_0_01_cag<-data.frame(maf_under_0_01_cag$SNP, maf_under_0_01_cag$R2, maf_under_0_01_cag$dataset, maf_under_0_01_cag$groupPerc)
names(tmp_maf_under_0_01_cag)<-c("SNP","R2", "dataset", "groupPerc")

maf_under_0_01_merged<-rbind(tmp_maf_under_0_01_hussin, tmp_maf_under_0_01_cag)

maf_under_0_01_merged$maf_group<-"MAF < 0.01"


mean_value_Hussin<-mean(maf_under_0_01_merged[which(maf_under_0_01_merged$dataset =="Impute-Merge" ),]$R2)
vline_Hussin<-(mean_value_Hussin+0.02)
mean_value_CAG<-mean(maf_under_0_01_merged[which(maf_under_0_01_merged$dataset =="Merge-Impute" ),]$R2)
vline_CAG<-(mean_value_CAG-0.02)

p0<-ggplot(maf_under_0_01_merged, aes(x=R2, color=dataset)) +  geom_density()+ geom_vline(xintercept = mean_value_CAG, color = "#F8766D") + geom_text(aes(x=vline_CAG, label=as.character(round(mean_value_CAG, digits = 4)), y=55), colour="#F8766D", angle=90) + geom_vline(xintercept = mean_value_Hussin, color = "#00BFC4") + geom_text(aes(x=vline_Hussin, label=as.character(round(mean_value_Hussin, digits = 4)), y=55), colour="#00BFC4", angle=90) +ylim(c(0,60)) + xlab("R2") + ylab("density") + ggtitle(paste0("MAF under 0.01 in WES (", as.character(nrow(tmp_maf_under_0_01_hussin)), " SNPs)")) + guides(fill=guide_legend(title="Method"))#+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/MAF_under_0_01_inCAGandJH_dist_R2.png",p0)



maf_under_0_01_merged$groupPerc<-factor(maf_under_0_01_merged$groupPerc,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=99%", "99%<Accuracy<=100%"))


p13<-ggplot(maf_under_0_01_merged, aes(x=R2, fill=as.factor(dataset))) + geom_histogram(position = "dodge")+ xlab("R2")+ ggtitle(paste0("MAF under 0.01 in WES (", as.character(nrow(tmp_maf_under_0_01_hussin)), " SNPs)")) + ylab("Number of variants") + facet_wrap(~maf_under_0_01_merged$groupPerc, scales = "free_y") + guides(fill=guide_legend(title="Method"))

ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/maf_001_test_Hussin_and_CAG_SNP_distribution_R2_hist.png", width = 20, height = 10 ,p13)





occurence_groupPerc_JH<-as.data.frame(table(tmp_maf_under_0_01_hussin$groupPerc))
names(occurence_groupPerc_JH)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_JH$method<-"Impute-Merge"


occurence_groupPerc_CAG<-as.data.frame(table(tmp_maf_under_0_01_cag$groupPerc))
names(occurence_groupPerc_CAG)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_CAG$method<-"Merge-Impute"

occurence_groupPerc_JH_CAG<-rbind(occurence_groupPerc_JH, occurence_groupPerc_CAG)

occurence_groupPerc_JH_CAG$Accuracy<-factor(occurence_groupPerc_JH_CAG$Accuracy,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=99%", "99%<Accuracy<=100%"))
occurence_groupPerc_JH_CAG$method<-factor(occurence_groupPerc_JH_CAG$method, levels=c("Merge-Impute", "Impute-Merge"))


png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/maf_001_Hussin_and_CAG_reimputed_SNP_accuracy.png", width = 2000, height = 2000)
p2<-ggplot(occurence_groupPerc_JH_CAG, aes(x=Accuracy, y=NumberOfSNPs, fill=as.factor(method))) +  geom_bar(stat="identity",position = "dodge" ) + ggtitle(paste0("MAF under 0.01 in WES (", as.character(nrow(tmp_maf_under_0_01_hussin)), " SNPs)"))+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) + geom_text(aes(label=NumberOfSNPs), position=position_dodge(width=0.9), size=15)+ guides(fill=guide_legend(title="method"))+ theme(text = element_text(size = 50))
print(p2)
dev.off()




#-------0.01<=MAF<0.05

list_under_0_05_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_192>=0.01 & modified_data_final_hussin$maf_192<0.05),]$SNP
list_under_0_05_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_192>=0.01 & modified_data_final_CAG$maf_192<0.05),]$SNP
common<-intersect(list_under_0_05_hussin, list_under_0_05_cag)


maf_under_0_05_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$SNP %in% common),]
maf_under_0_05_cag<-modified_data_final_CAG[which(modified_data_final_CAG$SNP %in% common),]
maf_under_0_05_hussin$dataset<-"Impute-Merge"
maf_under_0_05_cag$dataset<-"Merge-Impute"

tmp_maf_under_0_05_hussin<-data.frame(maf_under_0_05_hussin$SNP, maf_under_0_05_hussin$R2, maf_under_0_05_hussin$dataset, maf_under_0_05_hussin$groupPerc)
names(tmp_maf_under_0_05_hussin)<-c("SNP","R2", "dataset", "groupPerc")
tmp_maf_under_0_05_cag<-data.frame(maf_under_0_05_cag$SNP, maf_under_0_05_cag$R2, maf_under_0_05_cag$dataset, maf_under_0_05_cag$groupPerc)
names(tmp_maf_under_0_05_cag)<-c("SNP","R2", "dataset", "groupPerc")

maf_under_0_05_merged<-rbind(tmp_maf_under_0_05_hussin, tmp_maf_under_0_05_cag)
maf_under_0_05_merged$maf_group<-" 0.01 <= MAF < 0.05"




mean_value_Hussin<-mean(maf_under_0_05_merged[which(maf_under_0_05_merged$dataset =="Impute-Merge" ),]$R2)
vline_Hussin<-(mean_value_Hussin+0.02)
mean_value_CAG<-mean(maf_under_0_05_merged[which(maf_under_0_05_merged$dataset =="Merge-Impute" ),]$R2)
vline_CAG<-(mean_value_CAG-0.02)


p1<-ggplot(maf_under_0_05_merged, aes(x=R2, color=dataset)) +  geom_density() + geom_vline(xintercept = mean_value_CAG, color = "#F8766D") + geom_text(aes(x=vline_CAG, label=as.character(round(mean_value_CAG, digits = 4)), y=55), colour="#F8766D", angle=90) + geom_vline(xintercept = mean_value_Hussin, color = "#00BFC4") + geom_text(aes(x=vline_Hussin, label=as.character(round(mean_value_Hussin, digits = 4)), y=55), colour="#00BFC4", angle=90) +ylim(c(0,60))  + xlab("R2") + ylab("density") + ggtitle(paste0("MAF greater or equal to 0.01 and under 0.05 in WES (", as.character(nrow(tmp_maf_under_0_05_hussin)), " SNPs)"))+ guides(fill=guide_legend(title="Method")) #+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/MAF_under_0_05_inCAGandJH_dist_R2.png",p1)


maf_under_0_05_merged$groupPerc<-factor(maf_under_0_05_merged$groupPerc,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=99%", "99%<Accuracy<=100%"))


p13<-ggplot(maf_under_0_05_merged, aes(x=R2, fill=as.factor(dataset))) + geom_histogram(position = "dodge")+ xlab("R2")+ ggtitle(paste0("MAF greater or equal to 0.01 and under 0.05 in WES (", as.character(nrow(tmp_maf_under_0_05_hussin)), " SNPs)"))+ ylab("Number of variants") + facet_wrap(~maf_under_0_05_merged$groupPerc, scales = "free_y") + guides(fill=guide_legend(title="Method"))

ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/maf_001_005_test_Hussin_and_CAG_SNP_distribution_R2_hist.png", width = 20, height = 10 ,p13)




occurence_groupPerc_JH<-as.data.frame(table(tmp_maf_under_0_05_hussin$groupPerc))
names(occurence_groupPerc_JH)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_JH$method<-"Impute-Merge"


occurence_groupPerc_CAG<-as.data.frame(table(tmp_maf_under_0_05_cag$groupPerc))
names(occurence_groupPerc_CAG)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_CAG$method<-"Merge-Impute"

occurence_groupPerc_JH_CAG<-rbind(occurence_groupPerc_JH, occurence_groupPerc_CAG)

occurence_groupPerc_JH_CAG$Accuracy<-factor(occurence_groupPerc_JH_CAG$Accuracy,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=99%", "99%<Accuracy<=100%"))
occurence_groupPerc_JH_CAG$method<-factor(occurence_groupPerc_JH_CAG$method, levels=c("Merge-Impute", "Impute-Merge"))


png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/maf_005_Hussin_and_CAG_reimputed_SNP_accuracy.png", width = 2000, height = 2000)
p2<-ggplot(occurence_groupPerc_JH_CAG, aes(x=Accuracy, y=NumberOfSNPs, fill=as.factor(method))) +  geom_bar(stat="identity",position = "dodge" ) + ggtitle(paste0("MAF greater or equal to 0.01 and under 0.05 in WES (", as.character(nrow(tmp_maf_under_0_05_hussin)), " SNPs)"))+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) + geom_text(aes(label=NumberOfSNPs), position=position_dodge(width=0.9), size=15)+ guides(fill=guide_legend(title="Method"))+ theme(text = element_text(size = 50))
print(p2)
dev.off()



#-------0.05<=MAF<0.1

list_under_0_1_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_192>=0.05 & modified_data_final_hussin$maf_192<0.1),]$SNP
list_under_0_1_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_192>=0.05 & modified_data_final_CAG$maf_192<0.1),]$SNP
common<-intersect(list_under_0_1_hussin, list_under_0_1_cag)


maf_under_0_1_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$SNP %in% common),]
maf_under_0_1_cag<-modified_data_final_CAG[which(modified_data_final_CAG$SNP %in% common),]
maf_under_0_1_hussin$dataset<-"Impute-Merge"
maf_under_0_1_cag$dataset<-"Merge-Impute"

tmp_maf_under_0_1_hussin<-data.frame(maf_under_0_1_hussin$SNP, maf_under_0_1_hussin$R2, maf_under_0_1_hussin$dataset, maf_under_0_1_hussin$groupPerc)
names(tmp_maf_under_0_1_hussin)<-c("SNP","R2", "dataset", "groupPerc")
tmp_maf_under_0_1_cag<-data.frame(maf_under_0_1_cag$SNP, maf_under_0_1_cag$R2, maf_under_0_1_cag$dataset, maf_under_0_1_cag$groupPerc)
names(tmp_maf_under_0_1_cag)<-c("SNP","R2", "dataset", "groupPerc")

maf_under_0_1_merged<-rbind(tmp_maf_under_0_1_hussin, tmp_maf_under_0_1_cag)
maf_under_0_1_merged$maf_group<-" 0.05 <= MAF < 0.1"




mean_value_Hussin<-mean(maf_under_0_1_merged[which(maf_under_0_1_merged$dataset =="Impute-Merge" ),]$R2)
vline_Hussin<-(mean_value_Hussin+0.02)
mean_value_CAG<-mean(maf_under_0_1_merged[which(maf_under_0_1_merged$dataset =="Merge-Impute" ),]$R2)
vline_CAG<-(mean_value_CAG-0.02)

p2<-ggplot(maf_under_0_1_merged, aes(x=R2, color=dataset)) +  geom_density()+ geom_vline(xintercept = mean_value_CAG, color = "#F8766D") + geom_text(aes(x=vline_CAG, label=as.character(round(mean_value_CAG, digits = 4)), y=55), colour="#F8766D", angle=90) + geom_vline(xintercept = mean_value_Hussin, color = "#00BFC4") + geom_text(aes(x=vline_Hussin, label=as.character(round(mean_value_Hussin, digits = 4)), y=55), colour="#00BFC4", angle=90) + ylim(c(0,60))  + xlab("R2") + ylab("density") + ggtitle(paste0("MAF greater or equal to 0.05 and under 0.1 in WES (", as.character(nrow(tmp_maf_under_0_1_hussin)), " SNPs)"))+ guides(fill=guide_legend(title="Method")) #+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/MAF_under_0_1_inCAGandJH_dist_R2.png",p2)


maf_under_0_1_merged$groupPerc<-factor(maf_under_0_1_merged$groupPerc,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=99%", "99%<Accuracy<=100%"))


p13<-ggplot(maf_under_0_1_merged, aes(x=R2, fill=as.factor(dataset))) + geom_histogram(position = "dodge")+ xlab("R2")+ ggtitle(paste0("MAF greater or equal to 0.05 and under 0.1 in WES (", as.character(nrow(tmp_maf_under_0_1_hussin)), " SNPs)")) + ylab("Number of variants") + facet_wrap(~maf_under_0_1_merged$groupPerc, scales = "free_y") + guides(fill=guide_legend(title="Data set"))

ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/maf_01_test_Hussin_and_CAG_SNP_distribution_R2_hist.png", width = 20, height = 10 ,p13)




occurence_groupPerc_JH<-as.data.frame(table(tmp_maf_under_0_1_hussin$groupPerc))
names(occurence_groupPerc_JH)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_JH$method<-"Impute-Merge"


occurence_groupPerc_CAG<-as.data.frame(table(tmp_maf_under_0_1_cag$groupPerc))
names(occurence_groupPerc_CAG)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_CAG$method<-"Merge-Impute"

occurence_groupPerc_JH_CAG<-rbind(occurence_groupPerc_JH, occurence_groupPerc_CAG)

occurence_groupPerc_JH_CAG$Accuracy<-factor(occurence_groupPerc_JH_CAG$Accuracy,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=99%", "99%<Accuracy<=100%"))
occurence_groupPerc_JH_CAG$method<-factor(occurence_groupPerc_JH_CAG$method, levels=c("Merge-Impute", "Impute-Merge"))


png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/maf_01_Hussin_and_CAG_reimputed_SNP_accuracy.png", width = 2000, height = 2000)
p2<-ggplot(occurence_groupPerc_JH_CAG, aes(x=Accuracy, y=NumberOfSNPs, fill=as.factor(method))) +  geom_bar(stat="identity",position = "dodge" ) + ggtitle(paste0("MAF greater or equal to 0.05 and under 0.1 in WES (", as.character(nrow(tmp_maf_under_0_1_hussin)), " SNPs)"))+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) + geom_text(aes(label=NumberOfSNPs), position=position_dodge(width=0.9), size=15)+ guides(fill=guide_legend(title="Data set"))+ theme(text = element_text(size = 50))
print(p2)
dev.off()


#-------0.1<=MAF<0.25

list_under_0_25_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_192>=0.1 & modified_data_final_hussin$maf_192<0.25),]$SNP
list_under_0_25_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_192>=0.1 & modified_data_final_CAG$maf_192<0.25),]$SNP
common<-intersect(list_under_0_25_hussin, list_under_0_25_cag)


maf_under_0_25_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$SNP %in% common),]
maf_under_0_25_cag<-modified_data_final_CAG[which(modified_data_final_CAG$SNP %in% common),]
maf_under_0_25_hussin$dataset<-"Impute-Merge"
maf_under_0_25_cag$dataset<-"Merge-Impute"

tmp_maf_under_0_25_hussin<-data.frame(maf_under_0_25_hussin$SNP, maf_under_0_25_hussin$R2, maf_under_0_25_hussin$dataset, maf_under_0_25_hussin$groupPerc)
names(tmp_maf_under_0_25_hussin)<-c("SNP","R2", "dataset", "groupPerc")
tmp_maf_under_0_25_cag<-data.frame(maf_under_0_25_cag$SNP, maf_under_0_25_cag$R2, maf_under_0_25_cag$dataset, maf_under_0_25_cag$groupPerc)
names(tmp_maf_under_0_25_cag)<-c("SNP","R2", "dataset", "groupPerc")

maf_under_0_25_merged<-rbind(tmp_maf_under_0_25_hussin, tmp_maf_under_0_25_cag)
maf_under_0_25_merged$maf_group<-" 0.1 <= MAF < 0.25"



mean_value_Hussin<-mean(maf_under_0_25_merged[which(maf_under_0_25_merged$dataset =="Impute-Merge" ),]$R2)
vline_Hussin<-(mean_value_Hussin+0.02)
mean_value_CAG<-mean(maf_under_0_25_merged[which(maf_under_0_25_merged$dataset =="Merge-Impute" ),]$R2)
vline_CAG<-(mean_value_CAG-0.02)

p3<-ggplot(maf_under_0_25_merged, aes(x=R2, color=dataset)) +  geom_density()+ geom_vline(xintercept = mean_value_CAG, color = "#F8766D") + geom_text(aes(x=vline_CAG, label=as.character(round(mean_value_CAG, digits = 4)), y=55), colour="#F8766D", angle=90) + geom_vline(xintercept = mean_value_Hussin, color = "#00BFC4") + geom_text(aes(x=vline_Hussin, label=as.character(round(mean_value_Hussin, digits = 4)), y=55), colour="#00BFC4", angle=90)+ylim(c(0,60))+ xlab("R2") + ylab("density") + ggtitle(paste0("MAF greater or equal to 0.1 and under 0.25 in WES (", as.character(nrow(tmp_maf_under_0_25_hussin)), " SNPs)"))+ guides(fill=guide_legend(title="Method")) #+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/MAF_under_0_25_inCAGandJH_dist_R2.png",p3)


maf_under_0_25_merged$groupPerc<-factor(maf_under_0_25_merged$groupPerc,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=99%", "99%<Accuracy<=100%"))

p13<-ggplot(maf_under_0_25_merged, aes(x=R2, fill=as.factor(dataset))) + geom_histogram(position = "dodge")+ xlab("R2")+ ggtitle(paste0("MAF greater or equal to 0.1 and under 0.25 in WES (", as.character(nrow(tmp_maf_under_0_25_hussin)), " SNPs)")) + ylab("Number of variants") + facet_wrap(~maf_under_0_25_merged$groupPerc, scales = "free_y") + guides(fill=guide_legend(title="Data set"))

ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/maf_025_test_Hussin_and_CAG_SNP_distribution_R2_hist.png", width = 20, height = 10 ,p13)



occurence_groupPerc_JH<-as.data.frame(table(tmp_maf_under_0_25_hussin$groupPerc))
names(occurence_groupPerc_JH)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_JH$method<-"Impute-Merge"


occurence_groupPerc_CAG<-as.data.frame(table(tmp_maf_under_0_25_cag$groupPerc))
names(occurence_groupPerc_CAG)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_CAG$method<-"Merge-Impute"

occurence_groupPerc_JH_CAG<-rbind(occurence_groupPerc_JH, occurence_groupPerc_CAG)

occurence_groupPerc_JH_CAG$Accuracy<-factor(occurence_groupPerc_JH_CAG$Accuracy,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=99%", "99%<Accuracy<=100%"))
occurence_groupPerc_JH_CAG$method<-factor(occurence_groupPerc_JH_CAG$method, levels=c("Merge-Impute", "Impute-Merge"))


png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/maf_025_Hussin_and_CAG_reimputed_SNP_accuracy.png", width = 2000, height = 2000)
p2<-ggplot(occurence_groupPerc_JH_CAG, aes(x=Accuracy, y=NumberOfSNPs, fill=as.factor(method))) +  geom_bar(stat="identity",position = "dodge" )  + ggtitle(paste0("MAF greater or equal to 0.1 and under 0.25 in WES (", as.character(nrow(tmp_maf_under_0_25_hussin)), " SNPs)"))+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) + geom_text(aes(label=NumberOfSNPs), position=position_dodge(width=0.9), size=15)+ guides(fill=guide_legend(title="Data set"))+ theme(text = element_text(size = 50))
print(p2)
dev.off()


#-------0.25<=MAF<0.5

list_under_0_5_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_192>=0.25 & modified_data_final_hussin$maf_192<=0.5),]$SNP
list_under_0_5_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_192>=0.25 & modified_data_final_CAG$maf_192<=0.5),]$SNP
common<-intersect(list_under_0_5_hussin, list_under_0_5_cag)


maf_under_0_5_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$SNP %in% common),]
maf_under_0_5_cag<-modified_data_final_CAG[which(modified_data_final_CAG$SNP %in% common),]
maf_under_0_5_hussin$dataset<-"Impute-Merge"
maf_under_0_5_cag$dataset<-"Merge-Impute"

tmp_maf_under_0_5_hussin<-data.frame(maf_under_0_5_hussin$SNP, maf_under_0_5_hussin$R2, maf_under_0_5_hussin$dataset, maf_under_0_5_hussin$groupPerc)
names(tmp_maf_under_0_5_hussin)<-c("SNP","R2", "dataset", "groupPerc")
tmp_maf_under_0_5_cag<-data.frame(maf_under_0_5_cag$SNP, maf_under_0_5_cag$R2, maf_under_0_5_cag$dataset, maf_under_0_5_cag$groupPerc)
names(tmp_maf_under_0_5_cag)<-c("SNP","R2", "dataset", "groupPerc")

maf_under_0_5_merged<-rbind(tmp_maf_under_0_5_hussin, tmp_maf_under_0_5_cag)
maf_under_0_5_merged$maf_group<-" 0.25 <= MAF <= 0.5"



mean_value_Hussin<-mean(maf_under_0_5_merged[which(maf_under_0_5_merged$dataset =="Impute-Merge" ),]$R2)
vline_Hussin<-(mean_value_Hussin+0.02)
mean_value_CAG<-mean(maf_under_0_5_merged[which(maf_under_0_5_merged$dataset =="Merge-Impute" ),]$R2)
vline_CAG<-(mean_value_CAG-0.02)

p4<-ggplot(maf_under_0_5_merged, aes(x=R2, color=dataset)) +  geom_density()+ geom_vline(xintercept = mean_value_CAG, color = "#F8766D") + geom_text(aes(x=vline_CAG, label=as.character(round(mean_value_CAG, digits = 4)), y=55), colour="#F8766D", angle=90) + geom_vline(xintercept = mean_value_Hussin, color = "#00BFC4") + geom_text(aes(x=vline_Hussin, label=as.character(round(mean_value_Hussin, digits = 4)), y=55), colour="#00BFC4", angle=90)+ylim(c(0,60))  + xlab("R2") + ylab("density") + ggtitle(paste0("MAF greater or equal to 0.25 and equal or under 0.5 in WES (", as.character(nrow(tmp_maf_under_0_5_hussin)), " SNPs)")) + guides(fill=guide_legend(title="Method"))#+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/MAF_under_0_5_inCAGandJH_dist_R2.png",p4)



maf_under_0_5_merged$groupPerc<-factor(maf_under_0_5_merged$groupPerc,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=99%", "99%<Accuracy<=100%"))


p13<-ggplot(maf_under_0_5_merged, aes(x=R2, fill=as.factor(dataset))) + geom_histogram(position = "dodge")+ xlab("R2") + ggtitle(paste0("MAF greater or equal to 0.25 and under 0.5 in WES (", as.character(nrow(tmp_maf_under_0_5_hussin)), " SNPs)")) + ylab("Number of variants") + facet_wrap(~maf_under_0_5_merged$groupPerc, scales = "free_y") + guides(fill=guide_legend(title="Data set"))

ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/maf_05_test_Hussin_and_CAG_SNP_distribution_R2_hist.png", width = 20, height = 10 ,p13)



occurence_groupPerc_JH<-as.data.frame(table(tmp_maf_under_0_5_hussin$groupPerc))
names(occurence_groupPerc_JH)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_JH$method<-"Impute-Merge"


occurence_groupPerc_CAG<-as.data.frame(table(tmp_maf_under_0_5_cag$groupPerc))
names(occurence_groupPerc_CAG)<-c("Accuracy", "NumberOfSNPs")
occurence_groupPerc_CAG$method<-"Merge-Impute"

occurence_groupPerc_JH_CAG<-rbind(occurence_groupPerc_JH, occurence_groupPerc_CAG)

occurence_groupPerc_JH_CAG$Accuracy<-factor(occurence_groupPerc_JH_CAG$Accuracy,levels=c("Accuracy<=80%","80%<Accuracy<=95%", "95%<Accuracy<=99%", "99%<Accuracy<=100%"))
occurence_groupPerc_JH_CAG$method<-factor(occurence_groupPerc_JH_CAG$method, levels=c("Merge-Impute", "Impute-Merge"))


png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/maf_05_Hussin_and_CAG_reimputed_SNP_accuracy.png", width = 2000, height = 2000)
p2<-ggplot(occurence_groupPerc_JH_CAG, aes(x=Accuracy, y=NumberOfSNPs, fill=as.factor(method))) +  geom_bar(stat="identity",position = "dodge" )  + ggtitle(paste0("MAF greater or equal to 0.25 and under 0.5 in WES (", as.character(nrow(tmp_maf_under_0_5_hussin)), " SNPs)"))+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) + geom_text(aes(label=NumberOfSNPs), position=position_dodge(width=0.9), size=15)+ guides(fill=guide_legend(title="Data set"))+ theme(text = element_text(size = 50))
print(p2)
dev.off()





if(FALSE){
list_over_0_1_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_192>0.1),]$SNP
list_over_0_1_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_192>0.1),]$SNP
common<-intersect(list_over_0_1_hussin, list_over_0_1_cag)


maf_over_0_1_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$SNP %in% common),]
maf_over_0_1_cag<-modified_data_final_CAG[which(modified_data_final_CAG$SNP %in% common),]
maf_over_0_1_hussin$dataset<-"Impute-Merge"
maf_over_0_1_cag$dataset<-"Merge-Impute"

tmp_maf_over_0_1_hussin<-data.frame(maf_over_0_1_hussin$SNP, maf_over_0_1_hussin$R2, maf_over_0_1_hussin$dataset)
names(tmp_maf_over_0_1_hussin)<-c("SNP","R2", "dataset")
tmp_maf_over_0_1_cag<-data.frame(maf_over_0_1_cag$SNP, maf_over_0_1_cag$R2, maf_over_0_1_cag$dataset)
names(tmp_maf_over_0_1_cag)<-c("SNP","R2", "dataset")

maf_over_0_1_merged<-rbind(tmp_maf_over_0_1_hussin, tmp_maf_over_0_1_cag)



mean_value<-mean(maf_over_0_1_merged$R2)

p5<-ggplot(maf_over_0_1_merged, aes(x=R2, color=dataset)) +  geom_density()+ geom_vline(xintercept = mean_value, color = "black")+geom_text(aes(x=0.95, label=as.character(round(mean_value, digits = 4)), y=40), colour="black", angle=90) + xlab("R2") + ylab("density") + ggtitle("MAF over 0.1 in WES") #+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/MAF_over_01_inCAGandJH_dist_R2.png",p5)
#p1<-ggplot(maf_over_0_1_merged, aes(x=R2, color=dataset, fill=dataset)) +  geom_histogram(position=position_dodge(), alpha=0.5) + xlab("R2") + ggtitle("MAF over 0.1 in WES") #+ theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/MAF_over_01_inCAGandJH_dist_R2_hist.png",p1)
}






#---------------all in one----------------


maf_under_hussin<-modified_data_final_hussin
maf_under_cag<-modified_data_final_CAG
maf_under_hussin$dataset<-"Impute-Merge"
maf_under_cag$dataset<-"Merge-Impute"

tmp_maf_under_hussin<-data.frame(maf_under_hussin$SNP, maf_under_hussin$R2, maf_under_hussin$dataset, maf_under_hussin$groupPerc)
names(tmp_maf_under_hussin)<-c("SNP","R2", "dataset", "groupPerc")
tmp_maf_under_cag<-data.frame(maf_under_cag$SNP, maf_under_cag$R2, maf_under_cag$dataset, maf_under_cag$groupPerc)
names(tmp_maf_under_cag)<-c("SNP","R2", "dataset", "groupPerc")

maf_under_merged<-rbind(tmp_maf_under_hussin, tmp_maf_under_cag)
maf_under_merged$maf_group<-"All"



final_maf_under<-rbind(maf_under_merged, maf_under_0_01_merged, maf_under_0_05_merged, maf_under_0_1_merged, maf_under_0_25_merged, maf_under_0_5_merged)

final_maf_under$maf_group<-factor(final_maf_under$maf_group, levels=c("All", "MAF < 0.01" , " 0.01 <= MAF < 0.05", " 0.05 <= MAF < 0.1", " 0.1 <= MAF < 0.25", " 0.25 <= MAF <= 0.5"))


#means for each facet
CAG_final_maf_under<-final_maf_under[which(final_maf_under$dataset=="Merge-Impute"),]
JH_final_maf_under<-final_maf_under[which(final_maf_under$dataset=="Impute-Merge"),]
CAG_mean <- CAG_final_maf_under %>%   group_by(maf_group) %>%   summarize(mean_val = mean(R2))
JH_mean <- JH_final_maf_under %>%   group_by(maf_group) %>%   summarize(mean_val = mean(R2))





png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_JH_CAG/MAF_all_groups_Hussin_and_CAG_reimputed_SNP_accuracy.png", width = 1000, height = 600)
p2<-ggplot(final_maf_under,  aes(x=R2, fill=dataset)) +  geom_histogram(position = "dodge")+theme_classic()+ facet_wrap(~maf_group, scales = "free_y")+ geom_vline(data= CAG_mean, aes(xintercept=mean_val, color = "#F8766D")) + geom_vline(data= JH_mean, aes(xintercept=mean_val, color = "#00BFC4")) + xlab("R2") + ylab("Number of imputed SNPs") +guides(color = "none")+ guides(fill=guide_legend(title="Method"))+ theme(text = element_text(size = 20)) +theme(legend.title=element_blank())    

print(p2)
dev.off()










figure1 <- ggarrange(p0, p1, labels = c("A", "B"), ncol = 1, nrow = 2)
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/1_combined_MAF_under_inCAGandJH_dist_R2.png",figure1)

figure2 <- ggarrange(p2, p3, p4, labels = c( "C", "D", "E"), ncol = 1, nrow = 2)
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/2_combined_MAF_under_inCAGandJH_dist_R2.png",figure2)



p0_big<- p0 + theme(text = element_text(size = 20))  
p1_big<- p1 + theme(text = element_text(size = 20))  
p2_big<- p2 + theme(text = element_text(size = 20))  
p3_big<- p3 + theme(text = element_text(size = 20))  
p4_big<- p4 + theme(text = element_text(size = 20))  


png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/combined_MAF_under_inCAGandJH_dist_R2.png", width = 1000, height = 2500)
figure1 <- ggarrange(p0_big, p1_big, p2_big, p3_big, p4_big, labels = c("A", "B", "C", "D", "E"), ncol = 1, nrow = 5)
print(figure1)
dev.off()





















#-----------------------------------------investigate the low R2 in Hussin method------------------------

modified_data_final_hussin<-read.table("/lustre07/scratch/justinp/topmed_new/imputation_pipeline/WES_reimputation/CAG_MERGED/modified_data_final_hussin.Rdata", header=T)



hussin_low_r2<-modified_data_final_hussin[which(modified_data_final_hussin$R2<=0.05),]

p1<-ggplot(hussin_low_r2, aes(x=R2)) +  geom_density() + ylim(c(0,30))+ xlab("R2") + ylab("density") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_vs_WES_dist_R2_under0_05.png",p1)




p1<-ggplot(hussin_low_r2, aes(x=Position)) +  geom_histogram(color="black") + xlim(c(0,23))+ xlab("Position") + ylab("density") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_vs_WES_position_under0_05.png",p1)



cag_low_r2<-modified_data_final_CAG[which(modified_data_final_CAG$SNP %in% hussin_low_r2$SNP),]


#get their accuracy
#compare R2 with the same in CAG
hussin_low_r2$R2_CAG<-cag_low_r2$R2[match(hussin_low_r2$SNP, cag_low_r2$SNP)]
hussin_low_r2$goodPerc_CAG<-cag_low_r2$goodPerc[match(hussin_low_r2$SNP, cag_low_r2$SNP)]
hussin_low_r2$maf_192_cag<-cag_low_r2$maf_192[match(hussin_low_r2$SNP, cag_low_r2$SNP)]
hussin_low_r2$groupPerc_cag<-cag_low_r2$groupPerc[match(hussin_low_r2$SNP, cag_low_r2$SNP)]
hussin_low_r2$groupR2_cag<-cag_low_r2$groupR2[match(hussin_low_r2$SNP, cag_low_r2$SNP)]
hussin_low_r2$MAF_WES_minus_MAF_192<-abs(hussin_low_r2$maf_WES - hussin_low_r2$maf_192)
hussin_low_r2$MAF_topmed_minus_MAF_192<-abs(hussin_low_r2$maf_topmed - hussin_low_r2$maf_192)

#opposite, ones that are low in CAG R2_CAG<=0.05
hussin_low_r2_cag<-modified_data_final_hussin[which(modified_data_final_hussin$R2_CAG<=0.05),]

p1<-ggplot(hussin_low_r2_cag, aes(x=R2)) +  geom_histogram() + xlab("R2") + ylab("density") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/Hussin_vs_WES_dist_R2_under0_05_in_CAG.png",p1)


#split per chr
modified_data_final_hussin$chromosome<-floor(modified_data_final_hussin$Position)







#----------------DOSAGES-----------------
dosage_CAG<-read.table("/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_from_CaG_TopMed/match_with_Exome_samples/low_R2_JH_inCAG.matchExomesSamples.vcf.dosages",header=F)
id_dosage_CAG<-read.table("/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_from_CaG_TopMed/match_with_Exome_samples/ID.txt", header=F)
t<-c("chr", "pos", "SNP")
n<-c(t,id_dosage_CAG$V1)
names(dosage_CAG)<-n



dosage_JH<-read.table("/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_all_arrays_09122021/maf0.01.hwe0.000001.miss5perc/Imputed/match_with_Exome_samples/low_R2_JH_inJH.matchExomesSamples.vcf.dosages",header=F)
id_dosage_JH<-read.table("/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_all_arrays_09122021/maf0.01.hwe0.000001.miss5perc/Imputed/match_with_Exome_samples/ID.txt", header=F)
t<-c("chr", "pos", "SNP")
n<-c(t,id_dosage_JH$V1)
names(dosage_JH)<-n


dosage_CAG$all_0<-"not"
dosage_CAG$var_dosage<-0
dosage_CAG$mean_dosage<-0
for (row in 1:nrow((dosage_CAG))){

	occurence_0<-(as.numeric(dosage_CAG[row,4:195]))
        if(var(occurence_0) == 0 ){
		dosage_CAG[row,]$all_0<-"yes"
	}

	dosage_CAG[row,]$var_dosage<-var(as.numeric(dosage_CAG[row,4:195]))
	dosage_CAG[row,]$mean_dosage<-mean(as.numeric(dosage_CAG[row,4:195]))
}

dosage_JH$all_0<-"not"
dosage_JH$var_dosage<-0
dosage_JH$mean_dosage<-0
for (row in 1:nrow(dosage_JH)){
        
	occurence_0<-(as.numeric(dosage_JH[row,4:195]))
        if(var(occurence_0) == 0){
                dosage_JH[row,]$all_0<-"yes"
        }

	dosage_JH[row,]$var_dosage<-var(as.numeric(dosage_JH[row,4:195]))
        dosage_JH[row,]$mean_dosage<-mean(as.numeric(dosage_JH[row,4:195]))
}


#add dosage to all SNPs
modified_data_final_hussin$dosage_CAG<-dosage_CAG$mean_dosage[match(modified_data_final_hussin$SNP, dosage_CAG$SNP)]
modified_data_final_hussin$dosage_JH<-dosage_JH$mean_dosage[match(modified_data_final_hussin$SNP, dosage_JH$SNP)]
modified_data_final_hussin$var_dosage_CAG<-dosage_CAG$var_dosage[match(modified_data_final_hussin$SNP, dosage_CAG$SNP)]
modified_data_final_hussin$var_dosage_JH<-dosage_JH$var_dosage[match(modified_data_final_hussin$SNP, dosage_JH$SNP)]
modified_data_final_hussin$all_0_CAG<-dosage_CAG$all_0[match(modified_data_final_hussin$SNP, dosage_CAG$SNP)]
modified_data_final_hussin$all_0_Hussin<-dosage_JH$all_0[match(modified_data_final_hussin$SNP, dosage_JH$SNP)]




#write.table(modified_data_final_hussin, "/lustre07/scratch/justinp/topmed_new/imputation_pipeline/WES_reimputation/CAG_MERGED/v2_modified_data_final_hussin.Rdata", row.name=F, col.name=T)

modified_data_final_hussin<-read.table("/lustre07/scratch/justinp/topmed_new/imputation_pipeline/WES_reimputation/CAG_MERGED/v2_modified_data_final_hussin.Rdata", header=T)
#modified_data_final_hussin<-read.table("/lustre07/scratch/justinp/topmed_new/imputation_pipeline/WES_reimputation/CAG_MERGED/modified_data_final_hussin.Rdata", header=T)







#subset that are bad in both datasets
hussin_low_r2_cag<-modified_data_final_hussin[which(modified_data_final_hussin$R2_CAG<0.3 & modified_data_final_hussin$R2<=0.05),]
#18948 SNPs vs 

hussin_low_r2_cag$pos_on_chr<-0

for (row in 1:nrow(hussin_low_r2_cag)){

        position<-str_split(hussin_low_r2_cag[row,]$SNP, ":")[[1]][2]
        #print(position)

        hussin_low_r2_cag[row,]$pos_on_chr<-position

}


#LOW in JH and HIGH in CAG
hussin_low_r2_high_cag<-modified_data_final_hussin[which(modified_data_final_hussin$R2_CAG>=0.3 & modified_data_final_hussin$R2<=0.05),]
#36621 SNPs

hussin_low_r2_high_cag$pos_on_chr<-0

for (row in 1:nrow(hussin_low_r2_high_cag)){
	
        position<-str_split(hussin_low_r2_high_cag[row,]$SNP, ":")[[1]][2]
        #print(position)

	hussin_low_r2_high_cag[row,]$pos_on_chr<-position

}

#Dosage per individual per SNP------------------

#duplicate every SNPs (rows) for the number of individuals
list_ind_JH<-(names(dosage_JH))[4:195]
list_ind_CAG<-(names(dosage_CAG))[4:195]

dosage_final<-data.frame(NA, NA, NA, NA)
names(dosage_final)<-c("SNP", "ind", "dosage_JH", "dosage_CAG")
#loop through the SNPs
for (row in 1:nrow(dosage_JH)){
	#loop through the individuals
	for (ind in 1:length(list_ind_JH)){
		
		#print(list_ind_JH[ind])
                
		#get the good column JH
		ind_id<-list_ind_JH[ind]
		column_ind_JH<-match(c(ind_id), names(dosage_JH))
		column_ind_CAG<-match(c(ind_id), names(dosage_CAG))


		#get the SNP ID
		snp_id<-dosage_JH[row,]$SNP
		#get the JH dosage for this SNP

		dos_JH<-dosage_JH[which(dosage_JH$SNP==snp_id),column_ind_JH]
		dos_CAG<-dosage_CAG[which(dosage_CAG$SNP==snp_id),column_ind_CAG]

		row_to_add<-data.frame(snp_id,ind_id, dos_JH, dos_CAG)
		names(row_to_add)<-c("SNP", "ind", "dosage_JH", "dosage_CAG")
		#print(row_to_add)
		
		#expand the initial dataset
		dosage_final<-rbind(dosage_final, row_to_add)
	}


}
#remove empty first row
dosage_final <- dosage_final[-1,]



write.table(dosage_final, "/lustre07/scratch/justinp/topmed_new/imputation_pipeline/WES_reimputation/CAG_MERGED/Dosage_per_ind_per_SNP.Rdata", row.name=F, col.name=T)

dosage_final<-read.table("/lustre07/scratch/justinp/topmed_new/imputation_pipeline/WES_reimputation/CAG_MERGED/Dosage_per_ind_per_SNP.Rdata", header=T)



#--------------------extract one example for dosages-------------------------------------

#REORDER COLUMS AND SAMPLES ID DOSAGE
order_dosage_CAG<-names(dosage_JH)
ordered_dosage_CAG <- dosage_CAG[, order_dosage_CAG]
dosage_CAG<-ordered_dosage_CAG

#bad in both
bad_both<-hussin_low_r2_cag[which(hussin_low_r2_cag$R2_CAG<0.05 & hussin_low_r2_cag$R2<0.05),]
bad_both_SNP_id<-bad_both[2,2]
#first SNP id "chr1:100188744:A:G"
bad_both_CAG<-dosage_CAG[which(dosage_CAG$SNP==bad_both_SNP_id),]
bad_both_JH<-dosage_JH[which(dosage_JH$SNP==bad_both_SNP_id),]

write.table(bad_both_CAG, "/lustre07/scratch/justinp/topmed_new/imputation_pipeline/WES_reimputation/CAG_MERGED/bad_both_CAG.dosage", row.names=F)
write.table(bad_both_JH, "/lustre07/scratch/justinp/topmed_new/imputation_pipeline/WES_reimputation/CAG_MERGED/bad_both_JH.dosage", row.names=F)




#god in CAG and bad in J
bad_JH<-hussin_low_r2_high_cag[which(hussin_low_r2_high_cag$R2_CAG>0.9 & hussin_low_r2_high_cag$R2<0.05),]
bad_JH_SNP_id<-bad_JH[2,2]
#first SNP id "chr1:100188744:A:G"
bad_JH_CAG<-dosage_CAG[which(dosage_CAG$SNP==bad_JH_SNP_id),]
bad_JH_JH<-dosage_JH[which(dosage_JH$SNP==bad_JH_SNP_id),]


write.table(bad_JH_CAG, "/lustre07/scratch/justinp/topmed_new/imputation_pipeline/WES_reimputation/CAG_MERGED/bad_JH_CAG.dosage", row.names=F)
write.table(bad_JH_JH, "/lustre07/scratch/justinp/topmed_new/imputation_pipeline/WES_reimputation/CAG_MERGED/bad_JH_JH.dosage", row.names=F)








#--------------for each SNP number SNP with dosage bewteen 0 and 1 excluded

SNP_id<-c()
number_JH<-c()
number_CAG<-c()
for (row in 1:nrow(dosage_JH)){
	
	#print(row)
	row_JH<-as.numeric(dosage_JH[row,4:195])
	row_CAG<-as.numeric(dosage_CAG[row,4:195])
	number_JH<-c(number_JH,length(row_JH[which(row_JH!=0 & row_JH!=1 & row_JH!=2)]))
	number_CAG<-c(number_CAG,length(row_CAG[which(row_CAG!=0 & row_CAG!=1 & row_CAG!=2)]))
	#number_JH<-c(number_JH,length(row_JH[which(row_JH>0 & row_JH<1)]))
	#number_CAG<-c(number_CAG,length(row_CAG[which(row_CAG>0 & row_CAG<1)]))

	SNP_id<-c(SNP_id, dosage_JH[row,]$SNP)
	
}

data_number_dosage<-data.frame(SNP_id, number_JH, number_CAG)
data_number_dosage$diff_number<-data_number_dosage$number_JH - data_number_dosage$number_CAG
data_number_dosage$index<-seq(1, nrow(data_number_dosage))

data_number<-data_number_dosage[which(data_number_dosage$SNP_id %in% hussin_low_r2_high_cag$SNP),]
data_number$index<-seq(1, nrow(data_number))

#data_number<-data_number_dosage[which(data_number_dosage$SNP_id %in% hussin_low_r2$SNP),]



p1<-ggplot(data_number, aes(x=number_CAG, y=number_JH)) + geom_point(alpha=0.2) + ggtitle(paste0("SNPs with R2<0.05 in JH and R2>=0.3 in CAG ",nrow(data_number), " SNPs")) + xlab("Number of sample with 0<dosage<1 CAG") + ylab("Number of sample with 0<dosage<1 JH")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/number_of_not_0_1_2_Dosages_low_R2_JH_high_CAG.png", p1)
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/number_of_between_0_and_1_Dosages_low_R2_JH_high_CAG.png", p1)


p1<-ggplot(data_number, aes(x=index,y=diff_number)) + geom_point(alpha=0.2)+ geom_hline(yintercept=mean(data_number$diff_number), color="red") + geom_text(aes(x=36000, label=as.character((mean(data_number$diff_number)), digits = 4), y=(mean(data_number$diff_number)+2)), colour="red") + ggtitle(paste0("SNPs with R2<0.05 in JH and R2>=0.3 in CAG ",nrow(data_number), " SNPs")) + xlab("Index of SNP") + ylab("Number of sample with 0<dosage<1 JH-CAG")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/DIFF_number_of_note_0_1_2_Dosages_low_R2_JH_high_CAG.png", p1)
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/DIFF_number_of_between_0_and_1_Dosages_low_R2_JH_high_CAG.png", p1)





#dosage not mean--------------------------

#get the subset of low R2 in JH and high in CAG
dosage_JH_tmp<-dosage_JH[which(dosage_JH$SNP %in% hussin_low_r2_high_cag$SNP),]
dosage_CAG_tmp<-dosage_CAG[which(dosage_CAG$SNP %in% hussin_low_r2_high_cag$SNP),]



dosage_JH_vector<-c()
dosage_CAG_vector<-c()
for (col in 4:(ncol(dosage_CAG))){
#for (col in 4:(ncol(dosage_CAG)-3)){
        print(col)

        tmp<-dosage_JH_tmp[,col]
        dosage_JH_vector<-c(dosage_JH_vector, tmp)

        tmp<-dosage_CAG_tmp[,col]
        dosage_CAG_vector<-c(dosage_CAG_vector, tmp)

}

data_dosage<-data.frame(dosage_JH_vector, dosage_CAG_vector)


#data_dosage$dosage_diff<-abs(dosage_JH_vector - dosage_CAG_vector)
#data_dosage$group_dosage_diff<-0
#data_dosage[which(data_dosage$dosage_diff<0.01),]$group_dosage_diff<-"diff<0.01"
#data_dosage[which(data_dosage$dosage_diff>=0.01 & data_dosage$dosage_diff<0.02),]$group_dosage_diff<-"0.1<=diff<0.02"
#data_dosage[which(data_dosage$dosage_diff>0.02),]$group_dosage_diff<-"diff>0.02"
#data_dosage$group_dosage_diff<-factor(data_dosage$group_dosage_diff,levels=c("diff<0.01", "diff>0.02","0.1<=diff<0.02"))


#plot the distributions of the the dosages
p1<-ggplot(data_dosage, aes(x=dosage_CAG_vector, y=dosage_JH_vector)) + geom_point(alpha=0.2) +ggtitle(paste0("SNPs with R2<0.05 in JH and R2>=0.3 in CAG ",nrow(hussin_low_r2_high_cag), " SNPs")) + xlab("Dosage CAG") + ylab("Dosage JH")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/not_mean_Dosages_low_R2_JH_high_CAG.png", p1)




#write.table(hussin_low_r2_high_cag, "hussin_low_r2_high_cag.Rdata", row.names=F, col.names=T)
#hussin_low_r2_high_cag<-read.table("hussin_low_r2_high_cag.Rdata", header=T)

#---------------------dosages final plots-----------------------
hussin_low_r2_high_cag$dosage_diff<-abs(hussin_low_r2_high_cag$dosage_JH -hussin_low_r2_high_cag$dosage_CAG)
hussin_low_r2_high_cag$group_dosage_diff<-NA
hussin_low_r2_high_cag[which(hussin_low_r2_high_cag$dosage_diff<0.01),]$group_dosage_diff<-"diff<0.01"
hussin_low_r2_high_cag[which(hussin_low_r2_high_cag$dosage_diff>=0.01 & hussin_low_r2_high_cag$dosage_diff<0.02),]$group_dosage_diff<-"0.1<=diff<0.02"
hussin_low_r2_high_cag[which(hussin_low_r2_high_cag$dosage_diff>0.02),]$group_dosage_diff<-"diff>0.02"



hussin_low_r2_high_cag$group_dosage_diff<-factor(hussin_low_r2_high_cag$group_dosage_diff,levels=c("diff<0.01", "diff>0.02","0.1<=diff<0.02"))


p1<-ggplot(hussin_low_r2_high_cag, aes(x=dosage_CAG, y=dosage_JH, color=group_dosage_diff)) + geom_point(alpha=0.2) +ggtitle(paste0("SNPs with R2<0.05 in JH and R2>=0.3 in CAG ",nrow(hussin_low_r2_high_cag), " SNPs")) + xlab("Average Dosage CAG") + ylab("Average Dosage JH")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/Dosages_low_R2.png", p1)

#Zoom in 
data_dosage_0<-hussin_low_r2_high_cag[which(hussin_low_r2_high_cag$dosage_CAG<0.25 & hussin_low_r2_high_cag$dosage_JH<0.25),]
p1<-ggplot(data_dosage_0, aes(x=dosage_CAG, y=dosage_JH, color=group_dosage_diff)) + geom_point(alpha=0.2) +ggtitle(paste0("SNPs with R2<0.05 in JH and R2>=0.3 in CAG ",nrow(hussin_low_r2_high_cag), " SNPs")) + xlab("Average Dosage CAG") + ylab("Average Dosage JH")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/Dosages_0_low_R2.png", p1)

data_dosage_2<-hussin_low_r2_high_cag[which(hussin_low_r2_high_cag$dosage_CAG>1.75 & hussin_low_r2_high_cag$dosage_JH>1.75),]
p1<-ggplot(data_dosage_2, aes(x=dosage_CAG, y=dosage_JH, color=group_dosage_diff)) + geom_point(alpha=0.2) +ggtitle(paste0("SNPs with R2<0.05 in JH and R2>=0.3 in CAG ",nrow(hussin_low_r2_high_cag), " SNPs")) + xlab("Average Dosage CAG") + ylab("Average Dosage JH")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/Dosages_2_low_R2.png", p1)


#variance instead of mean

p1<-ggplot(hussin_low_r2_high_cag, aes(x=var_dosage_CAG, y=var_dosage_JH)) + geom_point(alpha=0.2) +ggtitle(paste0("SNPs with R2<0.05 in JH and R2>=0.3 in CAG ",nrow(hussin_low_r2_high_cag), " SNPs")) + xlab("Variance CAG") + ylab("Variance JH")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/Variance_dosage_low_R2.png", p1)





p1<-ggplot(hussin_low_r2_high_cag, aes(x=Position, y=dosage_diff)) + geom_point(alpha=0.2)+geom_hline(yintercept=mean(hussin_low_r2_high_cag$dosage_diff), color="red") +ggtitle(paste0("SNPs with R2<0.05 in JH and R2>=0.3 in CAG ",nrow(hussin_low_r2_high_cag), " SNPs")) + xlab("Position") + ylab("Absolute difference in dosages")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/Dosages_sub_low_R2.png", p1)






#FINAL ANALYSES--------------------------

#dist of the SNPs per chromomes
i=1
while(i<=22){ 
	
	print(i)
	data_tmp<-hussin_low_r2_high_cag[which(hussin_low_r2_high_cag$chr==i),]

	
	png(filename = paste0("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/chr",i,"_position_low_R2.png"), width = 2000, height = 500)
	p1<-ggplot(data_tmp, aes(x=as.numeric(pos_on_chr), y=as.numeric(pos_on_chr))) + geom_point(alpha=0.02) +ggtitle(paste0("Chr",i," SNPs with R2<0.05 in JH and R2>=0.3 in CAG ",nrow(data_tmp), " SNPs")) + xlab("Position") + ylab("Position")
	#p1<-ggplot(data_tmp, aes(x=as.numeric(pos_on_chr), y=as.numeric(pos_on_chr))) + geom_bin2d(aes(fill=stat(count))) +scale_fill_continuous(type = "viridis") +ggtitle(paste0("Chr",i," SNPs with R2<0.05 in JH and R2>=0.3 in CAG ",nrow(data_tmp), " SNPs")) + xlab("Position") + ylab("Position")
	#ggsave(paste0("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/chr",i,"_position_low_R2.png"), p1)
	print(p1)
	dev.off()

	i=i+1
}







#LOW JH + HIGH CAG
#dist of the low R2 only in JH in each chromosomesi

#add chr1 last position for each chromosome
tmp_chr1_last<-hussin_low_r2_high_cag[which(hussin_low_r2_high_cag$chromosome==1),]
row_to_add<-tmp_chr1_last[nrow(tmp_chr1_last) ,]  

i=1
while(i<=22){
	
	row_to_add$chromosome<-i
	tmp<-rbind(hussin_low_r2_high_cag, row_to_add)
	hussin_low_r2_high_cag<-tmp
	i=i+1
}


png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/v2Dist_position_in_genome_low_R2.png", width = 2000, height = 2000)
#with the end position of chr1 at the end of each chromosome (loop over)
#png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/Dist_position_in_genome_low_R2.png", width = 2000, height = 2000)
p0<-ggplot(hussin_low_r2_high_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(binwidth=((max(as.numeric(hussin_low_r2_high_cag$pos_on_chr))-min(as.numeric(hussin_low_r2_high_cag$pos_on_chr)))/75),colour="black", fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2>=0.3 in CAG (",nrow(hussin_low_r2_high_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))
#p0<-ggplot(hussin_low_r2_high_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(bins=75,colour="black", fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2>=0.3 in CAG (",nrow(hussin_low_r2_high_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))
#p0<-ggplot(hussin_low_r2_high_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(aes(y=..density..), colour="black", fill="white")+ geom_density(alpha=.2, fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2>=0.3 in CAG (",nrow(hussin_low_r2_high_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))
#print(p0)
print(p0)
dev.off()
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/Dist_position_in_genome_low_R2.png", p0)


png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/Dist_smooth_position_in_genome_low_R2.png", width = 2000, height = 2000)
#p0<-ggplot(hussin_low_r2_high_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(colour="black", fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2>=0.3 in CAG (",nrow(hussin_low_r2_high_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))
p0<-ggplot(hussin_low_r2_high_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(aes(y=..density..), colour="black", fill="white")+ geom_density(alpha=.2, fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2>=0.3 in CAG (",nrow(hussin_low_r2_high_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))
#print(p0)
print(p0)
dev.off()



tmp_chr1_last<-hussin_low_r2_cag[which(hussin_low_r2_cag$chromosome==1),]
row_to_add<-tmp_chr1_last[nrow(tmp_chr1_last) ,]

i=1
while(i<=22){

        row_to_add$chromosome<-i
        tmp<-rbind(hussin_low_r2_cag, row_to_add)
        hussin_low_r2_cag<-tmp
        i=i+1
}


#low CAG+JH
png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/v2Dist_position_in_genome_low_R2_JH_and_CAG.png", width = 2000, height = 2000)
p0<-ggplot(hussin_low_r2_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(binwidth=((max(as.numeric(hussin_low_r2_cag$pos_on_chr))-min(as.numeric(hussin_low_r2_cag$pos_on_chr)))/75),colour="black", fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2<0.3 in CAG (",nrow(hussin_low_r2_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))

#p0<-ggplot(hussin_low_r2_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(bins=75,colour="black", fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2<0.3 in CAG (",nrow(hussin_low_r2_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))
#p0<-ggplot(hussin_low_r2_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(aes(y=..density..), colour="black", fill="white")+ geom_density(alpha=.2, fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2<0.3 in CAG (",nrow(hussin_low_r2_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))
#p0<-ggplot(hussin_low_r2_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(aes(y=..density..), colour="black", fill="white")+ geom_density(alpha=.2, fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2<0.3 in CAG (",nrow(hussin_low_r2_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))
print(p0)
dev.off()







#low CAG+JH
#dist of the low R2 only in JH in each chromosomes
png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/Dist_smooth_position_in_genome_low_R2_JH_and_CAG.png", width = 2000, height = 2000)
#p0<-ggplot(hussin_low_r2_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(colour="black", fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2<0.3 in CAG (",nrow(hussin_low_r2_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))
p0<-ggplot(hussin_low_r2_cag, aes(x=as.numeric(pos_on_chr))) +geom_histogram(aes(y=..density..), colour="black", fill="white")+ geom_density(alpha=.2, fill="#00BFC4") +xlab("Position") + ggtitle(paste0("R2<0.05 in JH and R2<0.3 in CAG (",nrow(hussin_low_r2_cag),"SNPs)")) + facet_wrap(~chromosome,scales = "free_x") + theme( text = element_text(size = 25))
print(p0)
dev.off()

































#where in the genome
p1<-ggplot(hussin_low_r2 ,aes(x=Position, y=R2)) + geom_bin2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis") + xlab("R2 Hussin") 
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2/R2_from_position_in_genome_low_R2.png", p1)




#Accuracy Hussin
occurence_groupPerc<-as.data.frame(table(hussin_low_r2$groupPerc))
names(occurence_groupPerc)<-c("Accuracy_Hussin", "NumberOfSNPs")

png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2_hussin_CAG_vs_WES_reimputed_SNP_accuracy.png", width = 2000, height = 2000)
p2<-ggplot(occurence_groupPerc, aes(x=Accuracy_Hussin, y=NumberOfSNPs, fill=as.factor(Accuracy_Hussin))) +  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) + geom_text(aes(label=NumberOfSNPs), position=position_dodge(width=0.9), size=15)+   theme(legend.position = "none", text = element_text(size = 50))
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_vs_WES_omnireimputed_SNP_accuracy.png",p2)
print(p2)
dev.off()

#accuracy_CAG
occurence_groupPerc<-as.data.frame(table(hussin_low_r2$groupPerc_cag))
names(occurence_groupPerc)<-c("Accuracy_CAG", "NumberOfSNPs")

png(filename = "/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/low_R2_hussin_CAG_vs_WES_reimputed_SNP_accuracy_CAG.png", width = 2000, height = 2000)
p2<-ggplot(occurence_groupPerc, aes(x=Accuracy_CAG, y=NumberOfSNPs, fill=as.factor(Accuracy_CAG))) +  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1)) + geom_text(aes(label=NumberOfSNPs), position=position_dodge(width=0.9), size=15)+   theme(legend.position = "none", text = element_text(size = 50))
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_vs_WES_omnireimputed_SNP_accuracy.png",p2)
print(p2)
dev.off()



#maf wes - maf 192 hussin
p1<-ggplot(hussin_low_r2, aes(x=MAF_WES_minus_MAF_192)) +  geom_histogram() + xlab("MAF_WES - MAF_192_hussin") + ylab("density")# + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/test_low_R2_Hussin_vs_WES_maff_diff.png",p1)

p1<-ggplot(hussin_low_r2, aes(x=MAF_topmed_minus_MAF_192)) +  geom_histogram() + xlab("MAF_WES - MAF_192_hussin") + ylab("density")# + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/test_low_R2_Hussin_vs_topmed_maff_diff.png",p1)
















p1<-ggplot(hussin_low_r2, aes(x=R2_CAG)) +  geom_density() + ylim(c(0,30))+ xlab("R2") + ylab("density") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_vs_WES_dist_R2cag_under0_05.png",p1)





comp_r2<-ggplot(hussin_low_r2, aes(x=R2, y=R2_CAG)) + geom_bin_2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis") + geom_abline(col="red") +  xlab("R2_Hussin") + ylab("R2_CAG") +ggtitle("R2 under 0.05 in Hussin method (55569 SNPs)")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_low_R2.png", comp_r2)



comp_accuracy<-ggplot(hussin_low_r2, aes(x=goodPerc_CAG, y=goodPerc)) + geom_bin_2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis") + geom_abline(col="red") +  xlab("Accuracy CAG") + ylab("Accuracy Hussin") +ggtitle("R2 under 0.05 in Hussin method (55569 SNPs)")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_low_R2_accuracy.png", comp_accuracy)



pMAF_topmed_hussin<-ggplot(hussin_low_r2, aes(x=maf_Hussin, y=maf_CAG)) + geom_point(alpha=0.5)+ geom_abline(col="red") + xlab("MAF Hussin") + ylab("MAF CAG")+ggtitle("R2 under 0.05 in Hussin method (55569 SNPs)")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_low_R2_maf.png",pMAF_topmed_hussin)


pMAF_topmed_hussin<-ggplot(hussin_low_r2, aes(x=maf_192, y=maf_192_cag)) + geom_point(alpha=0.5)+ geom_abline(col="red")+xlim(c(0,0.5))+ylim(c(0,0.5)) + xlab("MAF_192 Hussin") + ylab("MAF_192 CAG")+ggtitle("R2 under 0.05 in Hussin method (55569 SNPs)")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_low_R2_maf_192.png",pMAF_topmed_hussin)





#compare maf vs WES maf
pMAF_WES_hussin<-ggplot(hussin_low_r2, aes(x=maf_WES, y=maf_Hussin)) + geom_point(alpha=0.5)+ geom_abline(col="red") +xlim(c(0,0.5))+ylim(c(0,0.5))+ xlab("MAF WES") + ylab("MAF Hussin")+ggtitle("R2 under 0.05 in Hussin method (55569 SNPs)")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_low_R2_maf_WES.png",pMAF_WES_hussin)


pMAF_WES_CAG<-ggplot(hussin_low_r2, aes(x=maf_WES, y=maf_CAG)) + geom_point(alpha=0.5)+ geom_abline(col="red") +xlim(c(0,0.5))+ylim(c(0,0.5))+ xlab("MAF WES") + ylab("MAF CAG")+ggtitle("R2 under 0.05 in Hussin method (55569 SNPs)")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/CAG_low_R2_maf_WES.png",pMAF_WES_CAG)



#compare maf_192 vs WES maf
pMAF_WES_hussin<-ggplot(hussin_low_r2, aes(x=maf_WES, y=maf_192, color=groupR2)) + geom_point(alpha=0.5)+ geom_abline(col="red")+xlim(c(0,0.5))+ylim(c(0,0.5)) + xlab("MAF WES") + ylab("MAF_192 Hussin")+ggtitle("R2 under 0.05 in Hussin method (55569 SNPs)")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_low_R2_maf_192_WES_color_R2.png",pMAF_WES_hussin)
#pMAF_WES_hussin<-ggplot(hussin_low_r2, aes(x=maf_WES, y=maf_192, color=groupPerc)) + geom_point(alpha=0.5)+ geom_abline(col="red")+xlim(c(0,0.5))+ylim(c(0,0.5)) + xlab("MAF WES") + ylab("MAF_192 Hussin")+ggtitle("R2 under 0.05 in Hussin method (55569 SNPs)")
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_low_R2_maf_192_WES.png",pMAF_WES_hussin)


pMAF_WES_CAG<-ggplot(hussin_low_r2, aes(x=maf_WES, y=maf_192_cag, color=groupR2_cag)) + geom_point(alpha=0.5)+ geom_abline(col="red") +xlim(c(0,0.5))+ylim(c(0,0.5))+ xlab("MAF WES") + ylab("MAF_192 CAG")+ggtitle("R2 under 0.05 in Hussin method (55569 SNPs)")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/CAG_low_R2_maf_192_WES_color_R2.png",pMAF_WES_CAG)
#pMAF_WES_CAG<-ggplot(hussin_low_r2, aes(x=maf_WES, y=maf_192_cag, color=groupPerc_cag)) + geom_point(alpha=0.5)+ geom_abline(col="red") +xlim(c(0,0.5))+ylim(c(0,0.5))+ xlab("MAF WES") + ylab("MAF_192 CAG")+ggtitle("R2 under 0.05 in Hussin method (55569 SNPs)")
#ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/CAG_low_R2_maf_192_WES.png",pMAF_WES_CAG)
























#----------------------------------------STRATIFICATION MAF-----------------------------------------------

#add info in a single file
modified_data_final_hussin$groupPerc_cag<-modified_data_final_CAG$groupPerc[match(modified_data_final_CAG$SNP, modified_data_final_hussin$SNP)]
modified_data_final_hussin$goodPerc_cag<-modified_data_final_CAG$goodPerc[match(modified_data_final_CAG$SNP, modified_data_final_hussin$SNP)]
modified_data_final_hussin$maf_192_cag<-modified_data_final_CAG$maf_192[match(modified_data_final_CAG$SNP, modified_data_final_hussin$SNP)]
modified_data_final_hussin$R2_cag<-modified_data_final_CAG$R2[match(modified_data_final_CAG$SNP, modified_data_final_hussin$SNP)]
modified_data_final_hussin$groupR2_cag<-modified_data_final_CAG$groupR2[match(modified_data_final_CAG$SNP, modified_data_final_hussin$SNP)]
modified_data_final_hussin$diff_accuracy<-(modified_data_final_hussin$goodPerc-modified_data_final_hussin$goodPerc_cag)



#-------------------Accuracy JH-CAG

p1<-ggplot(modified_data_final_hussin, aes(x=Position, y=diff_accuracy)) + geom_bin2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis") + xlab("Chromosome and Position") +ylab("JH accuracy - CAG accuracy")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/Hussin_vs_CAG_accuracy_substraction.png",p1)






#--------------------------density R2 vs Accuracy-----------------------------

tmp_R2<-data.frame(modified_data_final_hussin$SNP, modified_data_final_hussin$R2)
names(tmp_R2)<-c("SNP", "Exactitude")
tmp_R2$Exactitude<-tmp_R2$Exactitude*100
tmp_R2$Metric="R2 * 100"

tmp_accuracy<-data.frame(modified_data_final_hussin$SNP, modified_data_final_hussin$goodPerc)
names(tmp_accuracy)<-c("SNP", "Exactitude")
tmp_accuracy$Metric="Accuracy %"

data_exactitude<-rbind(tmp_R2, tmp_accuracy)

p1<-ggplot(data_exactitude, aes(x=Exactitude, color=Metric)) +  geom_density() + xlab("R2") + ylab("density") + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))+   theme(legend.position = "none")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_cag/Hussin_vs_WES_dist_R2_and_accuracy.png",p1)







#----------------group MAF< 0.01 -------------------




WESmaf_001_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_WES<0.01),]
WESmaf_001_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_WES<0.01),]

#----Chi square---


WESmaf_001_hussin$method="Impute-Merge"
WESmaf_001_cag$method="Merge-Impute"

WESmaf_001_hussin$perfect_accuracy=NA
WESmaf_001_cag$perfect_accuracy=NA


WESmaf_001_hussin[which(WESmaf_001_hussin$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="100%"
WESmaf_001_hussin[which(WESmaf_001_hussin$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="under 100%"

WESmaf_001_cag[which(WESmaf_001_cag$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="100%"
WESmaf_001_cag[which(WESmaf_001_cag$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="under 100%"

tmp_WESmaf_001_hussin<-WESmaf_001_hussin[,1:17]
tmp_WESmaf_001_cag<-WESmaf_001_cag[,1:17]

merged_data<-rbind(tmp_WESmaf_001_hussin, tmp_WESmaf_001_cag)

table(merged_data$method , merged_data$perfect_accuracy)

test<-chisq.test(merged_data$method , merged_data$perfect_accuracy, correct=FALSE)




nb_better_imputed<-nrow(WESmaf_001_hussin[which(WESmaf_001_hussin$diff_accuracy>0),])
nb_worse_imputed<-nrow(WESmaf_001_hussin[which(WESmaf_001_hussin$diff_accuracy<0),])
nb_equal_imputed<-nrow(WESmaf_001_hussin[which(WESmaf_001_hussin$diff_accuracy==0),])
total<-nrow(WESmaf_001_hussin)
total<-nrow(WESmaf_001_hussin)
perc_better_imputed_JH<-(nb_better_imputed/total)*100
perc_worse_imputed_JH<-(nb_worse_imputed/total)*100
perc_equal_imputed_JH<-(nb_equal_imputed/total)*100
summary<-data.frame(c("accuracy_higher_in_JH", "accuracy_lower_in_JH", "accuracy_equal_in_JH"), c(perc_better_imputed_JH, perc_worse_imputed_JH, perc_equal_imputed_JH))






group_1_accuracy<-ggplot(WESmaf_001_hussin, aes(x=goodPerc, y=goodPerc_cag)) + geom_point(alpha=0.5)+ geom_abline(col="red") +xlim(c(0,100)) +ylim(c(0,100))+ xlab("Accuracy JH") + ylab("Accuracy CAG") + ggtitle(" 0.01 <= MAF in WES < 0.05")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/group_1_accuracy.png",group_1_accuracy)


group_1_accuracy<-ggplot(WESmaf_001_hussin, aes(x=goodPerc, y=goodPerc_cag)) + geom_bin_2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis")+ geom_abline(col="red") +xlim(c(0,100)) +ylim(c(0,100))+ xlab("Accuracy JH") + ylab("Accuracy CAG") + ggtitle(" 0.01 <= MAF in WES < 0.05")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/group_1_accuracy_heatmap.png",group_1_accuracy)






#----------------group 0.01 <= MAF <0.05--------------------



WESmaf_001_005_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_WES>=0.01 & modified_data_final_hussin$maf_WES<0.05),]
WESmaf_001_005_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_WES>=0.01 & modified_data_final_CAG$maf_WES<0.05),]


#----Chi square---


WESmaf_001_005_hussin$method="Impute-Merge"
WESmaf_001_005_cag$method="Merge-Impute"

WESmaf_001_005_hussin$perfect_accuracy=NA
WESmaf_001_005_cag$perfect_accuracy=NA


WESmaf_001_005_hussin[which(WESmaf_001_005_hussin$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="100%"
WESmaf_001_005_hussin[which(WESmaf_001_005_hussin$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="under 100%"

WESmaf_001_005_cag[which(WESmaf_001_005_cag$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="100%"
WESmaf_001_005_cag[which(WESmaf_001_005_cag$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="under 100%"

tmp_WESmaf_001_005_hussin<-WESmaf_001_005_hussin[,1:17]
tmp_WESmaf_001_005_cag<-WESmaf_001_005_cag[,1:17]

merged_data<-rbind(tmp_WESmaf_001_005_hussin, tmp_WESmaf_001_005_cag)

table(merged_data$method , merged_data$perfect_accuracy)

test<-chisq.test(merged_data$method , merged_data$perfect_accuracy, correct=FALSE)







nb_better_imputed<-nrow(WESmaf_001_005_hussin[which(WESmaf_001_005_hussin$diff_accuracy>0),])
nb_worse_imputed<-nrow(WESmaf_001_005_hussin[which(WESmaf_001_005_hussin$diff_accuracy<0),])
nb_equal_imputed<-nrow(WESmaf_001_005_hussin[which(WESmaf_001_005_hussin$diff_accuracy==0),])
total<-nrow(WESmaf_001_005_hussin)
total<-nrow(WESmaf_001_005_hussin)
perc_better_imputed_JH<-(nb_better_imputed/total)*100
perc_worse_imputed_JH<-(nb_worse_imputed/total)*100
perc_equal_imputed_JH<-(nb_equal_imputed/total)*100
summary<-data.frame(c("accuracy_higher_in_JH", "accuracy_lower_in_JH", "accuracy_equal_in_JH"), c(perc_better_imputed_JH, perc_worse_imputed_JH, perc_equal_imputed_JH)) 






group_1_accuracy<-ggplot(WESmaf_001_005_hussin, aes(x=goodPerc, y=goodPerc_cag)) + geom_point(alpha=0.5)+ geom_abline(col="red") +xlim(c(0,100)) +ylim(c(0,100))+ xlab("Accuracy JH") + ylab("Accuracy CAG") + ggtitle(" 0.01 <= MAF in WES < 0.05")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/group_1_accuracy.png",group_1_accuracy)


group_1_accuracy<-ggplot(WESmaf_001_005_hussin, aes(x=goodPerc, y=goodPerc_cag)) + geom_bin_2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis")+ geom_abline(col="red") +xlim(c(0,100)) +ylim(c(0,100))+ xlab("Accuracy JH") + ylab("Accuracy CAG") + ggtitle(" 0.01 <= MAF in WES < 0.05")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/group_1_accuracy_heatmap.png",group_1_accuracy)






#----------------group 0.05 <= MAF <0.1--------------------



WESmaf_005_01_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_WES>=0.05 & modified_data_final_hussin$maf_WES<0.1),]
WESmaf_005_01_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_WES>=0.05 & modified_data_final_CAG$maf_WES<0.1),]


#----Chi square---


WESmaf_005_01_hussin$method="Impute-Merge"
WESmaf_005_01_cag$method="Merge-Impute"

WESmaf_005_01_hussin$perfect_accuracy=NA
WESmaf_005_01_cag$perfect_accuracy=NA


WESmaf_005_01_hussin[which(WESmaf_005_01_hussin$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="100%"
WESmaf_005_01_hussin[which(WESmaf_005_01_hussin$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="under 100%"

WESmaf_005_01_cag[which(WESmaf_005_01_cag$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="100%"
WESmaf_005_01_cag[which(WESmaf_005_01_cag$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="under 100%"

tmp_WESmaf_005_01_hussin<-WESmaf_005_01_hussin[,1:17]
tmp_WESmaf_005_01_cag<-WESmaf_005_01_cag[,1:17]

merged_data<-rbind(tmp_WESmaf_005_01_hussin, tmp_WESmaf_005_01_cag)

table(merged_data$method , merged_data$perfect_accuracy)

test<-chisq.test(merged_data$method , merged_data$perfect_accuracy, correct=FALSE)




nb_better_imputed<-nrow(WESmaf_005_01_hussin[which(WESmaf_005_01_hussin$diff_accuracy>0),])
nb_worse_imputed<-nrow(WESmaf_005_01_hussin[which(WESmaf_005_01_hussin$diff_accuracy<0),])
nb_equal_imputed<-nrow(WESmaf_005_01_hussin[which(WESmaf_005_01_hussin$diff_accuracy==0),])
total<-nrow(WESmaf_005_01_hussin)
total<-nrow(WESmaf_005_01_hussin)
perc_better_imputed_JH<-(nb_better_imputed/total)*100
perc_worse_imputed_JH<-(nb_worse_imputed/total)*100
perc_equal_imputed_JH<-(nb_equal_imputed/total)*100


summary<-data.frame(c("accuracy_higher_in_JH", "accuracy_lower_in_JH", "accuracy_equal_in_JH"), c(perc_better_imputed_JH, perc_worse_imputed_JH, perc_equal_imputed_JH))



group_1_accuracy<-ggplot(WESmaf_005_01_hussin, aes(x=goodPerc, y=goodPerc_cag)) + geom_point(alpha=0.5)+ geom_abline(col="red") +xlim(c(0,100)) +ylim(c(0,100))+ xlab("Accuracy JH") + ylab("Accuracy CAG") + ggtitle(" 0.05 <= MAF in WES < 0.1")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/group_2_accuracy.png",group_1_accuracy)


group_1_accuracy<-ggplot(WESmaf_005_01_hussin, aes(x=goodPerc, y=goodPerc_cag)) + geom_bin_2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis")+ geom_abline(col="red") +xlim(c(0,100)) +ylim(c(0,100))+ xlab("Accuracy JH") + ylab("Accuracy CAG") + ggtitle(" 0.05 <= MAF in WES < 0.1")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/group_2_accuracy_heatmap.png",group_1_accuracy)






#----------------group 0.1 <= MAF <0.25--------------------



WESmaf_01_025_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_WES>=0.1 & modified_data_final_hussin$maf_WES<0.25),]
WESmaf_01_025_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_WES>=0.1 & modified_data_final_CAG$maf_WES<0.25),]


#----Chi square---


WESmaf_01_025_hussin$method="Impute-Merge"
WESmaf_01_025_cag$method="Merge-Impute"

WESmaf_01_025_hussin$perfect_accuracy=NA
WESmaf_01_025_cag$perfect_accuracy=NA


WESmaf_01_025_hussin[which(WESmaf_01_025_hussin$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="100%"
WESmaf_01_025_hussin[which(WESmaf_01_025_hussin$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="under 100%"

WESmaf_01_025_cag[which(WESmaf_01_025_cag$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="100%"
WESmaf_01_025_cag[which(WESmaf_01_025_cag$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="under 100%"

tmp_WESmaf_01_025_hussin<-WESmaf_01_025_hussin[,1:17]
tmp_WESmaf_01_025_cag<-WESmaf_01_025_cag[,1:17]

merged_data<-rbind(tmp_WESmaf_01_025_hussin, tmp_WESmaf_01_025_cag)

table(merged_data$method , merged_data$perfect_accuracy)

test<-chisq.test(merged_data$method , merged_data$perfect_accuracy, correct=FALSE)



nb_better_imputed<-nrow(WESmaf_01_025_hussin[which(WESmaf_01_025_hussin$diff_accuracy>0),])
nb_worse_imputed<-nrow(WESmaf_01_025_hussin[which(WESmaf_01_025_hussin$diff_accuracy<0),])
nb_equal_imputed<-nrow(WESmaf_01_025_hussin[which(WESmaf_01_025_hussin$diff_accuracy==0),])
total<-nrow(WESmaf_01_025_hussin)
total<-nrow(WESmaf_01_025_hussin)
perc_better_imputed_JH<-(nb_better_imputed/total)*100
perc_worse_imputed_JH<-(nb_worse_imputed/total)*100
perc_equal_imputed_JH<-(nb_equal_imputed/total)*100


summary<-data.frame(c("accuracy_higher_in_JH", "accuracy_lower_in_JH", "accuracy_equal_in_JH"), c(perc_better_imputed_JH, perc_worse_imputed_JH, perc_equal_imputed_JH))



group_1_accuracy<-ggplot(WESmaf_01_025_hussin, aes(x=goodPerc, y=goodPerc_cag)) + geom_point(alpha=0.5)+ geom_abline(col="red") +xlim(c(0,100)) +ylim(c(0,100))+ xlab("Accuracy JH") + ylab("Accuracy CAG") + ggtitle(" 0.1 <= MAF in WES < 0.25")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/group_3_accuracy.png",group_1_accuracy)


group_1_accuracy<-ggplot(WESmaf_01_025_hussin, aes(x=goodPerc, y=goodPerc_cag)) + geom_bin_2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis")+ geom_abline(col="red") +xlim(c(0,100)) +ylim(c(0,100))+ xlab("Accuracy JH") + ylab("Accuracy CAG") + ggtitle(" 0.1 <= MAF in WES < 0.25")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/group_3_accuracy_heatmap.png",group_1_accuracy)





#----------------group 0.25 <= MAF <0.5--------------------



WESmaf_025_05_hussin<-modified_data_final_hussin[which(modified_data_final_hussin$maf_WES>=0.25 & modified_data_final_hussin$maf_WES<=0.5),]
WESmaf_025_05_cag<-modified_data_final_CAG[which(modified_data_final_CAG$maf_WES>=0.25 & modified_data_final_CAG$maf_WES<=0.5),]


#----Chi square---


WESmaf_025_05_hussin$method="Impute-Merge"
WESmaf_025_05_cag$method="Merge-Impute"

WESmaf_025_05_hussin$perfect_accuracy=NA
WESmaf_025_05_cag$perfect_accuracy=NA


WESmaf_025_05_hussin[which(WESmaf_025_05_hussin$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="100%"
WESmaf_025_05_hussin[which(WESmaf_025_05_hussin$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="under 100%"

WESmaf_025_05_cag[which(WESmaf_025_05_cag$groupPerc=="98%<Accuracy<=100%"),]$perfect_accuracy="100%"
WESmaf_025_05_cag[which(WESmaf_025_05_cag$groupPerc!="98%<Accuracy<=100%"),]$perfect_accuracy="under 100%"

tmp_WESmaf_025_05_hussin<-WESmaf_025_05_hussin[,1:17]
tmp_WESmaf_025_05_cag<-WESmaf_025_05_cag[,1:17]

merged_data<-rbind(tmp_WESmaf_025_05_hussin, tmp_WESmaf_025_05_cag)

table(merged_data$method , merged_data$perfect_accuracy)

test<-chisq.test(merged_data$method , merged_data$perfect_accuracy, correct=FALSE)




nb_better_imputed<-nrow(WESmaf_025_05_hussin[which(WESmaf_025_05_hussin$diff_accuracy>0),])
nb_worse_imputed<-nrow(WESmaf_025_05_hussin[which(WESmaf_025_05_hussin$diff_accuracy<0),])
nb_equal_imputed<-nrow(WESmaf_025_05_hussin[which(WESmaf_025_05_hussin$diff_accuracy==0),])
total<-nrow(WESmaf_025_05_hussin)
total<-nrow(WESmaf_025_05_hussin)
perc_better_imputed_JH<-(nb_better_imputed/total)*100
perc_worse_imputed_JH<-(nb_worse_imputed/total)*100
perc_equal_imputed_JH<-(nb_equal_imputed/total)*100


summary<-data.frame(c("accuracy_higher_in_JH", "accuracy_lower_in_JH", "accuracy_equal_in_JH"), c(perc_better_imputed_JH, perc_worse_imputed_JH, perc_equal_imputed_JH))



group_1_accuracy<-ggplot(WESmaf_025_05_hussin, aes(x=goodPerc, y=goodPerc_cag)) + geom_point(alpha=0.5)+ geom_abline(col="red") +xlim(c(0,100)) +ylim(c(0,100))+ xlab("Accuracy JH") + ylab("Accuracy CAG") + ggtitle(" 0.25 <= MAF in WES < 0.5")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/group_4_accuracy.png",group_1_accuracy)


group_1_accuracy<-ggplot(WESmaf_005_01_hussin, aes(x=goodPerc, y=goodPerc_cag)) + geom_bin_2d(aes(fill=stat(log(count)))) +scale_fill_continuous(type = "viridis")+ geom_abline(col="red") +xlim(c(0,100)) +ylim(c(0,100))+ xlab("Accuracy JH") + ylab("Accuracy CAG") + ggtitle(" 0.25 <= MAF in WES < 0.5")
ggsave("/lustre07/scratch/justinp/topmed_new/CAG_new/WES/CAG_MERGE/MERGE_hussin/group_4_accuracy_heatmap.png",group_1_accuracy)


