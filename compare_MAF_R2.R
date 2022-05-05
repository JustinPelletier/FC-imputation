library(ggplot2)




topmed<-read.table("/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_all_arrays_09122021/maf0.01.hwe0.000001.miss5perc/Imputed/merging/R2_computation/noArchi/probsMetrics_files/chr1_filtered_all_chr_R2_0.0_Hussin.vcf.maf_r2", header=F)
names(topmed)<-c("SNP", "MAF", "R2")


hrc<-read.table("/lustre06/project/6065672/shared/Cartagene/Genotypes_and_Phenotypes/Imputation_all_arrays_09122021/maf0.01.hwe0.000001.miss5perc/Imputed/HRC/merge/probsMetrics_files/chr1_filtered_all_chr_R2_0.0_Hussin.vcf.maf_r2", header=F)
names(hrc)<-c("SNP", "MAF", "R2")

hrc$panel<-"HRC"
topmed$panel<-"TOPMed"


data<-rbind(hrc, topmed)




#dist in rare variants of R2
rare_hrc<-hrc[which(hrc$MAF<0.01),]
rare_topmed<-topmed[which(topmed$MAF<0.01),]
rare_data<-rbind(rare_hrc, rare_topmed)
nrow(rare_topmed[which(rare_topmed$R2<0.3),])/nrow(rare_topmed)



png("/lustre07/scratch/justinp/topmed_new/imputation_pipeline/HRC_vs_TOPMed/distribution_HRC_TOPMED_low_MAF.png", width=2300, height=1200)
p1<-ggplot(rare_data, aes(x=R2, fill=panel)) + geom_histogram(color="black")+theme_classic() +ylab("Number of imputed variants")+facet_wrap(~panel, scales = "free_y") + theme(text = element_text(size = 50)) + theme(legend.title=element_blank())
print(p1)
#ggsave("/lustre07/scratch/justinp/topmed_new/imputation_pipeline/HRC_vs_TOPMed/comparison_HRC_TOPMED_MAF_and_R2.png", p1)
dev.off()







#by maf group

mean_maf<-c()
mean_r2<-c()
data_set<-c()

mean_maf<-c(mean_maf, mean(hrc[which(hrc$MAF<0.01 & hrc$MAF>0),]$MAF))
mean_r2<-c(mean_r2, mean(hrc[which(hrc$MAF<0.01 & hrc$MAF>0),]$R2))
data_set<-c(data_set, "HRC")
mean_maf<-c(mean_maf, mean(hrc[which(hrc$MAF>=0.01&hrc$MAF<0.05),]$MAF))
mean_r2<-c(mean_r2, mean(hrc[which(hrc$MAF>=0.01&hrc$MAF<0.05),]$R2))
data_set<-c(data_set, "HRC")
mean_maf<-c(mean_maf, mean(hrc[which(hrc$MAF>=0.05&hrc$MAF<0.1),]$MAF))
mean_r2<-c(mean_r2, mean(hrc[which(hrc$MAF>=0.05&hrc$MAF<0.1),]$R2))
data_set<-c(data_set, "HRC")
mean_maf<-c(mean_maf, mean(hrc[which(hrc$MAF>=0.1&hrc$MAF<0.25),]$MAF))
mean_r2<-c(mean_r2, mean(hrc[which(hrc$MAF>=0.1&hrc$MAF<0.25),]$R2))
data_set<-c(data_set, "HRC")
mean_maf<-c(mean_maf, mean(hrc[which(hrc$MAF>=0.25&hrc$MAF<=0.5),]$MAF))
mean_r2<-c(mean_r2, mean(hrc[which(hrc$MAF>=0.25&hrc$MAF<=0.5),]$R2))
data_set<-c(data_set, "HRC")


mean_maf<-c(mean_maf, mean(topmed[which(topmed$MAF<0.01 & topmed$MAF>0),]$MAF))
mean_r2<-c(mean_r2, mean(topmed[which(topmed$MAF<0.01 & topmed$MAF>0),]$R2))
data_set<-c(data_set, "TOPMed")
mean_maf<-c(mean_maf, mean(topmed[which(topmed$MAF>=0.01&topmed$MAF<0.05),]$MAF))
mean_r2<-c(mean_r2, mean(topmed[which(topmed$MAF>=0.01&topmed$MAF<0.05),]$R2))
data_set<-c(data_set, "TOPMed")
mean_maf<-c(mean_maf, mean(topmed[which(topmed$MAF>=0.05&topmed$MAF<0.1),]$MAF))
mean_r2<-c(mean_r2, mean(topmed[which(topmed$MAF>=0.05&topmed$MAF<0.1),]$R2))
data_set<-c(data_set, "TOPMed")
mean_maf<-c(mean_maf, mean(topmed[which(topmed$MAF>=0.1&topmed$MAF<0.25),]$MAF))
mean_r2<-c(mean_r2, mean(topmed[which(topmed$MAF>=0.1&topmed$MAF<0.25),]$R2))
data_set<-c(data_set, "TOPMed")
mean_maf<-c(mean_maf, mean(topmed[which(topmed$MAF>=0.25&topmed$MAF<=0.5),]$MAF))
mean_r2<-c(mean_r2, mean(topmed[which(topmed$MAF>=0.25&topmed$MAF<=0.5),]$R2))
data_set<-c(data_set, "TOPMed")




data_group<-data.frame(mean_maf, mean_r2, data_set)
data_group$colour<-0
data_group[which(data_group$data_set=="TOPMed"),]$colour<-"#619CFF"
data_group[which(data_group$data_set=="HRC"),]$colour<-"#00BA38"


png("/lustre07/scratch/justinp/topmed_new/imputation_pipeline/HRC_vs_TOPMed/comparison_HRC_TOPMED_MAF_and_R2.png", width=2300, height=1200)
p1<-ggplot(data_group, aes(x=mean_maf, y=mean_r2, color=data_set)) + geom_point(size=12) + geom_line(aes(color=data_set),size=3) +theme_classic() +xlab("MAF") + ylab("R2")+ labs(fill = "")+ theme(text = element_text(size = 75)) + theme(legend.title=element_blank())+scale_color_manual(values=c("#00BA38", "#619CFF"))
print(p1)
#ggsave("/lustre07/scratch/justinp/topmed_new/imputation_pipeline/HRC_vs_TOPMed/comparison_HRC_TOPMED_MAF_and_R2.png", p1)
dev.off()


#% of good SNP

mean_maf<-c()
mean_r2<-c()
data_set<-c()

mean_maf<-c(mean_maf, mean(hrc[which(hrc$MAF<0.01),]$MAF))
mean_r2<-c(mean_r2, mean(hrc[which(hrc$MAF<0.01),]$R2))
data_set<-c(data_set, "HRC")
mean_maf<-c(mean_maf, mean(hrc[which(hrc$MAF>=0.01&hrc$MAF<0.05),]$MAF))
mean_r2<-c(mean_r2, mean(hrc[which(hrc$MAF>=0.01&hrc$MAF<0.05),]$R2))
data_set<-c(data_set, "HRC")
mean_maf<-c(mean_maf, mean(hrc[which(hrc$MAF>=0.05&hrc$MAF<0.1),]$MAF))
mean_r2<-c(mean_r2, mean(hrc[which(hrc$MAF>=0.05&hrc$MAF<0.1),]$R2))
data_set<-c(data_set, "HRC")
mean_maf<-c(mean_maf, mean(hrc[which(hrc$MAF>=0.1&hrc$MAF<0.25),]$MAF))
mean_r2<-c(mean_r2, mean(hrc[which(hrc$MAF>=0.1&hrc$MAF<0.25),]$R2))
data_set<-c(data_set, "HRC")
mean_maf<-c(mean_maf, mean(hrc[which(hrc$MAF>=0.25&hrc$MAF<=0.5),]$MAF))
mean_r2<-c(mean_r2, mean(hrc[which(hrc$MAF>=0.25&hrc$MAF<=0.5),]$R2))
data_set<-c(data_set, "HRC")


mean_maf<-c(mean_maf, mean(topmed[which(topmed$MAF<0.01),]$MAF))
mean_r2<-c(mean_r2, mean(topmed[which(topmed$MAF<0.01),]$R2))
data_set<-c(data_set, "TOPMed")
mean_maf<-c(mean_maf, mean(topmed[which(topmed$MAF>=0.01&topmed$MAF<0.05),]$MAF))
mean_r2<-c(mean_r2, mean(topmed[which(topmed$MAF>=0.01&topmed$MAF<0.05),]$R2))
data_set<-c(data_set, "TOPMed")
mean_maf<-c(mean_maf, mean(topmed[which(topmed$MAF>=0.05&topmed$MAF<0.1),]$MAF))
mean_r2<-c(mean_r2, mean(topmed[which(topmed$MAF>=0.05&topmed$MAF<0.1),]$R2))
data_set<-c(data_set, "TOPMed")
mean_maf<-c(mean_maf, mean(topmed[which(topmed$MAF>=0.1&topmed$MAF<0.25),]$MAF))
mean_r2<-c(mean_r2, mean(topmed[which(topmed$MAF>=0.1&topmed$MAF<0.25),]$R2))
data_set<-c(data_set, "TOPMed")
mean_maf<-c(mean_maf, mean(topmed[which(topmed$MAF>=0.25&topmed$MAF<=0.5),]$MAF))
mean_r2<-c(mean_r2, mean(topmed[which(topmed$MAF>=0.25&topmed$MAF<=0.5),]$R2))
data_set<-c(data_set, "TOPMed")






#% that passes the R3>0.3 per maf class
(nrow((hrc[which(hrc$MAF<0.01 & hrc$R2>=0.3),]))/nrow(hrc[which(hrc$MAF<0.01),]))*100
(nrow((topmed[which(topmed$MAF<0.01 & topmed$R2>=0.3),]))/nrow(topmed[which(topmed$MAF<0.01),]))*100
(nrow((hrc[which(hrc$MAF>=0.01 & hrc$R2>=0.3),]))/nrow(hrc[which(hrc$MAF>=0.01),]))*100
(nrow((topmed[which(topmed$MAF>=0.01 & topmed$R2>=0.3),]))/nrow(topmed[which(topmed$MAF>=0.01),]))*100



















nrow(topmed[which(topmed$MAF<0.01),])
nrow(hrc[which(hrc$MAF<0.01),])









#p<-ggplot(data, aes(x=MAF, y=R2, group=panel)) + geom_line(aes(color=panel)) + scale_fill_discrete(name = "Imputation")
#ggsave("/lustre07/scratch/justinp/topmed_new/imputation_pipeline/HRC_vs_TOPMed/comparison_HRC_TOPMED_MAF_and_R2.png", p)











