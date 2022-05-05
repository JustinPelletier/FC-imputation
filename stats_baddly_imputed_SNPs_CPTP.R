library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

#CHR     pos     ID             REF     alt      ?       filter    info     R2     mafCAG   mafATL    mafATP    mafBCGP    mafOHS


CAG<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/badly_imputed_in_only_CAG.vcf.tmp", header=F)
names(CAG)<-c("chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "R2", "maf_ATL", "maf_ATP", "maf_BCGP", "maf_OHS", "maf_CAG")
#names(CAG)<-c("chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "maf", "ATL", "ATP", "BCGP", "OHS")
chr_new<-str_replace_all(CAG$chr, "chr", "")
chr_pos<-paste(chr_new, CAG$pos, sep = ".", collapse = NULL)
CAG$chr_pos<-chr_pos

OHS<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/badly_imputed_in_only_OHS.vcf.tmp", header=F)
names(OHS)<-c("chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "R2", "maf_CAG","maf_ATL", "maf_ATP", "maf_BCGP", "maf_OHS")
#names(CAG)<-c("chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "maf","CAG", "ATL", "ATP", "BCGP")
chr_new<-str_replace_all(OHS$chr, "chr", "")
chr_pos<-paste(chr_new, OHS$pos, sep = ".", collapse = NULL)
OHS$chr_pos<-chr_pos


BCGP<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/badly_imputed_in_only_BCGP.vcf.tmp", header=F)
names(BCGP)<-c("chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "R2", "maf_CAG","maf_ATL", "maf_ATP","maf_OHS", "maf_BCGP")
#names(CAG)<-c("chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "maf", "CAG", "ATL", "ATP", "OHS")
chr_new<-str_replace_all(BCGP$chr, "chr", "")
chr_pos<-paste(chr_new, BCGP$pos, sep = ".", collapse = NULL)
BCGP$chr_pos<-chr_pos


ATL<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/badly_imputed_in_only_ATL.vcf.tmp", header=F)
names(ATL)<-c("chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "R2", "maf_CAG", "maf_ATP", "maf_BCGP", "maf_OHS", "maf_ATL")
#names(CAG)<-c("chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "maf", "CAG", "ATP", "BCGP", "OHS")
chr_new<-str_replace_all(ATL$chr, "chr", "")
chr_pos<-paste(chr_new, ATL$pos, sep = ".", collapse = NULL)
ATL$chr_pos<-chr_pos


ATP<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/badly_imputed_in_only_ATP.vcf.tmp", header=F)
names(ATP)<-c("chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "R2", "maf_CAG", "maf_ATL", "maf_BCGP", "maf_OHS", "maf_ATP")
#names(CAG)<-c("chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "maf", "CAG", "ATL", "BCGP", "OHS")
chr_new<-str_replace_all(ATP$chr, "chr", "")
chr_pos<-paste(chr_new, ATP$pos, sep = ".", collapse = NULL)
ATP$chr_pos<-chr_pos






#------------------------------------------------Distribution of R2 in baddly imputed SNPs---------------------------------------------------------------------#



#plot_CAG<-ggplot(CAG, aes(x=R2)) +  geom_bar(aes(y = (..count..)/sum(..count..))) + xlab("R2") + ylab("Percentage") #+ ylim(c(0,200000))  #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
plot_CAG<-ggplot(CAG, aes(x=R2)) + geom_histogram(color="black", fill="#00BCF4") + xlab("R2") + ylab("Count") + ylim(c(0,160000)) +ggtitle("CAG")  #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/CAG_dist_R2_badly_imputed.png", plot_CAG)



plot_OHS<-ggplot(OHS, aes(x=R2)) + geom_histogram(color="black", fill="orchid4") +xlab("R2") + ylab("Count") + ylim(c(0,160000))  +ggtitle("OHS") #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/OHS_dist_R2_badly_imputed.png", plot_OHS)



plot_BCGP<-ggplot(BCGP, aes(x=R2)) + geom_histogram(color="black", fill="seagreen4") + xlab("R2") + ylab("Count") + ylim(c(0,160000)) +ggtitle("BCGP")  #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/BCGP_dist_R2_badly_imputed.png", plot_BCGP)



plot_ATL<-ggplot(ATL, aes(x=R2)) + geom_histogram(color="black", fill="firebrick2") + xlab("R2") + ylab("Count") + ylim(c(0,160000))  +ggtitle("ATL") #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/ATL_dist_R2_badly_imputed.png", plot_ATL)



plot_ATP<-ggplot(ATP, aes(x=R2)) + geom_histogram(color="black", fill="orange1")   + xlab("R2") + ylab("Count") + ylim(c(0,160000)) +ggtitle("ATP") #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/ATP_dist_R2_badly_imputed.png", plot_ATP)





#----------------------------------------------MAF plots---------------------------------------------------------------


CAG <- as_tibble(CAG)
colnames(CAG)
col_order <- c("chr_pos","chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "R2", "maf_ATL", "maf_ATP", "maf_BCGP", "maf_OHS", "maf_CAG")
CAG2 <- data.frame(CAG[, col_order])

#need to duplicate the columns into row and add dataset for all the MAF columns to put un a geom_bar
data_CAG<-NA
for (i in c(11,12,13,14,15)){
	
	
	print(i)
	data_tmp<-CAG2[,1:10]
	vector_tmp<-CAG2[,i]
	data_CAG_tmp<-cbind(data_tmp, vector_tmp)
	data_CAG<-rbind(data_CAG, data_CAG_tmp)

}
#add dataset column to the dataframe
data_set<-c(rep("ATL", nrow(CAG2)),rep("ATP", nrow(CAG2)),rep("BCGP", nrow(CAG2)),rep("OHS", nrow(CAG2)),rep("CAG", nrow(CAG2)))
data_CAG2 = data_CAG[-1,]
data_CAG2$dataset<-data_set



OHS <- as_tibble(OHS)
colnames(OHS)
col_order <- c("chr_pos","chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "R2", "maf_ATL", "maf_ATP", "maf_BCGP", "maf_OHS", "maf_CAG")
OHS2 <- data.frame(OHS[, col_order])

#need to duplicate the columns into row and add dataset for all the MAF columns to put un a geom_bar
data_OHS<-NA
for (i in c(11,12,13,14,15)){


        print(i)
        data_tmp<-OHS2[,1:10]
        vector_tmp<-OHS2[,i]
        data_OHS_tmp<-cbind(data_tmp, vector_tmp)
        data_OHS<-rbind(data_OHS, data_OHS_tmp)

}
#add dataset column to the dataframe
data_set<-c(rep("ATL", nrow(OHS2)),rep("ATP", nrow(OHS2)),rep("BCGP", nrow(OHS2)),rep("OHS", nrow(OHS2)),rep("CAG", nrow(OHS2)))
data_OHS2 = data_OHS[-1,]
data_OHS2$dataset<-data_set



BCGP <- as_tibble(BCGP)
colnames(BCGP)
col_order <- c("chr_pos","chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "R2", "maf_ATL", "maf_ATP", "maf_BCGP", "maf_OHS", "maf_CAG")
BCGP2 <- data.frame(BCGP[, col_order])

#need to duplicate the columns into row and add dataset for all the MAF columns to put un a geom_bar
data_BCGP<-NA
for (i in c(11,12,13,14,15)){


        print(i)
        data_tmp<-BCGP2[,1:10]
        vector_tmp<-BCGP2[,i]
        data_BCGP_tmp<-cbind(data_tmp, vector_tmp)
        data_BCGP<-rbind(data_BCGP, data_BCGP_tmp)

}
#add dataset column to the dataframe
data_set<-c(rep("ATL", nrow(BCGP2)),rep("ATP", nrow(BCGP2)),rep("BCGP", nrow(BCGP2)),rep("OHS", nrow(BCGP2)),rep("CAG", nrow(BCGP2)))
data_BCGP2 = data_BCGP[-1,]
data_BCGP2$dataset<-data_set



ATL <- as_tibble(ATL)
colnames(ATL)
col_order <- c("chr_pos","chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "R2", "maf_ATL", "maf_ATP", "maf_BCGP", "maf_OHS", "maf_CAG")
ATL2 <- data.frame(ATL[, col_order])

#need to duplicate the columns into row and add dataset for all the MAF columns to put un a geom_bar
data_ATL<-NA
for (i in c(11,12,13,14,15)){


        print(i)
        data_tmp<-ATL2[,1:10]
        vector_tmp<-ATL2[,i]
        data_ATL_tmp<-cbind(data_tmp, vector_tmp)
        data_ATL<-rbind(data_ATL, data_ATL_tmp)

}
#add dataset column to the dataframe
data_set<-c(rep("ATL", nrow(ATL2)),rep("ATP", nrow(ATL2)),rep("BCGP", nrow(ATL2)),rep("OHS", nrow(ATL2)),rep("CAG", nrow(ATL2)))
data_ATL2 = data_ATL[-1,]
data_ATL2$dataset<-data_set


ATP <- as_tibble(ATP)
colnames(ATP)
col_order <- c("chr_pos","chr", "pos", "SNP", "ref", "alt", "useless", "filter", "info", "R2", "maf_ATL", "maf_ATP", "maf_BCGP", "maf_OHS", "maf_CAG")
ATP2 <- data.frame(ATP[, col_order])

#need to duplicate the columns into row and add dataset for all the MAF columns to put un a geom_bar
data_ATP<-NA
for (i in c(11,12,13,14,15)){


        print(i)
        data_tmp<-ATP2[,1:10]
        vector_tmp<-ATP2[,i]
        data_ATP_tmp<-cbind(data_tmp, vector_tmp)
        data_ATP<-rbind(data_ATP, data_ATP_tmp)

}
#add dataset column to the dataframe
data_set<-c(rep("ATL", nrow(ATP2)),rep("ATP", nrow(ATP2)),rep("BCGP", nrow(ATP2)),rep("OHS", nrow(ATP2)),rep("CAG", nrow(ATP2)))
data_ATP2 = data_ATP[-1,]
data_ATP2$dataset<-data_set





plot_ATP<-ggplot(data_ATP2, aes(x=vector_tmp, color=dataset)) + geom_histogram(fill="white",position="dodge")+ xlab("maf") + ylab("count") +ggtitle("Badly imputed in ATP")  #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/ATP_dist_maf_badly_imputed.png", plot_ATP)


plot_ATL<-ggplot(data_ATL2, aes(x=vector_tmp, color=dataset)) + geom_histogram(fill="white",position="dodge")+ xlab("maf") + ylab("count") +ggtitle("Badly imputed in ATL")  #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/ATL_dist_maf_badly_imputed.png", plot_ATL)



plot_BCGP<-ggplot(data_BCGP2, aes(x=vector_tmp, color=dataset)) + geom_histogram(fill="white",position="dodge")+ xlab("maf") + ylab("count") +ggtitle("Badly imputed in BCGP")  #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/BCGP_dist_maf_badly_imputed.png", plot_BCGP)



plot_OHS<-ggplot(data_OHS2, aes(x=vector_tmp, color=dataset)) + geom_histogram(fill="white",position="dodge")+ xlab("maf") + ylab("count") +ggtitle("Badly imputed in OHS")  #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/OHS_dist_maf_badly_imputed.png", plot_OHS)



plot_CAG<-ggplot(data_CAG2, aes(x=vector_tmp, color=dataset)) + geom_histogram(fill="white",position="dodge")+ xlab("maf") + ylab("count") +ggtitle("Badly imputed in CAG")  #+ facet_wrap(~data_final$groupPerc,scales = "free_y") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=1))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/CAG_dist_maf_badly_imputed.png", plot_CAG)


















#----------------------------------------------Number of SNPs in each cohort-------------------------------------------------------------------------------




good_except_CAG<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/03_08_badly_imputed_in_only_CAG.vcf", header=F)
#good_except_CAG<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/common_to_4_chips_R2_0.3_all_cohorts_except_CAG.ID", header=F)
#good_except_CAG<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/common_to_4_chips_R2_0.8_all_cohorts_except_CAG.ID", header=F)



good_except_OHS<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/03_08_badly_imputed_in_only_OHS.vcf", header=F)
#good_except_OHS<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/common_to_4_chips_R2_0.3_all_cohorts_except_OHS.ID", header=F)
#good_except_OHS<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/common_to_4_chips_R2_0.8_all_cohorts_except_OHS.ID", header=F)



good_except_BCGP<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/03_08_badly_imputed_in_only_BCGP.vcf", header=F)
#good_except_BCGP<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/common_to_4_chips_R2_0.3_all_cohorts_except_BCGP.ID", header=F)
#good_except_BCGP<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/common_to_4_chips_R2_0.8_all_cohorts_except_BCGP.ID", header=F)



good_except_ATL<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/03_08_badly_imputed_in_only_ATL.vcf", header=F)
#good_except_ATL<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/common_to_4_chips_R2_0.3_all_cohorts_except_ATL.ID", header=F)
#good_except_ATL<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/common_to_4_chips_R2_0.8_all_cohorts_except_ATL.ID", header=F)



good_except_ATP<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/03_08_badly_imputed_in_only_ATP.vcf", header=F)
#good_except_ATP<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/common_to_4_chips_R2_0.3_all_cohorts_except_ATP.ID", header=F)
#good_except_ATP<-read.table("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/common_to_4_chips_R2_0.8_all_cohorts_except_ATP.ID", header=F)



total_CAG<-nrow(good_except_CAG)
total_OHS<-nrow(good_except_OHS)
total_BCGP<-nrow(good_except_BCGP)
total_ATL<-nrow(good_except_ATL)
total_ATP<-nrow(good_except_ATP)

only_CAG<-nrow(CAG)
only_OHS<-nrow(OHS)
only_BCGP<-nrow(BCGP)
only_ATL<-nrow(ATL)
only_ATP<-nrow(ATP)

cohort<-c("CAG", "OHS", "BCGP", "ATL", "ATP")
total<-c(total_CAG, total_OHS, total_BCGP, total_ATL, total_ATP)
baddly_imputed<-c(only_CAG, only_OHS, only_OHS, only_BCGP, only_ATL)

data_total<-data.frame(cohort,total)
data_total$type="Well imputed in the 4 other cohorts (R2>=0.8)"
names(data_total)<-c("Cohort", "Sites", "type")
data_badly<-data.frame(cohort, baddly_imputed)
data_badly$type="Badly imputed in this cohort (R2<0.8)"
names(data_badly)<-c("Cohort", "Sites", "type")

data<-rbind(data_total, data_badly)


test<-ggplot(data_total, aes(fill=Cohort, x=Cohort, y=Sites)) + geom_bar(position=position_dodge(width = 0.2), stat="identity") + xlab("CanPath Sub-cohort") + ylab("Number of Sites (SNP+indels)") + guides(fill=guide_legend(title="")) + theme(legend.position = "none") + geom_text(aes(label=Sites), position=position_dodge(width=0.2), size = 10) + theme(text = element_text(size = 40))# + scale_fill_manual(values=c("red", "purple", "blue"))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/stats_badly_imputed_R2_03_08.png",height=12, width=20, test)

#no numbers on top
test<-ggplot(data_total, aes(fill=Cohort, x=Cohort, y=Sites)) + geom_bar(position=position_dodge(width = 0.2), stat="identity") +theme_classic()+ xlab("CanPath sub-cohort") + ylab("Number of variants") + guides(fill=guide_legend(title="")) + theme(text = element_text(size = 40)) + theme(legend.position = "none")# + scale_fill_manual(values=c("red", "purple", "blue"))
ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/no_number_stats_badly_imputed_R2_03_08.png",height=12, width=20, test)
#ggsave("/lustre07/scratch/justinp/topmed_new/CPTP_new/Difficulties_imputation/stats_badly_imputed_R2_03.png",height=12, width=20, test)
















