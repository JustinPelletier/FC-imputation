library(ggplot2)

data_CAG<-read.table("nb_lines_CAG_miss005.txt", header=T)
#data_CAG<-read.table("nb_lines_CAG.txt", header=T)
#data_CAG<-read.table("modified_nb_lines_CAG.txt", header=T)
#data_CAG<-read.table("modified_nb_lines_CAG.txt", header=T)
data_ini_CAG<-read.table("nb_lines_initial_CAG_miss005.txt", header=T)
#data_ini_CAG<-read.table("nb_lines_initial_CAG.txt", header=T)
#data_ini_CAG<-read.table("modified_nb_lines_initial_CAG.txt", header=T)
#data_ini_CAG<-read.table("modified_nb_lines_initial_CAG.txt", header=T)
data_CAG_old<-read.table("final_nb_lines_CAG_old_miss005.txt", header=T)
#data_CAG_old<-read.table("final_nb_lines_CAG_old.txt", header=T)
#data_CAG_old<-read.table("modified_final_nb_lines_CAG_old.txt", header=T)
#data_CAG_old<-read.table("modified_final_nb_lines_CAG_old.txt", header=T)


#remove column
data_CAG$ini_SNP<-NULL
data_CAG_old$ini_SNP<-NULL


data_CAG$r2<-as.character(data_CAG$r2)
data_ini_CAG$r2<-as.character(data_ini_CAG$r2)
data_CAG_old$r2<-as.character(data_CAG_old$r2)
#data_CAG$chip <- factor(data_CAG$chip, levels = c("",""))

data_CAG$ref<-"TOPMed imputation"
data_ini_CAG$ref<-"Genotyped/QCed"
data_CAG_old$ref<-"HRC imputation"


test_CAG<-rbind(data_CAG, data_CAG_old, data_ini_CAG)
#for liftover figure
#test_CAG<-read.table("nb_lines_CAG_before_after_lift_over.txt", header=T)

png(filename = "CAG_lift_over_nb_SNP.png", width = 5000, height = 2500)
p1<-ggplot(test_CAG, aes(fill=assembly, x=reorder(chip, SNP), y=SNP)) + geom_bar(position="identity", stat="identity") + xlab("Puce de génotypage") + ylab("Nombre de SNPs") + guides(fill=guide_legend(title="Cohorte de référence")) + facet_wrap(~test_CAG$r2) + theme(text = element_text(size = 78))
print(p1)

dev.off()



png(filename = "test_modified_CAG_nb_SNP_initial.png", width = 5000, height = 2500)
p4<-ggplot(test_CAG, aes(fill=ref, x=reorder(chip, SNP), y=SNP)) + geom_bar(position="identity", stat="identity",alpha=0.7) + xlab("Puce de génotypage") + ylab("Nombre de SNPs") + guides(fill=guide_legend(title="Cohorte de référence")) + facet_wrap(~test_CAG$r2) + theme(text = element_text(size = 78))
print(p4)

dev.off()


#write.table(test_CAG, "CAG_nb_line_from_script_r.txt", sep = "	", row.names = F, col.names = TRUE)
test_CAG_superpose<-read.table("CAG_nb_line_from_script_r.txt", header=T)

#test superposition
png(filename = "test_superposition_CAG_nb_SNP_initial.png", width = 5000, height = 2500)
p5<-ggplot(test_CAG_superpose, aes(fill=ref, x=reorder(chip, SNP), y=SNP)) + geom_bar(position="identity", stat="identity",alpha=0.7) + xlab("Puce de génotypage") + ylab("Nombre de SNPs") + guides(fill=guide_legend(title="Cohorte de référence")) + theme(text = element_text(size = 78))
print(p5)

dev.off()




#-----------------------------FINAL PLOTS------------------------------------------



p1<-ggplot(data_CAG, aes(fill=r2, x=reorder(chip, SNP), y=SNP)) + geom_bar(position="dodge", stat="identity") + xlab("Genotyping array") + ylab("Number of sites") + guides(fill=guide_legend(title="R squared")) #+ scale_fill_viridis(discrete = T)
ggsave("test_CAG_nb_SNP.png", p1)



p2<-ggplot(data_CPTP, aes(fill=r2, x=reorder(chip, SNP), y=SNP)) + geom_bar(position="dodge", stat="identity") + xlab("Genotyping chip") + ylab("Number of SNPs") + guides(fill=guide_legend(title="R squared")) #+ scale_fill_viridis(discrete = T)
ggsave("test_CPTP_nb_SNP.png", p2)



p3<-ggplot(test_CAG, aes(fill=r2, x=reorder(chip, SNP), y=SNP)) + geom_bar(position="dodge", stat="identity") + xlab("Genotyping microarray") + ylab("Number of sites") + guides(fill=guide_legend(title="R squared")) + facet_wrap(~test_CAG$ref)
ggsave("test_CAG_old_nb_SNP.png", p3)


test_CAG$r2<-factor(test_CAG$r2,levels=c("R2<0.3", "0.3<=R2<0.8", "R2>=0.8"))
test_CAG$chip<-factor(test_CAG$chip,levels=c("GSA_760","AXIOM", "OMNI", "GSA_4224", "GSA_5300", "GSA_17K" ))


#png(filename = "test_CAG_old_nb_SNP.png", width = 4000, height = 2500)
png(filename = "v2_test_CAG_old_nb_SNP.png", width = 5000, height = 2500)
#p4<-ggplot(test_CAG, aes(fill=ref, x=reorder(chip, SNP), y=SNP)) + geom_bar(position="dodge", stat="identity") + xlab("Genotyping chip") + ylab("Number of SNPs") + guides(fill=guide_legend(title="Reference panel")) + facet_wrap(~test_CAG$r2) + theme(text = element_text(size = 78)) 
p4<-ggplot(test_CAG, aes(fill=ref, x=chip, y=SNP)) + geom_bar(position="dodge", stat="identity")+theme_classic() + theme(axis.text.x = element_text(angle = 35, vjust = 0.9, hjust=1)) + xlab("Genotyping Array") + ylab("Number of sites") + guides(fill=guide_legend(title="")) + facet_wrap(~test_CAG$r2,scales = "free_y") + theme(text = element_text(size = 78)) 
print(p4)
#ggsave("test_CAG_old_nb_SNP.png",width = NA,  height = NA, p4)
dev.off()








#percentage of increase average between HRC and TOPmed
chips<-unique(test_CAG$chip)
average_increase_08<-c()
average_increase_03<-c()
average_hrc_03<-c()
average_hrc_08<-c()


for (i in 1:length(chips)){

	print(chips[i])
		
	increase_03<-(as.numeric(test_CAG[which(test_CAG$chip==chips[i] & test_CAG$r2=="R2>=0.3"& test_CAG$ref=="TOPMed imputation"),]$SNP)-as.numeric(test_CAG[which(test_CAG$chip==chips[i] & test_CAG$r2=="R2>=0.3" & test_CAG$ref=="HRC imputation"),]$SNP))
	increase_08<-(as.numeric(test_CAG[which(test_CAG$chip==chips[i] & test_CAG$r2=="R2>=0.8"& test_CAG$ref=="TOPMed imputation"),]$SNP)-as.numeric(test_CAG[which(test_CAG$chip==chips[i] & test_CAG$r2=="R2>=0.8" & test_CAG$ref=="HRC imputation"),]$SNP))


	hrc_03<-as.numeric(test_CAG[which(test_CAG$chip==chips[i] & test_CAG$r2=="R2>=0.3"& test_CAG$ref=="HRC imputation"),]$SNP)
	hrc_08<-as.numeric(test_CAG[which(test_CAG$chip==chips[i] & test_CAG$r2=="R2>=0.8"& test_CAG$ref=="HRC imputation"),]$SNP)
	
	average_hrc_03<-c(average_hrc_03, hrc_03)
	average_hrc_08<-c(average_hrc_08, hrc_08)


	average_increase_08<-c(average_increase_08, increase_08)
	average_increase_03<-c(average_increase_03, increase_03)
	


}

mean_average_increase_03<-mean(average_increase_03)
mean_average_increase_08<-mean(average_increase_08)

mean_hrc_03<-mean(average_hrc_03)
mean_hrc_08<-mean(average_hrc_08)

perc_increase_topmed_03<-(mean_average_increase_03/mean_hrc_03)
perc_increase_topmed_08<-(mean_average_increase_08/mean_hrc_08)







