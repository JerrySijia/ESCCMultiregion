library(entropy)
library(ggplot2)
library(reshape2)

setwd("E:\\project\\escc_multiregion\\analysis\\figure_12_pyclone_count_barplot")
pyclone_frame <- read.table(file = "ESCC.clone_count_statx.tsv",sep = "\t",header = T)

total_clone_count <- sum(pyclone_frame$clone_count)

median(pyclone_frame$clone_count)

range(pyclone_frame$clone_count)



#================entropy========================
all_clone_prev_frame = read.table(file = "escc.clonal_prev.tsv",sep ="\t",header = T,stringsAsFactors = F)
all_clone_prev_frame <- all_clone_prev_frame[c("timepoint","clonal_prev")]
sample_entropy_frame <- aggregate(x = all_clone_prev_frame$clonal_prev, by = list(all_clone_prev_frame$timepoint),FUN = entropy.empirical)

colnames(sample_entropy_frame) <- c("sample_id","entropy")
write.table(sample_entropy_frame,file = "escc.entropy_matrix.tsv",sep = "\t",row.names = F,col.names = T,quote = F)





#=============pyclone clone count barplot=================

sample_list_frame <- read.table(file = "E:\\project\\escc_multiregion\\all_organized_tables\\escc_20200811_by_patient.all_info_list.xls",
           sep ="\t",
           header = T,
           stringsAsFactors = F)

update_sample_id_array <- sample_list_frame$updated_sample_id
names(update_sample_id_array) <- sample_list_frame$patient_id

pyclone_count_frame <- read.table(file = "escc.pyclone_stat.tsv",sep = "\t",header = T,stringsAsFactors = F)

pyclone_count_frame$updated_sample_id <- update_sample_id_array[pyclone_count_frame$sample_id]
pyclone_count_frame$sample_id <- NULL

pyclone_count_frame$subclone_count <- pyclone_count_frame$total_clone - pyclone_count_frame$subclone_on_single_branch_count
ordered_sample_id_array <- (pyclone_count_frame[order(pyclone_count_frame$total_clone_count,decreasing = T),])$updated_sample_id
pyclone_count_frame$total_clone_count <- NULL
melt_pyclone_count_frame <- melt(pyclone_count_frame,id=c("updated_sample_id"))

ggplot(melt_pyclone_count_frame,
       aes(x = updated_sample_id,y = value,fill = factor(variable))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(size = 12,angle = 60, hjust = 1,vjust = 1),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid=element_blank(),
        legend.position = 'top',
        panel.border= element_rect(colour = "black",fill = NA),panel.background = element_blank()) +
  guides(fill=guide_legend(title=NULL)) +
  xlab("") +
  ylab("Subclone Number") +
  scale_x_discrete(limits = ordered_sample_id_array) +
  scale_fill_manual(values = c("subclone_count" = "#F8937E","subclone_on_single_branch_count" = "#F9DEC9"),
                      labels=c("Subclone","Subclone on single branch"))
















