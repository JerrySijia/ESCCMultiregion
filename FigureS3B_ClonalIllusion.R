library(ggplot2)
library(reshape2)
library(plyr)

setwd("E:\\project\\escc_multiregion\\analysis\\figure_13_clonal_illusion")

sample_info_frame <- read.table(file = "E:\\project\\escc_multiregion\\sample_landscape\\updated_sample_id_index_list.xls",sep = "\t",header = T,stringsAsFactors = F)
updated_patient_id <- unique(sample_info_frame$updated_patient_id)
names(updated_patient_id) <- unique(sample_info_frame$patient_id)


clonal_illusion_frame <- read.table(file = "escc.clonal_illusion_stats.xls",sep = "\t",header = TRUE,stringsAsFactors = FALSE)

clonal_illusion_frame$no_illusion_count <- clonal_illusion_frame$subclonal_snv_count - clonal_illusion_frame$illusion_count

print(median(clonal_illusion_frame$illusion_count / clonal_illusion_frame$subclonal_snv_count))
print(range(clonal_illusion_frame$illusion_count/clonal_illusion_frame$subclonal_snv_count))

melt_clonal_illusion_frame <- subset(clonal_illusion_frame,select = -subclonal_snv_count)

melt_clonal_illusion_frame <- melt(melt_clonal_illusion_frame,id = "patient_id")

melt_clonal_illusion_frame$updated_patient_id <- updated_patient_id[melt_clonal_illusion_frame$patient_id]


ggplot(melt_clonal_illusion_frame,aes(x = updated_patient_id,y = value,fill = variable)) + geom_bar(position = position_stack(reverse = TRUE),stat = "identity") + 
    scale_fill_manual(values = c("#2694ab","#C1C1C1"),labels=c("Clonal illusion", "No clonal illusion")) +
    scale_x_discrete(limits = melt_clonal_illusion_frame$updated_patient_id[order(clonal_illusion_frame$subclonal_snv_count,decreasing = TRUE)]) + 
    xlab("") +
    ylab("Subclonal mutations") +
    guides(fill=guide_legend(title=NULL)) +
    theme(axis.text.x = element_text(size = 12,angle = 60, hjust = 1,vjust = 1),
          axis.title.y = element_text(size=16),
          axis.text.y = element_text(size=12),
          panel.grid=element_blank(),
          legend.position = 'top',
          panel.border= element_rect(colour = "black",fill = NA),panel.background = element_blank())

