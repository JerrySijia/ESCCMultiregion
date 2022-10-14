library(ggplot2)
library(reshape)

setwd("E:\\project\\escc_multiregion\\analysis\\figure_20_lohhla_summary")

neo_lohhla_frame <- read.table(file = "escc.neo_lohhla_stat.tsv", sep = "\t", header = T,stringsAsFactors = F)

patient_order_array <- unique(neo_lohhla_frame$patient_id)[order(neo_lohhla_frame$num[neo_lohhla_frame$type == 'clone_loss'],decreasing = T)]
cast_neo_lohhla_frame <- cast(neo_lohhla_frame,patient_id~type)
write.table(cast_neo_lohhla_frame,file = "neo_lohhla_frame.tsv",sep = "\t",row.names = F,col.names = T,quote = F)

ggplot(neo_lohhla_frame,aes(x = patient_id,y = num,fill = type)) + 
    geom_bar(stat = "identity",position = position_stack(reverse = T)) +
    theme(axis.text.x = element_text(size = 12,angle = 60, hjust = 1,vjust = 1),
          axis.title.y = element_text(size=16),
          axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=12),
          panel.grid=element_blank(),
          legend.position = 'top',
          panel.border= element_rect(colour = "black",fill = NA),panel.background = element_blank()) +
    scale_x_discrete(limits = patient_order_array) +
    scale_fill_manual(values = c("clone_loss" = "#388AB7","subclone_loss" = "#F0863E"),guide=FALSE) +
    xlab("") +
    ylab("Number of predicted binders to lost HLA allele")
    



















