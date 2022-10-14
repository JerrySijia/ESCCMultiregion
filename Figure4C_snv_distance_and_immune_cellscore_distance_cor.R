
library(ggplot2)

setwd("E:\\project\\escc_multiregion\\immune\\Danaher")
cellscores_frame <- read.table(file = "escc.cellscore_matrix.txt",sep = "\t",header = TRUE)
snv_frame <- read.table(file = "escc_snv_matrix.txt",sep = "\t",header = TRUE)


cellscores_dist <- dist(cellscores_frame,method = "euclidean")
snv_dist <- dist(t(snv_frame),method = "euclidean")

cellscores_melt_frame <- data.frame(t(combn(rownames(cellscores_frame),2)), as.numeric(cellscores_dist))
snv_melt_frame <- data.frame(t(combn(colnames(snv_frame),2)), as.numeric(snv_dist))

snv_cellscore_merge_frame = merge(cellscores_melt_frame,snv_melt_frame,by = c("X1","X2"))

snv_cellscore_merge_frame = snv_cellscore_merge_frame[(substr(merge_frame$X1,0,6) == substr(merge_frame$X2,0,6)),]

colnames(snv_cellscore_merge_frame) <- c("sample_id_A","sample_id_B","cellscore_distance","snv_distance")

#write.table(snv_cellscore_merge_frame,file = "escc.snv_cellscore_matrix.tsv",sep = "\t",col.names = T,row.names = F,quote = F)

stat<-cor.test(snv_cellscore_merge_frame$cellscore_distance,snv_cellscore_merge_frame$snv_distance,method = "spearman")



ggplot(snv_cellscore_merge_frame,aes(x = cellscore_distance,y = snv_distance)) + 
    geom_point(col = "steelblue",size = 2,alpha=.5,shape = 19) + 
    theme(panel.grid=element_blank(),panel.background = element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour = "black")) +
    scale_x_continuous(name="Pairwise immune distance") + 
    scale_y_continuous(name="Pairwise genomic distance") +
    annotate('text',x = 2,y = 23,label = paste0("Spearman's rho = " ,signif(stat$estimate,2)))
    








