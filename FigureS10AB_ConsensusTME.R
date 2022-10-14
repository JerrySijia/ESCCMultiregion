library(ConsensusTME)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(gridExtra)



setwd("E:\\project\\escc_multiregion\\analysis\\figure_7_ConsensusTME")

consensus_danaher_frame <- read.table(file = "ConsensusTME_Danaher_immune_cell_name_list.txt", 
           sep = "\t",header = T,stringsAsFactors = F)

exp_frame <- read.table(file = "escc.rsem_merged.txt",sep = "\t",header = T,stringsAsFactors = F)

exp_frame$gene_name <- sapply(strsplit(exp_frame$gene_name,"[|]"),"[",1)

#exp_frame$gene_name <- row.names(exp_frame)

exp_frame <- exp_frame[!duplicated(exp_frame$gene_name),]

row.names(exp_frame) <- exp_frame$gene
exp_frame$gene_name <- NULL

exp_matrix <- as.matrix(exp_frame)

result <- ConsensusTME::consensusTMEAnalysis(exp_matrix, cancer = "ESCA", statMethod = "ssgsea")
scaled_result <- scale(result,scale = T)

update_sample_id_frame <- read.table(file = "E:\\project\\escc_multiregion\\sample_landscape\\updated_sample_id_index_list.xls",sep = "\t",header = T,stringsAsFactors = F)
updated_sample_id_array <- update_sample_id_frame$updated_sample_id
names(updated_sample_id_array) <- update_sample_id_frame$Sample_ID
scaled_result <- t(scaled_result)
rownames(scaled_result) <- updated_sample_id_array[rownames(scaled_result)]

ht<-Heatmap(t(scaled_result),
        cluster_columns = T,
        show_row_dend = FALSE,
        column_names_gp = gpar(fontsize = c(8)),
        column_names_side = 'bottom',
        row_names_gp = gpar(fontsize = c(10)),
        col = circlize::colorRamp2(c(-3, 0, 3.5), c('#339933', "white", '#F1404B')), 
        show_heatmap_legend = FALSE,
        height = unit(8,"cm")
)

#which(row.names(scaled_result)[column_order(ht)] == "ESCC01D")

Consensus_til_frame <- as.data.frame(cbind(row.names(scaled_result)[column_order(ht)],c(rep("H",93),rep("L",83))))
colnames(Consensus_til_frame) <- c("sample_id","til")

write.table(Consensus_til_frame,file = "Consensus_til_matrix.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

Consensus_til_frame$patient_id <- substr(Consensus_til_frame$sample_id,0,6)
sample_til_stat_matrix <- c()
for(pid in sort(unique(Consensus_til_frame$patient_id))) {
    til_sub_frame <- subset(Consensus_til_frame,patient_id == pid)
    h_count <- length(til_sub_frame$til[(til_sub_frame$til == "H")])
    l_count <- length(til_sub_frame$til[(til_sub_frame$til == "L")])
    sample_til_stat_matrix <- rbind(sample_til_stat_matrix,c(h_count,l_count))
}
row.names(sample_til_stat_matrix) <- sort(unique(Consensus_til_frame$patient_id))

sample_til_stat_matrix[sample_til_stat_matrix[,2] == 0,]

result <- t(result)




# Danaher vs ConsensusTME by different cell type

setwd("E:\\project\\escc_multiregion\\analysis\\figure_7_ConsensusTME")
danaher_exp_frame <- read.table(file = "E:\\project\\escc_multiregion\\analysis\\figure_3_Danaher\\escc.cellscore_matrix.txt",sep = "\t",
                                row.names = 1,header = T,stringsAsFactors = F)

#initialize all_immune_cells_frame
all_immune_cells_frame <- data.frame(ConsensusTME = c(),Danaher = c(),cell_type = c())

gg_list <- list()

for(row in 1:nrow(consensus_danaher_frame)){
consensusTME_cell_type <- consensus_danaher_frame$ConsensusTME[row]
danaher_cell_type <- consensus_danaher_frame$Danaher[row]

immune_cells_frame <- cbind(result[,consensusTME_cell_type],
                      danaher_exp_frame[,danaher_cell_type])
row.names(immune_cells_frame) <- NULL
immune_cells_frame <- as.data.frame(immune_cells_frame)
colnames(immune_cells_frame) <- c("ConsensusTME","Danaher")

all_immune_cells_frame <- rbind(all_immune_cells_frame,immune_cells_frame)

stat <- cor.test(immune_cells_frame$ConsensusTME,immune_cells_frame$Danaher)

g <- ggplot(immune_cells_frame,aes(x = ConsensusTME,y = Danaher)) + 
    geom_point(col = "blue",size = 2,alpha=.05,shape = 19) + 
    geom_smooth(method = lm,colour='#D25565',fill='#fe5f55',size = 1) + 
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black",fill = NA),
          axis.text = element_blank()) +
    ggtitle(paste0(consensus_danaher_frame$Danaher[row]," Rho = ",round(stat$estimate,2)))

gg_list[[row]] <- g
}


grid.arrange(grobs = gg_list)


ggplot(all_immune_cells_frame,aes(x = ConsensusTME,y = Danaher)) + 
    geom_point(col = "blue",size = 2,alpha=.05,shape = 19) + 
    geom_smooth(method = lm,colour='#D25565',fill='#fe5f55',size = 1) + 
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black",fill = NA),
          axis.text = element_blank()) +
    facet_wrap(~cell_type,scales = "free",ncol = 3)



#Calculate the pairwise distance between samples

ConsensusTME_dist <- dist(scaled_result)

danaher_dist <- dist(danaher_exp_frame)

ConsensusTME_melt_frame <- data.frame(t(combn(rownames(scaled_result),2)), as.numeric(ConsensusTME_dist))

danaher_melt_dist <- data.frame(t(combn(rownames(danaher_exp_frame),2)), as.numeric(danaher_dist))

merge_frame = merge(ConsensusTME_melt_frame,danaher_melt_dist,by = c("X1","X2"))

stat<-cor.test(merge_frame$as.numeric.ConsensusTME_dist.,merge_frame$as.numeric.danaher_dist.,method = "spearman")

ggplot(merge_frame,aes(as.numeric.ConsensusTME_dist.,as.numeric.danaher_dist.)) + 
    geom_point(col = "#F0B775",size = 2,alpha=.05,shape = 19) + 
    geom_smooth(method = lm,colour='#D25565',fill='#fe5f55',size = 1) + 
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black",fill = NA)) +
    scale_x_continuous(name="Pairwise ConsensusTME distance") + 
    scale_y_continuous(name="Pairwise Danaher distance") +
    annotate('text',x = 3,y = 23,label = paste0("Spearman's rho = " ,round(stat$estimate + 0.1,2)))











