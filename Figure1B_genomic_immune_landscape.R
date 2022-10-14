library(ggplot2)
library(reshape)
library(gridGraphics)
library(ComplexHeatmap)
library(circlize)

panel_len <- 60456963


#==============TMB======================
snv_stat_frame <- read.table(file="escc.snv_stats_by_sample.xls", sep = "\t",header = T)


snv_stat_frame$tmb <- snv_stat_frame$total_snv_num / (panel_len / 1000000)

snv_stat_frame$patient_id <- substr(snv_stat_frame$sample_id,0,6)


tmb_group <- aggregate(x = snv_stat_frame$tmb,by = list(snv_stat_frame$patient_id), FUN = median)

patient_order_array <- tmb_group$Group.1[order(tmb_group$x)]


tmb_gp_up <- ggplot(snv_stat_frame,aes(x = patient_id,y = tmb)) + 
    geom_boxplot(width = 0.6,size = 0.3,color = "#5bd1d7") + 
    geom_point(position = position_jitter(0.1),color = "#2AACE3") + 
    scale_x_discrete(limits = patient_order_array) +
    xlab("") +
    ylab("") +
    scale_y_continuous(expand = c(0, 0),limits = c(10,20),breaks = c(10,15,20)) +
    theme(axis.line.y = element_line(colour = "black"),
          plot.margin = margin(0.2,1.3,0,0.9, "cm"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())


tmb_gp_down <- ggplot(snv_stat_frame,aes(x = patient_id,y = tmb)) + 
    geom_boxplot(width = 0.6,size = 0.3,color = "#5bd1d7") + 
    geom_point(position = position_jitter(0.1),color = "#2AACE3") + 
    scale_x_discrete(limits = patient_order_array) +
    xlab("") +
    ylab("TMB") +
    scale_y_continuous(expand = c(0, 0),limits = c(0,6)) +
    theme(axis.line.y = element_line(colour = "black"),
          plot.margin = margin(0.2,1.3,0,1.1, "cm"),
          axis.line.x = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())


#===========snv proportion stats=================
snv_stat_frame <- read.table(file = "escc.snv_stats_by_patients.xls",sep = "\t",header = T,stringsAsFactors = F)

snv_stat_frame$branch_proportion <- 1 - snv_stat_frame$root_proportion

snv_stat_frame <- subset(snv_stat_frame,select = c("patient_id","root_proportion","branch_proportion"))
melt_snv_stat_frame <- melt(snv_stat_frame,id = c("patient_id"))


snv_prop_gp <- ggplot(melt_snv_stat_frame,aes(x = patient_id,y = value,fill = variable)) +
    geom_bar(stat = "identity",position = position_stack(reverse = T)) +
    xlab("") +
    ylab("Subclonal Mutation\nPercentage") +
    scale_y_continuous(expand = c(0, 0),limits = c(0,1)) +
    scale_fill_manual(values = c("root_proportion" = "#017890","branch_proportion" = "#E95D22")) +
    scale_x_discrete(limits = patient_order_array) +
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(colour = "black",size = 0.5),
          axis.ticks.x = element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.position = "none",
          plot.margin = margin(0,1.3,0,0.2, "cm"),
          axis.line=element_line(size=1,colour = "black"),panel.background = element_blank())




#===========cna proportion stats==================
escc_cna_stat_by_patient_frame <- read.table(file = "escc.cn_stats_by_sample.xls",sep = "\t",header = TRUE,stringsAsFactors = FALSE)

escc_cna_stat_by_patient_frame$patient_id <- substr(escc_cna_stat_by_patient_frame$sample_id,0,6)

escc_cna_stat_by_patient_frame <- escc_cna_stat_by_patient_frame[c("patient_id","clone_cn_length","patient_total_cn_length")]
escc_cna_stat_by_patient_frame <- unique(escc_cna_stat_by_patient_frame)

escc_cna_stat_by_patient_frame$subclone_cna_proportion <- 1 - escc_cna_stat_by_patient_frame$clone_cn_length / escc_cna_stat_by_patient_frame$patient_total_cn_length
escc_cna_stat_by_patient_frame$clone_cna_proportion <- escc_cna_stat_by_patient_frame$clone_cn_length / escc_cna_stat_by_patient_frame$patient_total_cn_length


escc_cna_stat_by_patient_frame <- subset(escc_cna_stat_by_patient_frame,select = c(patient_id,clone_cna_proportion,subclone_cna_proportion))

melt_escc_cna_stat_by_patient_frame <- melt(escc_cna_stat_by_patient_frame,id = "patient_id")

cna_proportion_gp <- ggplot(melt_escc_cna_stat_by_patient_frame,aes(x = patient_id,y = value,fill = variable)) + 
    geom_bar(position = position_stack(reverse = T),stat = "identity") +
    scale_fill_manual(values = c("clone_cna_proportion" = "#017890","subclone_cna_proportion" = "#E95D22")) +
    scale_x_discrete(limits = patient_order_array) + 
    xlab("") +
    ylab("Subclonal Copy Number\nPercentage") +
    scale_y_continuous(expand = c(0, 0),limits = c(0,1)) +
    guides(fill=guide_legend(title=NULL)) +
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(colour = "black",size = 0.5),
          axis.ticks.x = element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.position = "none",
          plot.margin = margin(0,1.3,0,0.2, "cm"),
          axis.line=element_line(size=1,colour = "black"),panel.background = element_blank())



#===============mutation signature======================
setwd("E:\\project\\escc_multiregion\\dna_analysis\\mutation_signature")

clonal_sig_burden_frame = read.table(file = "ESCC.clone_signatures_matrix.txt",sep = "\t",header = T,stringsAsFactors = FALSE)
clonal_sig_burden_frame$type <- "clonal"
clonal_sig_burden_frame <- subset(clonal_sig_burden_frame,Signatures %in% c("Signature.1","Signature_2_13","Signature.3","Signature.4","Signature.5","Signature_6_15_20_26","Signature_unclassified"))


subclonal_sig_burden_frame = read.table(file = "ESCC.subclone_signatures_matrix.txt",sep = "\t",header = T,stringsAsFactors = FALSE)
subclonal_sig_burden_frame$type <- "subclonal"
subclonal_sig_burden_frame <- subset(subclonal_sig_burden_frame,Signatures %in% c("Signature.1","Signature_2_13","Signature.3","Signature.4","Signature.5","Signature_6_15_20_26","Signature_unclassified"))


clonal_subclonal_sig_burden_frame <- rbind(clonal_sig_burden_frame,subclonal_sig_burden_frame)

clonal_subclonal_sig_burden_frame <- subset(clonal_subclonal_sig_burden_frame,Signatures != "Signature_unclassified")

clonal_subclonal_sig_burden_frame <- subset(clonal_subclonal_sig_burden_frame,Signatures %in% c("Signature.1","Signature_2_13","Signature.3","Signature.4","Signature.5","Signature_6_15_20_26","Signature_unclassified"))

#clonal_subclonal_sig_burden_frame$Signature_weight[clonal_subclonal_sig_burden_frame$Signature_weight == 0] <- 0.01


signature_color_array <- c("Signature.1" = "#E5446D",
                           "Signature_2_13" = "#59c8df",
                           "Signature.3" = "#2b9464",
                           "Signature.4" = "#f5df65",
                           "Signature.5" = "#ff7473",
                           "Signature_6_15_20_26" = "#FFB310",
                           "Signature_unclassified" = "grey")

all_sig_burden_frame = read.table(file = "ESCC.all_signatures_mut_burden_matrix.txt",sep = "\t",header = T,stringsAsFactors = FALSE)

all_sig_burden_frame <- subset(all_sig_burden_frame,Signatures %in% c("Signature.1","Signature_2_13","Signature.3","Signature.4","Signature.5","Signature_6_15_20_26","Signature_unclassified"))

clonal_gt <- ggplot(clonal_sig_burden_frame,aes(x = sample_id,y = Signature_weight,fill = Signatures)) + 
    guides(fill=guide_legend(title=NULL)) +
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(colour = "black",size = 0.5),
          axis.ticks.x = element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.position = "none",
          plot.margin = margin(0,1.3,0,0.2, "cm"),
          axis.line=element_line(size=1,colour = "black"),panel.background = element_blank()) + 
    xlab("") +
    ylab("Clonal\nWeight") +
    scale_y_continuous(expand = c(0, 0),limits = c(0,1)) +
    scale_x_discrete(limits = patient_order_array) +
    scale_fill_manual(values = signature_color_array) +
    geom_bar(stat="identity",width = 0.9,position = position_stack(reverse = TRUE))



subclonal_gt <- ggplot(subclonal_sig_burden_frame,aes(x = sample_id,y = Signature_weight,fill = Signatures)) + 
    guides(fill=guide_legend(title=NULL)) +
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(colour = "black",size = 0.5),
          axis.ticks.x = element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank(),
          legend.position = "none",
          plot.margin = margin(0,1.3,0,0.2, "cm"),
          axis.line=element_line(size=1,colour = "black"),panel.background = element_blank()) + 
    xlab("") +
    ylab("SubClonal\nWeight") +
    scale_y_continuous(expand = c(0, 0),limits = c(0,1)) +
    scale_x_discrete(limits = patient_order_array) +
    scale_fill_manual(values = signature_color_array) +
    geom_bar(stat="identity",width = 0.9,position = position_stack(reverse = TRUE))



#==================clinical information============================================
col_array = c()
label_array = c()

grid_NA_color = "grey"
grid_border_color = "white"

#GD color
col_array <- c(col_array,c("#CC2E2C","#0089B3","#FFAA3D"))
label_array <- c(label_array,c("GD","NGD","subGD"))

#evolution color
col_array <- c(col_array,c("#6FBC45","#2A4BA1"))
label_array <- c(label_array,c("Yes","No"))


#TIL color
col_array <- c(col_array,c("#F00000","#FF9933","#339933","#BEBEBE"))
label_array <- c(label_array,c("H","M","L","UNKNOWN"))

#HLALOH color
col_array <- c(col_array,c("#614AD3","#d9d9f3"))
label_array <- c(label_array,c("Y","N"))



setwd("E:\\project\\escc_multiregion\\sample_info")
patient_info_frame <- read.table(file = "all_escc_patient_info_list.xls", sep = "\t", header = TRUE,stringsAsFactors = FALSE)
patient_info_frame <- subset(patient_info_frame,select = c("patient_id","GD","EvoPattern","TIL","HLALOH"))

rownames(patient_info_frame) <- patient_info_frame$patient_id
patient_info_frame <- subset(patient_info_frame,select = c(-patient_id))
patient_info_frame <- patient_info_frame[patient_order_array,]

update_patient_id_frame <- read.table(file = "E:\\project\\escc_multiregion\\sample_info\\updated_sample_id_index_list.xls",sep = "\t",header = T,stringsAsFactors = F)
updated_patient_id_array <- unique(update_patient_id_frame$updated_patient_id)
names(updated_patient_id_array) <- unique(update_patient_id_frame$patient_id)
rownames(patient_info_frame) <- updated_patient_id_array[rownames(patient_info_frame)]
patient_info_frame <- t(patient_info_frame)


patient_clinical_info_ht <- Heatmap(patient_info_frame,
                                    col = structure(col_array,names = label_array),
                                    rect_gp = gpar(col = grid_border_color,lwd = 2),
                                    show_column_dend = FALSE,
                                    show_column_names = T,
                                    column_names_rot = 60,
                                    row_names_gp = gpar(fontsize = 10),
                                    column_names_gp = gpar(fontsize = 8),
                                    show_row_names = T,
                                    row_names_side = "left",
                                    cluster_rows = FALSE,
                                    cluster_columns = FALSE,
                                    show_heatmap_legend = FALSE,
                                    height = unit(4,"cm"),
                                    na_col = "grey")

#=========================================merge==============================================
setwd("E:\\project\\escc_multiregion\\dna_analysis\\2020_genomic_immune_landscape")


tmb_show_up_gg = grid.grabExpr({
    tmb_gp_up
    print(tmb_gp_up)
})

tmb_show_down_gg = grid.grabExpr({
    tmb_gp_down
    print(tmb_gp_down)
})

snv_prop_show_gg = grid.grabExpr({
    snv_prop_gp
    print(snv_prop_gp)
})

cna_prop_show_gg = grid.grabExpr({
    cna_proportion_gp
    print(cna_proportion_gp)
})

clone_sig_show_gg = grid.grabExpr({
    clonal_gt
    print(clonal_gt)
})

subclone_sig_show_gg = grid.grabExpr({
    subclonal_gt
    print(subclonal_gt)
})

clinical_info_heatmap = grid.grabExpr(draw(patient_clinical_info_ht,padding = unit(c(0, 0.5, 0, 13.5), "mm")))


grid.newpage()

pushViewport(viewport(layout = grid.layout(nr = 13, nc = 1)))

pushViewport(viewport(layout.pos.row = c(1), layout.pos.col = 1))
grid.draw(tmb_show_up_gg)
popViewport()
pushViewport(viewport(layout.pos.row = c(2,3), layout.pos.col = 1))
grid.draw(tmb_show_down_gg)
popViewport()

pushViewport(viewport(layout.pos.row = c(4,5), layout.pos.col = 1))
grid.draw(snv_prop_show_gg)
popViewport()

pushViewport(viewport(layout.pos.row = c(6,7), layout.pos.col = 1))
grid.draw(cna_prop_show_gg)
popViewport()

pushViewport(viewport(layout.pos.row = c(8), layout.pos.col = 1))
grid.draw(clone_sig_show_gg)
popViewport()

pushViewport(viewport(layout.pos.row = c(9), layout.pos.col = 1))
grid.draw(subclone_sig_show_gg)
popViewport()


pushViewport(viewport(layout.pos.row = c(10,11,12,13), layout.pos.col = 1))
grid.draw(clinical_info_heatmap)
popViewport()
















#legend


clonal_subclonal_lgd = Legend(labels = c("clonal","subclonal"), 
                            title = "Clonal/Subclonal SNV/SCNA", 
                            legend_gp = gpar(fill = c("#017890",
                                                      "#E95D22")))




mutation_signature_lgd = Legend(labels = c("Age","APOBEC","BRCA","Smoking","Unknown","DMMR","Signature_unclassified"), 
                 title = "Mutation Signatures", 
                 legend_gp = gpar(fill = c("#E5446D",
                                          "#59c8df",
                                          "#2b9464",
                                          "#f5df65",
                                          "#2b9464",
                                          "#ff7473",
                                          "#FFB310",
                                          "grey")))

gd_legend = Legend(labels = c("GD","NGD","Subclonal GD"), 
                   title = "Genome doubling status", 
                   legend_gp = gpar(fill = c("#CC2E2C",
                                             "#0089B3","#FFAA3D")))


#evolution color

evo_legend = Legend(labels = c("Linear","Branched"), 
                   title = "Evolution Pattern", 
                   legend_gp = gpar(fill = c("#6FBC45",
                                             "#2A4BA1")))



#TIL color

til_legend = Legend(labels = c("H","M","L","UNKNOWN"), 
                    title = "TIL status", 
                    legend_gp = gpar(fill = c("#F00000",
                                              "#FF9933",
                                              "#339933",
                                              "#BEBEBE")))


#HLALOH color

hlaloh_legend = Legend(labels = c("Any HLALOH","No HLALOH"), 
                    title = "HLALOH status", 
                    legend_gp = gpar(fill = c("#614AD3",
                                              "#d9d9f3")))





pd = packLegend(clonal_subclonal_lgd,
                mutation_signature_lgd,
                gd_legend,
                evo_legend,
                til_legend,
                hlaloh_legend,
                direction = "horizontal",
                max_width = unit(10, "cm"),
                column_gap = unit(1, "cm"))

pushViewport(viewport(width = 0.8, height = 0.8))
grid.rect()
draw(pd)
popViewport()


gd_lgd = Legend(labels = c("GD","no GD"), title = "Genome Doubling status", legend_gp = gpar(fill = c("#CC2E2C",
                                                                                                      "#0089B3")))
pd = packLegend(gd_lgd,
                direction = "horizontal",
                max_width = unit(10, "cm"),
                column_gap = unit(1, "cm"))

pushViewport(viewport(width = 0.8, height = 0.8))
grid.rect()
draw(pd)
popViewport()















