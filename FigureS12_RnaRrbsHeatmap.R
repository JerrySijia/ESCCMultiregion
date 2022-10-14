
library(ComplexHeatmap)
library(stringr)
library(circlize)
library(ggplot2)



setwd("E:\\project\\escc_multiregion\\tri-omics\\triomics_cluster_heatmap")

sample_order_frame <- read.table(file = "sample_order.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)

#sample_order_frame <- subset(sample_order_frame,sample_id %in% colnames(tpm_frame))

sample_order_frame$patient_id <- substr(sample_order_frame$sample_id,1,6)

#RNA

#row.names(tpm_frame) <- tpm_frame$geneName
tpm_frame <- read.table(file = "escc_all.TPM_std_all_matrix.txt",sep = "\t",header = TRUE)
tpm_frame <- tpm_frame[sample_order_frame$sample_id]
tpm_frame <- scale(tpm_frame,center = TRUE,scale = TRUE)

#hc <- hclust(dist(tpm_frame[,1:61],method = "euclidean"),method = "ward.D")
#group_A_order<-dendsort(hc)$order



rna_ht <- Heatmap(tpm_frame[5000:5050,],
                show_row_names = FALSE,
                show_column_names = FALSE,
                border = TRUE,
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                col = colorRamp2(c(-4,0, 5), c("#4A1F67", "#3A6589", "#ECE129")),
                #row_dend_reorder = TRUE,
                show_heatmap_legend = FALSE,
                column_dend_reorder = FALSE,
                #column_split = factor(sample_order_frame$patient_id,levels = unique(sample_order_frame$patient_id)),
                show_row_dend = FALSE,
                cluster_column_slices = FALSE,
                #column_gap = unit(c(0.5), "mm"),
                column_title = NULL)


rrbs_frame <- read.table(file = "escc.all_methlevel_bin_std_frame.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
rrbs_frame <- rrbs_frame / 100
rrbs_frame <- rrbs_frame[,sample_order_frame$sample_id]

#rrbs_frame <- scale(rrbs_frame,center = TRUE,scale = TRUE)


rrbs_ht <- Heatmap(rrbs_frame[1000:1050,],
            show_row_names = FALSE,
            border = TRUE,
            col = colorRamp2(c(0, 0.5, 1.1), c("#3A54A4", "white", "#EE2123")),
            cluster_columns = FALSE,
            cluster_rows = T,
            row_names_gp = gpar(fontsize = 5),
            column_dend_reorder = FALSE,
            column_names_gp = gpar(fontsize = 5),
            #column_order = colnames(tpm_frame),
            show_heatmap_legend = FALSE,
            row_dend_reorder = TRUE,
            show_column_names = FALSE,
            #column_gap = unit(c(0.5), "mm"),
            #column_split = factor(sample_order_frame$patient_id,levels = unique(sample_order_frame$patient_id)),
            column_title = NULL,
           # border = TRUE,
            show_row_dend = FALSE)


#====================================

leng<-23
ploidy <- 2

patient_margin = 0.0
sample_margin = 0

sample_band_width = 0.8

plot_cnv_band_width = 1

### chromosome length (1-22), centromere position
chrLen<-c(249250621,
          243199373,
          198022430,
          191154276,
          180915260,
          171115067,
          159138663,
          146364022,
          141213431,
          135534747,
          135006516,
          133851895,
          115169878,
          107349540,
          102531392,
          90354753,
          81195210,
          78077248,
          59128983,
          63025520,
          48129895,
          51304566,
          0)
names(chrLen)<-c(1:22)
cumChrLen<-c(0,cumsum(chrLen))[-length(chrLen)-1]
names(cumChrLen)<-c(1:22)
chrpos <- cumChrLen[1:leng-1]+chrLen[1:leng-1]/2


#cnv part
cnv <- read.table("mescc.cnv.txt",sep="\t",header = TRUE)
cnv <- subset(cnv,end - start > 1000000)
cnv$patient_id <- substr(cnv$sample,1,6)


sample_clust <- read.table(file = "sample_order.txt",sep = "\t",header = TRUE)
cnv <- subset(cnv,cnv$sample %in% sample_clust$sample_id)
cnv <- merge(sample_clust,cnv,by.x = "sample_id",by.y = "sample",sort = FALSE)

cnv <- subset(cnv,cn != "neutral")

patient_id_array <- as.array(unique(as.character(cnv$patient_id)))

sample_id_array <- as.array(unique(as.character(cnv$sample)))

#plot
plot_height = length(sample_id_array) * 0.5 + (length(patient_id_array) - 1) * 0.2
plot_width = 22 * 2 * 0.5

sample_number <- length(sample_id_array)
patient_number <- length(patient_id_array)

cols<-c("#EE2223",
        "#FFB5C5",
        "#F4F0EF",
        "#1874CD",
        "orange")

# 4: "amplification"; 3: "gain"; 2: "neutral"; 1: "loss";-1: "upde"
names(cols)<-c('amplification','gain','neutral','loss','neutral_loh')
if(FALSE){
    par(mar=c(2,1,2,3))
    plot(c(cumChrLen[leng],cumChrLen[leng]),
         c(0,(sample_number - patient_number) * sample_margin + (patient_number - 1) * patient_margin + sample_number),
         col="white",
         xlim=c(0,cumChrLen[leng]),
         xaxt="n",
         yaxt="n",
         xlab="",
         ylab="",
         bty="n")
}
#rect(0,0,cumChrLen[leng],(sample_number - patient_number) * sample_margin + (patient_number - 1) * patient_margin + sample_number,lwd = 2)

plotted_sample_num = 0
plotted_patient_num = 0

cnv_x1 <- c()
cnv_x2 <- c()
cnv_y1 <- c()
cnv_y2 <- c()
cnv_col <- c()

all_patient_frame <- data.frame()

for(i in 1:patient_number){
    patient_frame <- subset(cnv,patient_id == patient_id_array[i])
    patient_sample_id_array <- as.array(unique(as.character(patient_frame$sample_id)))
    
    cnv_y1 <- c(cnv_y1,0 - 1000000)
    cnv_y2 <- c(cnv_y2,cumChrLen[leng] + 100000)
    cnv_x1 <- c(cnv_x1,plotted_sample_num * sample_band_width + patient_margin * plotted_patient_num)
    cnv_x2 <- c(cnv_x2,plotted_sample_num * sample_band_width + patient_margin * plotted_patient_num + length(patient_sample_id_array) * sample_band_width)
    
    for (s_id in patient_sample_id_array) {
        
        patient_sample_frame <- subset(patient_frame,sample_id == s_id)
        patient_sample_frame$start <- patient_sample_frame$start+cumChrLen[as.character(patient_sample_frame$chr)]
        patient_sample_frame$end <- patient_sample_frame$end+cumChrLen[as.character(patient_sample_frame$chr)]
        patient_sample_frame$col <- cols[as.character(patient_sample_frame$cn)]
        if(FALSE){
            for( j in nrow(patient_sample_frame):1){
                
                cnv_y1 <- c(cnv_y1,patient_sample_frame$start[j])
                cnv_y2 <- c(cnv_y2,patient_sample_frame$end[j])
                cnv_x1 <- c(cnv_x1,plotted_sample_num + patient_margin * plotted_patient_num + sample_margin)
                cnv_x2 <- c(cnv_x2,plotted_sample_num + patient_margin * plotted_patient_num + plot_cnv_band_width + sample_margin)
                cnv_col <- c(cnv_col,patient_sample_frame$col[j])
                cnv_plot_frame <- data.frame(x1 = cnv_x1,x2 = cnv_x2,y1 = cnv_y1,y2 = cnv_y2)
                
            }
        }
        #cat(plotted_sample_num + patient_margin * plotted_patient_num + sample_margin,"\n")
        #cat(plotted_sample_num + patient_margin * plotted_patient_num + plot_cnv_band_width + sample_margin,"\n")
        patient_sample_frame$x1 <- plotted_sample_num * sample_band_width + patient_margin * plotted_patient_num + sample_margin
        patient_sample_frame$x2 <- plotted_sample_num * sample_band_width + patient_margin * plotted_patient_num + sample_band_width + sample_margin
        
        patient_sample_frame <- patient_sample_frame[order(patient_sample_frame$chr,patient_sample_frame$start,decreasing = TRUE),]
        all_patient_frame <- rbind(all_patient_frame,patient_sample_frame)
        plotted_sample_num = plotted_sample_num + 1
    }
    
    plotted_patient_num = plotted_patient_num + 1
}
band_plot_frame <- data.frame(x1 = cnv_x1,x2 = cnv_x2,y1 = cnv_y1,y2 = cnv_y2)
cnv_gp<-ggplot() + 
    geom_rect(band_plot_frame,mapping = aes(xmin=cnv_x1, xmax=cnv_x2, ymin=cnv_y1, ymax=cnv_y2),fill = "#F4F0EF") + 
    geom_rect(all_patient_frame,mapping = aes(xmin=x1, xmax=x2, ymin=start, ymax=end),fill = all_patient_frame$col) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "black",fill = NA),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.spacing = unit(0, 'mm'),
          plot.margin=unit(c(1,12 , 1, 1.9), "mm")) +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))




#sample cluster

reorder_sample_array <- sample_order_frame$sample_id

patient_array <- unique(sample_order_frame$patient_id)

sample_matrix <- matrix(0,length(patient_array),length(reorder_sample_array))

dimnames(sample_matrix) <- list(patient_array,reorder_sample_array)

for(patient in patient_array){
    for(sample in colnames(sample_matrix)){
        if(str_detect(sample,patient)){
            sample_matrix[patient,sample] = 1
        }
    }
}


update_sample_id_frame <- read.table(file = "E:\\project\\escc_multiregion\\sample_info\\updated_sample_id_index_list.xls",sep = "\t",header = T,stringsAsFactors = F)
update_patient_id_array <- unique(update_sample_id_frame$updated_patient_id)
names(update_patient_id_array) <- unique(update_sample_id_frame$patient_id)
row.names(sample_matrix) <- update_patient_id_array[row.names(sample_matrix)]


#sample_matrix <- sample_matrix[,reorder_sample_array]

sample_ht <- Heatmap(sample_matrix,
                     show_column_dend = FALSE,
                     border = TRUE,
                     column_names_gp = gpar(fontsize = c(6)),
                     rect_gp = gpar(lwd = 1,col = "white"),
                     col = structure(c("grey","grey"), names = c("1", "0")),
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     show_column_names = FALSE,
                     column_names_side = "bottom",
                     show_row_names = TRUE,
                     column_order = reorder_sample_array,
                     show_row_dend = FALSE,
                     row_names_gp = gpar(fontsize = c(6)),
                     #column_split = factor(sample_order_frame$patient_id,levels = unique(sample_order_frame$patient_id)),
                     #column_gap = unit(c(0.1), "mm"),
                     show_heatmap_legend = FALSE,
                     column_title = NULL,
                     height = unit(10,"cm"),
                     cell_fun = function(j,i,x,y,w,h,col){
                         r = min(unit.c(w, h)) * 1.5
                         #grid.rect(x = x,y = y,width = w, height = h, gp = gpar(fill = "grey",col = "white",alpha = 1))
                         if(sample_matrix[i,j] == 1){
                            grid.circle(x, y,r , gp = gpar(fill = "#339933", col = "grey"))
                         }
                         if(j == 1){
                             start_index <- which(sample_matrix[i,reorder_sample_array]==1)[1]
                             end_index <- which(sample_matrix[i,reorder_sample_array]==1)[length(which(sample_matrix[i,reorder_sample_array]==1))]
                             #cat(start_index,end_index,j,"\n")
                             start_x <- (start_index - 1) * w + 0.5 * w
                             end_x <- (end_index - 1) * w + 0.5 * w
                             grid.segments(start_x,y,end_x,y,gp=gpar(lwd = 2,col = "#339933"))
                         }
                     })



rna_heatmap = grid.grabExpr(draw(rna_ht,padding = unit(c(1, 2.9, 1, 12), "mm")))
rrbs_heatmap = grid.grabExpr(draw(rrbs_ht,padding = unit(c(1, 2.9, 1, 12), "mm")))

if(F){
cnv_gp = grid.grabExpr({
    cnv_gp
    print(cnv_gp)
})
}
sample_heatmap = grid.grabExpr(draw(sample_ht,padding = unit(c(1, 2.9, 1, 1.9), "mm")))

grid.newpage()

pushViewport(viewport(layout = grid.layout(nr = 10, nc = 1)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(rna_heatmap)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(rrbs_heatmap)
popViewport()

pushViewport(viewport(layout.pos.row = c(3), layout.pos.col = 1))
grid.draw(cnv_gp)
popViewport()

pushViewport(viewport(layout.pos.row = c(6,7), layout.pos.col = 1))
grid.draw(sample_heatmap)
popViewport()
















if(F){
dev.new()


rna_col_fun = colorRamp2(c(-4,0, 5), c("#4A1F67", "#3A6589", "#ECE129"))
rna_lgd = Legend(col_fun = rna_col_fun, title = "RNA expression",at = c(-4,0, 5),labels = c("low","median", "high"))

rrbs_col_fun = colorRamp2(c(0, 0.5, 1), c("#3A54A4", "white", "#EE2123"))
rrbs_lgd = Legend(col_fun = rrbs_col_fun, title = "Methylation level",at = c(0, 0.5, 1),labels = c("low","median", "high"))

cna_col_fun = colorRamp2(c(0, 0.5, 1,2), c("#EE2223",
                                         "#FFB5C5",
                                         "#F4F0EF",
                                         "#1874CD"))
cna_lgd = Legend(col_fun = cna_col_fun, title = "Copy number",at = c(0, 0.5, 1, 2),labels = c("amplification","gain", "neutral","loss"))



pd = packLegend(rna_lgd,rrbs_lgd,cna_lgd,max_width = unit(10, "cm"),column_gap = unit(1, "cm"))

pushViewport(viewport(width = 0.8, height = 0.8))
grid.rect()
draw(pd)
popViewport()
}








