library(ggplot2)
library(vcd)
library(reshape2)


setwd("E:\\project\\escc_multiregion\\tri-omics\\gene_quadrants_in_cn_status_radar")

gene_cn_quadrants_frame <- read.table(file = "escc.gene_cn_quadrants_matrix.xls",sep = "\t",header = TRUE)

gene_cn_quadrants_table <- table(gene_cn_quadrants_frame$cn_status,gene_cn_quadrants_frame$quadrants)

cn_status_sum_matrix <- rowSums(gene_cn_quadrants_table)
quadrants_sum_matrix <- colSums(gene_cn_quadrants_table)

all_result_frame <- data.frame()
for (q in colnames(gene_cn_quadrants_table)) {
    for(cn_status in row.names(gene_cn_quadrants_table)){
        cn_quadrants_table <- matrix(c(gene_cn_quadrants_table[cn_status,q],
               cn_status_sum_matrix[cn_status] - gene_cn_quadrants_table[cn_status,q],
               quadrants_sum_matrix[q] - gene_cn_quadrants_table[cn_status,q],
               sum(gene_cn_quadrants_table) - sum(quadrants_sum_matrix[q]) - sum(cn_status_sum_matrix[cn_status]) + gene_cn_quadrants_table[cn_status,q]
               ),nrow = 2,dimnames = list(c(cn_status,paste("not_",cn_status,sep = "")),c(q,paste("not_",q,sep = ""))))
        print(cn_quadrants_table)
        stats <- fisher.test(cn_quadrants_table)
        result_frame <- data.frame("quadrants" = c(q),"scna_status" = c(cn_status),"oddratio" = c(stats$estimate),"pvalue" = c(stats$p.value))
        all_result_frame <- rbind(all_result_frame,result_frame)
    }
}
all_result_frame$logp <- -log10(all_result_frame$pvalue)

ggplot(all_result_frame,aes(x = oddratio,y = logp,col = quadrants,shape = scna_status)) + geom_point(size = 3) + 
    scale_x_continuous(breaks = c(0.9,0.95,1,1.05,1.1),labels = c("0.5","0.75","1","1.25","1.5")) +
    scale_y_continuous(breaks = c(25,50,75,100,125,150),labels = c("0.001","10-4","10-5","10-6","10-7","10-8")) +
    geom_vline(xintercept = 1,linetype = "dashed") +
    geom_hline(yintercept = 1.30103,linetype = "dashed") +
    xlab("Odds ratio") +
    ylab("Fisher¡¯s exact P") + 
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(size = 1)
          )

