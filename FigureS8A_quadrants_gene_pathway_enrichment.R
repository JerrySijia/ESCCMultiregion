library(ggplot2)
library(ReactomePA)
library(clusterProfiler)


setwd("E:\\project\\escc_multiregion\\tri-omics\\gene_quadrants_pathway_enrichment_analysis")


quadrants_gene_frame <- read.table(file = "escc.rna_heterogeneity_quadrants_matrix.xls",sep = "\t",header = TRUE)
q1_gene_cn_quadrants_frame <- subset(quadrants_gene_frame,rna_quadrants == "Q1")
q2_gene_cn_quadrants_frame <- subset(quadrants_gene_frame,rna_quadrants == "Q2")
q3_gene_cn_quadrants_frame <- subset(quadrants_gene_frame,rna_quadrants == "Q3")
q4_gene_cn_quadrants_frame <- subset(quadrants_gene_frame,rna_quadrants == "Q4")

q1_kegg <- bitr(q1_gene_cn_quadrants_frame$geneName, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db");

q1_enrich_result <- enrichPathway(q1_kegg$ENTREZID,pvalueCutoff = 0.05)
dotplot(q1_enrich_result,showCategory=10)


q2_kegg <- bitr(q2_gene_cn_quadrants_frame$geneName, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db");

q2_enrich_result <- enrichPathway(q2_kegg$ENTREZID,pvalueCutoff = 0.05)
barplot(q2_enrich_result,showCategory=10)

q3_kegg <- bitr(q3_gene_cn_quadrants_frame$geneName, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db");

q3_enrich_result <- enrichPathway(q3_kegg$ENTREZID,pvalueCutoff = 0.05)
barplot(q3_enrich_result,showCategory=10)

q4_kegg <- bitr(q4_gene_cn_quadrants_frame$geneName, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db");

q4_enrich_result <- enrichPathway(q4_kegg$ENTREZID,pvalueCutoff = 0.05)
barplot(q4_enrich_result,showCategory=10)


