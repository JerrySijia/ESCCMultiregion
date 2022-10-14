library(stringr)

setwd("E:\\project\\escc_multiregion\\dna_analysis\\msai_analysis")

chrLen <- c(249250621,
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
            51304566)
  
names(chrLen)<-c(1:22)

#abline(v = c(0, as.list(cumChrLen)), lty = 3) #

msai_snp_frame <- read.table(file = "ESCC34.msai_snp.tsv",sep = "\t",header = T)
sample_id_array <- c()

for(col in colnames(msai_snp_frame)){
  if(str_detect(col,"reads1")){
    sample_id_array <- c(sample_id_array,substr(col,0,7))
  }
}

msai_snp_frame_by_pos <- split(msai_snp_frame,list(msai_snp_frame$msai_chrom,msai_snp_frame$msai_start))

par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(length(sample_id_array),1), xaxt='n')

reference_sample_id <- "ESCC34C"
msai_snp_frame$aaf_reads1 <- "yes"
msai_snp_frame$aaf_reads1[msai_snp_frame[[paste0(reference_sample_id,"_tumor_reads1")]] < msai_snp_frame[[paste0(reference_sample_id,"_tumor_reads2")]]] <- "no"
chrom <- "6"
msai_snp_frame <- msai_snp_frame[msai_snp_frame$msai_chrom == chrom,]

x_start <- min(msai_snp_frame$msai_start)
x_end <- max(msai_snp_frame$msai_end)


for(sample_id in sample_id_array){
  sample_msai_snp_frame <- subset(msai_snp_frame,select = c("chrom","pos",paste0(sample_id,"_tumor_reads1"),paste0(sample_id,"_tumor_reads2"),"aaf_reads1"))
  sample_msai_snp_frame$aaf <- sample_msai_snp_frame[[paste0(sample_id,"_tumor_reads1")]] / (sample_msai_snp_frame[[paste0(sample_id,"_tumor_reads1")]] + sample_msai_snp_frame[[paste0(sample_id,"_tumor_reads2")]])
  sample_msai_snp_frame$aaf[sample_msai_snp_frame$aaf_reads1 == "no"] <- 1 - sample_msai_snp_frame$aaf[sample_msai_snp_frame$aaf_reads1 == "no"]
  sample_msai_snp_frame$baf <- 1 - sample_msai_snp_frame$aaf
  plot(ylab = "", 
       type = "n", 
       xaxt = "n", 
       #yaxt = "n",
       x = c(x_start,x_end), 
       y = c(-0.02, 1.1), 
       las = 1)
  points(sample_msai_snp_frame$pos,y = sample_msai_snp_frame$aaf,cex = 0.5,col = "#E63729")
  points(sample_msai_snp_frame$pos,y = sample_msai_snp_frame$baf,cex = 0.5,col = "#6891AA")
}
