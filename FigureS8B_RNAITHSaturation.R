library(ggplot2)

setwd("E:\\project\\escc_multiregion\\rna_analysis\\rna_ith_saturation")

patient_id_update_frame <- read.table(file = "E:\\project\\escc_multiregion\\sample_info\\updated_sample_id_index_list.xls", sep = "\t", header = T,stringsAsFactors = F)
updated_patient_id_array <- unique(patient_id_update_frame$updated_patient_id)
names(updated_patient_id_array) <- unique(patient_id_update_frame$patient_id)



ith_frame = read.table(file = "escc.all_patient_wise_ith_values_saturation_matrix.xls",sep = "\t",header = TRUE)
ith_frame$updated_patient_id <- updated_patient_id_array[ith_frame$patient_id]

ggplot(ith_frame,aes(x = sample_count,y = ith_score)) + 
    geom_point(col = "#41b6e6") + 
    geom_line(aes(x = sample_count,y = mean_value),size = 1,col = "#ff585d") +
    geom_line(aes(x = sample_count,y = mean_value + std_value),size = 1,col = "#ffb549") +
    geom_line(aes(x = sample_count,y = mean_value - std_value),size = 1,col = "#ffb549") +
    theme(panel.border = element_rect(fill = NA),strip.background = element_blank()) + 
    facet_wrap(updated_patient_id ~ .,nrow = 5) + 
    xlab("Number of Samples") + 
    ylab("Gene Expression ITH")



































