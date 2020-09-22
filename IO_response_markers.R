### IO_response_markers.R ###
# William Lee @ September 2020

bdir = "~/Box/Huang_lab/manuscripts/BCellProject"
setwd(bdir)

### extraction of canonical cell markers ###

# read in immune cell frequency features (xCell)
xCell_freq = "~/Box/Huang_lab/Huang_lab_data/xCell_Aran_GenomeBio2017/xCell_TCGA_RSEM.txt"
xCell_freq_df = read.table(sep="\t", header=T, file=xCell_freq, stringsAsFactors=FALSE)

# read in expression file
cell_mark_expr = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv.gz"
cell_mark_expr_df = read.table(sep="\t", header=T, file=cell_mark_expr, stringsAsFactors=FALSE)

# read in subtype file
subtype_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/Subtype\ Assignments.xlsx"
subtype = data.frame(readxl::read_xlsx(subtype_f))
colnames(subtype)[1] = "bcr_patient_barcode"

# get BCR patient barcodes, transpose xCell_freq_df, merge with subtype file
colnames(xCell_freq_df) <- gsub('\\.', '-', colnames(xCell_freq_df))
colnames(xCell_freq_df) <- substr(colnames(xCell_freq_df), 1, 12)
xCell_freq_df_t <- as.data.frame(t(xCell_freq_df[-1]))
colnames(xCell_freq_df_t) <- xCell_freq_df[, 1]
xCell_freq_df_t$bcr_patient_barcode <- rownames(xCell_freq_df_t); rownames(xCell_freq_df_t) <- NULL
xCell_freq_merge = merge(subtype, xCell_freq_df_t, by="bcr_patient_barcode") # 8225 x 68

# identify each COAD as COAD_MSI or COAD_MSS
xCell_freq_merge$DISEASE[xCell_freq_merge$DISEASE=="COAD" & xCell_freq_merge$SUBTYPE=="MSI"] = "COAD_MSI"
xCell_freq_merge$DISEASE[xCell_freq_merge$DISEASE=="COAD" & xCell_freq_merge$SUBTYPE!="MSI"] = "COAD_MSS"

# only keep cancers that have ORR + TMB information and don't arise from lymphatic tissue
xCell_freq_merge = xCell_freq_merge[xCell_freq_merge$DISEASE!="CHOL" & xCell_freq_merge$DISEASE!="DLBC" &
                                    xCell_freq_merge$DISEASE!="KIRC" & xCell_freq_merge$DISEASE!="KIRP" &
                                    xCell_freq_merge$DISEASE!="LAML" & xCell_freq_merge$DISEASE!="LGG"  &
                                    xCell_freq_merge$DISEASE!="PCPG" & xCell_freq_merge$DISEASE!="READ" &
                                    xCell_freq_merge$DISEASE!="STAD" & xCell_freq_merge$DISEASE!="THCA" & 
                                    xCell_freq_merge$DISEASE!="TGCT" & xCell_freq_merge$DISEASE!="THYM" &
                                    xCell_freq_merge$DISEASE!="UCEC",]

##########################################################################################################

# get average CD19 expression by cancer #
cd19=cell_mark_expr_df[cell_mark_expr_df$gene_id=="CD19|930",]

# get BCR patient barcodes, transpose cd19, merge with xCell_freq_merge
colnames(cd19) <- gsub('\\.', '-', colnames(cd19))
colnames(cd19) <- substr(colnames(cd19), 1, 12)
cd19_t <- as.data.frame(t(cd19[-1]))
colnames(cd19_t) <- cd19[, 1]
cd19_t$bcr_patient_barcode <- rownames(cd19_t); rownames(cd19_t) <- NULL
cd19_merge = merge(xCell_freq_merge, cd19_t, by="bcr_patient_barcode") # 5370 x 69

# generate CD19 summary df (cd19_df) 
cd19_df <- data.frame(matrix(ncol = 3, nrow = 0))
for (cancer in unique(cd19_merge$DISEASE)) {
  
  cd19_temp = cd19_merge[cd19_merge$DISEASE==cancer,]
  samp_size = nrow(cd19_temp)
  mean_expr = mean(cd19_temp$`CD19|930`)
  
  temp = cbind(cancer, samp_size, mean_expr)
  
  cd19_df = rbind(cd19_df, temp)
  
}

colnames(cd19_df) <- c("cancer", "sample_size", "avg_CD19_expr")

tt=cd19_df[order(cd19_df$cancer, decreasing=FALSE),]
tn = "out/CD19_expression_by_cancer.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

##########################################################################################################

# get average CD27 expression by cancer #
cd27=cell_mark_expr_df[cell_mark_expr_df$gene_id=="CD27|939",]

# get BCR patient barcodes, transpose cd27, merge with xCell_freq_merge
colnames(cd27) <- gsub('\\.', '-', colnames(cd27))
colnames(cd27) <- substr(colnames(cd27), 1, 12)
cd27_t <- as.data.frame(t(cd27[-1]))
colnames(cd27_t) <- cd27[, 1]
cd27_t$bcr_patient_barcode <- rownames(cd27_t); rownames(cd27_t) <- NULL
cd27_merge = merge(xCell_freq_merge, cd27_t, by="bcr_patient_barcode") # 5370 x 69

# generate CD27 summary df (cd27_df) 
cd27_df <- data.frame(matrix(ncol = 3, nrow = 0))
for (cancer in unique(cd27_merge$DISEASE)) {
  
  cd27_temp = cd27_merge[cd27_merge$DISEASE==cancer,]
  samp_size = nrow(cd27_temp)
  mean_expr = mean(cd27_temp$`CD27|939`)
  
  temp = cbind(cancer, samp_size, mean_expr)
  
  cd27_df = rbind(cd27_df, temp)
  
}

colnames(cd27_df) <- c("cancer", "sample_size", "avg_CD27_expr")


tt=cd27_df[order(cd27_df$cancer, decreasing=FALSE),]
tn = "out/CD27_expression_by_cancer.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

##########################################################################################################

# get average CD40 expression by cancer #
cd40=cell_mark_expr_df[cell_mark_expr_df$gene_id=="CD40|958",]

# get BCR patient barcodes, transpose cd40, merge with xCell_freq_merge
colnames(cd40) <- gsub('\\.', '-', colnames(cd40))
colnames(cd40) <- substr(colnames(cd40), 1, 12)
cd40_t <- as.data.frame(t(cd40[-1]))
colnames(cd40_t) <- cd40[, 1]
cd40_t$bcr_patient_barcode <- rownames(cd40_t); rownames(cd40_t) <- NULL
cd40_merge = merge(xCell_freq_merge, cd40_t, by="bcr_patient_barcode") # 5370 x 69

# generate CD40 summary df (cd40_df) 
cd40_df <- data.frame(matrix(ncol = 3, nrow = 0))
for (cancer in unique(cd40_merge$DISEASE)) {
  
  cd40_temp = cd40_merge[cd40_merge$DISEASE==cancer,]
  samp_size = nrow(cd40_temp)
  mean_expr = mean(cd40_temp$`CD40|958`)
  
  temp = cbind(cancer, samp_size, mean_expr)
  
  cd40_df = rbind(cd40_df, temp)
  
}

colnames(cd40_df) <- c("cancer", "sample_size", "avg_CD40_expr")

tt=cd40_df[order(cd40_df$cancer, decreasing=FALSE),]
tn = "out/CD40_expression_by_cancer.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

##########################################################################################################

# get average CD8A expression by cancer #
cd8A=cell_mark_expr_df[cell_mark_expr_df$gene_id=="CD8A|925",]

# get BCR patient barcodes, transpose cd8A, merge with xCell_freq_merge
colnames(cd8A) <- gsub('\\.', '-', colnames(cd8A))
colnames(cd8A) <- substr(colnames(cd8A), 1, 12)
cd8A_t <- as.data.frame(t(cd8A[-1]))
colnames(cd8A_t) <- cd8A[, 1]
cd8A_t$bcr_patient_barcode <- rownames(cd8A_t); rownames(cd8A_t) <- NULL
cd8A_merge = merge(xCell_freq_merge, cd8A_t, by="bcr_patient_barcode") # 5370 x 69

# generate CD8A summary df (cd8A_df) 
cd8A_df <- data.frame(matrix(ncol = 3, nrow = 0))
for (cancer in unique(cd8A_merge$DISEASE)) {
  
  cd8A_temp = cd8A_merge[cd8A_merge$DISEASE==cancer,]
  samp_size = nrow(cd8A_temp)
  mean_expr = mean(cd8A_temp$`CD8A|925`)
  
  temp = cbind(cancer, samp_size, mean_expr)
  
  cd8A_df = rbind(cd8A_df, temp)
  
}

colnames(cd8A_df) <- c("cancer", "sample_size", "avg_CD8A_expr")


tt=cd8A_df[order(cd8A_df$cancer, decreasing=FALSE),]
tn = "out/CD8A_expression_by_cancer.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

##########################################################################################################

# get average CD45RO expression by cancer #
cd45RO=cell_mark_expr_df[cell_mark_expr_df$gene_id=="PTPRC|5788",]

# get BCR patient barcodes, transpose cd45RO, merge with xCell_freq_merge
colnames(cd45RO) <- gsub('\\.', '-', colnames(cd45RO))
colnames(cd45RO) <- substr(colnames(cd45RO), 1, 12)
cd45RO_t <- as.data.frame(t(cd45RO[-1]))
colnames(cd45RO_t) <- cd45RO[, 1]
cd45RO_t$bcr_patient_barcode <- rownames(cd45RO_t); rownames(cd45RO_t) <- NULL
cd45RO_merge = merge(xCell_freq_merge, cd45RO_t, by="bcr_patient_barcode") # 5370 x 69

# generate CD45RO summary df (cd45RO_df) 
cd45RO_df <- data.frame(matrix(ncol = 3, nrow = 0))
for (cancer in unique(cd45RO_merge$DISEASE)) {
  
  cd45RO_temp = cd45RO_merge[cd45RO_merge$DISEASE==cancer,]
  samp_size = nrow(cd45RO_temp)
  mean_expr = mean(cd45RO_temp$`PTPRC|5788`)
  
  temp = cbind(cancer, samp_size, mean_expr)
  
  cd45RO_df = rbind(cd45RO_df, temp)
  
}

colnames(cd45RO_df) <- c("cancer", "sample_size", "avg_CD45RO_expr")

tt=cd45RO_df[order(cd45RO_df$cancer, decreasing=FALSE),]
tn = "out/CD45RO_expression_by_cancer.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

##########################################################################################################

### cross correlation of markers and immune cell frequency features (xCell) ###

library(PerformanceAnalytics)

cd27_select = select(cd27_merge, "bcr_patient_barcode", "CD27|939") # B- and T-cell marker

# B-cell markers #
cd40_select = select(cd40_merge, "bcr_patient_barcode", "CD40|958")
B_cell_merge = merge(cd19_merge, cd27_select, by="bcr_patient_barcode")
B_cell_merge = merge(B_cell_merge, cd40_select, by="bcr_patient_barcode")
B_cell_cc = select(B_cell_merge, "CD19|930", "CD27|939", "CD40|958", "B-cells", "Memory B-cells", "naive B-cells")
pdf("out/BCell_cross_correlation.pdf") # open pdf file
chart.Correlation(B_cell_cc, method="pearson", histogram=TRUE, pch=16) # generate plot
dev.off() # close pdf file

# T-cell markers #
cd45RO_select = select(cd45RO_merge, "bcr_patient_barcode", "PTPRC|5788")
T_cell_merge = merge(cd8A_merge, cd27_select, by="bcr_patient_barcode")
T_cell_merge = merge(T_cell_merge, cd45RO_select, by="bcr_patient_barcode")
T_cell_cc = select(T_cell_merge, "CD8A|925", "CD27|939", "PTPRC|5788", "CD8+ T-cells", "CD8+ Tem")
pdf("out/TCell_cross_correlation.pdf") # open pdf file
chart.Correlation(T_cell_cc, method="pearson", histogram=TRUE, pch=16) # generate plot
dev.off() # close pdf file