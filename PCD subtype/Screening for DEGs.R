
#This code is applicable for identifying differentially expressed genes between different groups

Find_DEG <- function(Data_C= NULL,
                     Data_T= NULL,
                     clinical_sample = NULL,
                     Threshold_LogFC = 1,
                     Threshold_P = 0.05,
                     Threshold_mean = 1){
  time_strat <- Sys.time()
  DEG_Result <- as.data.frame(matrix(0,nrow = length(rownames(Data_C)),ncol = 7))
  colnames(DEG_Result) <- c('Gene',
                            'Mean_C',
                            'Mean_T',
                            'LogFC',
                            'P.value',
                            'P.adjust',
                            'Status')
  DEG_Result$Gene <- rownames(Data_C)
  print('Start differential analysis')
  for (i in 1:length(DEG_Result$Gene)) {
    Control <- as.numeric(Data_C[i,])
    Test <- as.numeric(Data_T[i,])
    DEG_Result$Mean_C[i] <- mean(Control)
    DEG_Result$Mean_T[i] <- mean(Test)
    DEG_Result$LogFC[i] <- log2(DEG_Result$Mean_T[i]/DEG_Result$Mean_C[i])
    P <- wilcox.test(Test,Control)
    DEG_Result$P.value[i] <- P$p.value
  }
  #Calculate the adjusted P-value
  DEG_Result$P.adjust <- p.adjust(DEG_Result$P.value,method = 'BH')
  print('Label status')
  DEG_Result$Status <- 'None'
  a <- which(DEG_Result$LogFC > Threshold_LogFC & DEG_Result$P.adjust < Threshold_P)
  DEG_Result$Status[a] <- 'Up'
  a <- which(DEG_Result$LogFC < Threshold_LogFC*(-1) & DEG_Result$P.adjust < Threshold_P)
  DEG_Result$Status[a] <- 'Down'
  print('Differential analysis completed! Using time')
  print(Sys.time()-time_strat)
  # If the differential genes with mean values less than the threshold in both normal and tumor samples are assigned a value of None
  a <- which(DEG_Result$Mean_C < 1 & DEG_Result$Mean_T <1)
  DEG_Result$Status[a] <- 'None'
  cat('There are',length(DEG_Result$Status[DEG_Result$Status == 'Up']),'genes with increased expression','\n')
  cat('There are',length(DEG_Result$Status[DEG_Result$Status == 'Down']),'genes with reduced expression','\n')
  return(DEG_Result)
}

#Upload your own gene expression file and preprocess it
TCGA_TPM <- read.csv("D:/Desktop/research/TCGA_TPM.csv",row.names = 1)
TCGA_TPM_huanyuan <- 2^TCGA_TPM-0.001
TCGA_TPM_huanyuan[TCGA_TPM_huanyuan < 0] <- 0
row_zeros <- apply(TCGA_TPM_huanyuan, 1, function(x) sum(x == 0))
genes_to_keep <- (row_zeros / ncol(TCGA_TPM_huanyuan) <= 0.3)
#Remove genes that were not expressed in over 30% of the samples
exp_delete30 <- TCGA_TPM_huanyuan[genes_to_keep, ]
#Extract the expression level of cell death genes to be analyzed
common_1226 <- read.csv("D:/Desktop/research/01 data/common_all_1226.csv")
exp_delete30_1226 <- exp_delete30[rownames(exp_delete30) %in% common_1226$common_genes,]


#Load tumor and non tumor annotation files
subtype <- read.csv("D:/Desktop/research/01 data/TCGA-LIHC/tumor_normal529.csv")
Sample_A <- subtype$Samples[subtype$type == 'normal']
Sample_B <- subtype$Samples[subtype$type == 'tumor']
Sample_B <- Sample_B[-1]

#Sample grouping
merge_RNA_ComBat_A <- exp_delete30_1226[,colnames(exp_delete30_1226) %in% Sample_A]
merge_RNA_ComBat_B <- exp_delete30_1226[,colnames(exp_delete30_1226) %in% Sample_B]

#Run Function
merge_sutype_DEG <- Find_DEG(Data_C = merge_RNA_ComBat_A,
                             Data_T = merge_RNA_ComBat_B,
                             clinical_sample = clinical_sample,
                             Threshold_LogFC = 1,#Set threshold
                             Threshold_P = 0.05,
                             Threshold_mean = 1)

#Output the result file
setwd("D:/Desktop/research/05 groupDEGs")
write.csv(merge_sutype_DEG,"tumor_nomral_DEGs.csv")

