#This code is used to screen genes related to survival using univariate Cox regression

#Input gene expression and sample survival data, rows represent samples
td <- read.csv("D:/Desktop/research/ICGC_log2_input.csv",sep = ",")

library(survival)  
library(multtest)
pFilter=0.05 
outResult=data.frame() 
sigGenes=c("surstat","surtime") 
for(i in colnames(td[,3:ncol(td)])){ 
  tdcox <- coxph(Surv(surtime, surstat) ~ td[,i], data = td)
  tdcoxSummary = summary(tdcox) 
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] 
  
  if(pvalue<pFilter){ 
    sigGenes=c(sigGenes,i)
    p_adjusted <- p.adjust(pvalue, method = "BH") 
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"],
                          p_adjusted=p_adjusted
                    )
    )
  }
}

write.csv(outResult,file="D:/Desktop/research19/unicox_result.csv",row.names=F,quote=F)

UniCoxSurSigGeneExp=td[,sigGenes] 
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)
write.csv(UniCoxSurSigGeneExp,file="D:/Desktop/research19/01 data/TCGA/TCGA_unicox_keep.csv",row.names=F,quote=F)


