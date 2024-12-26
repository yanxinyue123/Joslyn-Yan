#This code is suitable for drawing survival curve graphs after consensus clustering

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

setwd("D:/Desktop/research19/02 cluster/")
datas <- read.csv("D:/Desktop/research19/group.csv",header = T,row.names = 1)
datas$time <- datas$time/30
fitd <- survdiff(Surv(time, event) ~ group,data = datas)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(time, event)~ group,data = datas,type = "kaplan-meier",error = "greenwood",conf.type = "plain")
# Paired Survival Analysis
ps <- pairwise_survdiff(Surv(time, event)~ group,data = datas,p.adjust.method = "none") # Correction is not used here. If correction is needed, none can be replaced with BH

#mycol <- brewer.pal(n = 10, "Paired")[c(2,4,6,8)]
mycol <- c("#00A087","#4DBBD5","#E64832","#3C5488")

names(fit$strata) <- gsub("group=", "", names(fit$strata))

p <- ggsurvplot(fit = fit,
                conf.int = FALSE, 
                risk.table = TRUE, 
                risk.table.col = "strata",
                palette = mycol, 
                data = datas,
                xlim = c(0,120), 
                size = 1,
                break.time.by = 12, 
                legend.title = "",
                xlab = "Time (months)",
                ylab = "Overall survival",
                risk.table.y.text = FALSE,
                tables.height = 0.3) 
## overall pvalue
## P<0.001——“<0.001”
#p.lab <- paste0("log-rank test P",ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))

#Scientific notation
options(scipen=999)  
options(digits=2) 
p.lab <- paste0("log-rank test P",
                paste0(" = ",sprintf("%.2e", p.val))
)

p$plot <- p$plot + annotate("text",
                            x = 0, y = 0.3, 
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
## Add pairing table
addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                         round(ps$p.value, 3))))
addTab[is.na(addTab)] <- "-"
df <- tibble(x = 0, y = 0, tb = list(addTab))
p$plot <- p$plot +
  geom_table(data = df,
             aes(x = x, y = y, label = tb),
             table.rownames = TRUE)
##Generate images
pdf("TCGA-survival.pdf", width = 4.5, height = 6)
print(p)
dev.off()

