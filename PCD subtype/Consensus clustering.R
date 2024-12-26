#This code is suitable for consensus clustering and contour coefficient evaluation of samples#

silhouette_SimilarityMatrix<-function(group, similarity_matrix)
{
  similarity_matrix=as.matrix(similarity_matrix)
  similarity_matrix<-(similarity_matrix+t(similarity_matrix))/2
  diag(similarity_matrix)=0
  normalize <- function(X) X / rowSums(X)
  similarity_matrix<-normalize(similarity_matrix)
  
  n <- length(group)
  if(!all(group == round(group))) stop("'group' must only have integer codes")
  cluster_id <- sort(unique(group <- as.integer(group)))
  k <- length(cluster_id)
  if(k <= 1 || k >= n)
    return(NA)
  doRecode <- (any(cluster_id < 1) || any(cluster_id > k))
  if(doRecode)
    group <- as.integer(fgroup <- factor(group))
  cluster_id <- sort(unique(group))
  
  wds <- matrix(NA, n,3, dimnames =list(names(group), c("cluster","neighbor","sil_width")))  
  for(j in 1:k)
  { 
    index <- (group == cluster_id[j])
    Nj <- sum(index)
    wds[index, "cluster"] <- cluster_id[j]
    dindex <- rbind(apply(similarity_matrix[!index, index, drop = FALSE], 2,
                          function(r) tapply(r, group[!index], mean)))
    maxC <- apply(dindex, 2, which.max)
    wds[index,"neighbor"] <- cluster_id[-j][maxC]
    s.i <- if(Nj > 1) {
      a.i <- colSums(similarity_matrix[index, index])/(Nj - 1)
      b.i <- dindex[cbind(maxC, seq(along = maxC))]
      ifelse(a.i != b.i, (a.i - b.i) / pmax(b.i, a.i), 0)
    } else 0
    wds[index,"sil_width"] <- s.i
  }
  attr(wds, "Ordered") <- FALSE
  class(wds) <- "silhouette"
  wds
}

setwd("D:/Desktop/research19/02 cluster")

#Consensus clustering#
library(ConsensusClusterPlus)

data <- TCGA_TPM_log2x1_21415_368[rownames(TCGA_TPM_log2x1_21415_368) %in% union_189$DEG.uni.cox,]
data <- as.matrix(data)
cluster <- ConsensusClusterPlus(data,maxK=6,reps=3000,pItem=0.8,pFeature=1,
                                title="km-euclidean-6",clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot="png",writeTable = TRUE)

pdf("silhouette.pdf")  
sil <- silhouette_SimilarityMatrix(cluster[[3]]$consensusClass, cluster[[3]]$consensusMatrix)

# plot(sil)
colors <- c("#B42D34", "#1F4781","#00A087") 

plot(sil, col = colors)
dev.off()