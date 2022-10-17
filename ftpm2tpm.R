countToTpm <- function(counts, effLen)
countToTpm <- function(counts, effLen)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}
 
countToFpkm <- function(counts, effLen)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
 
fpkmToTpm <- function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
 
countToEffCounts <- function(counts, len, effLen)
{
    counts * (len / effLen)
}

expMatrix<-read.table("C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/B16_mouse_FPKM.txt",header = T,row.names = 1)
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[10000:10003,]
colSums(tpms) #根据TPM的特征进行检查，看每列加和是否一致
write.table(tpms,"C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/B16_mouse_TPM.txt",quote = FALSE, sep = "\t",row.names = TRUE,col.names = TRUE)