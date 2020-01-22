##############################################
# CovSel simulations
#
# January 2020
# Author: K Luijken
##############################################

# S1, seeds
set.seed(1)
seedsS1 <- matrix(sample(1:10000000, size=nrow(S1)*nsim), nrow=nrow(S1))
saveseedsS1 <- lapply(1:nrow(seedsS1), FUN= function(x) saveRDS(seedsS1[x,],file = paste0(filepath,"Siminput/S1/seedsS1_",x,".rds")))
set.seed(11)
bsseedsS1 <- array(replicate(nrow(S1)*nsim,sample(1:1000000000, size = B, replace=F)),dim=c(nrow=B,ncol=nsim,height=nrow(S1)))
savebsseedsS1 <- lapply(1:(dim(bsseedsS1)[3]), FUN= function(x) saveRDS(bsseedsS1[,,x],file = paste0(filepath,"Siminput/S1/bsseedsS1_",x,".rds")))


# S2, seeds
set.seed(2)
seedsS2 <- matrix(sample(1:10000000, size=nrow(S2)*nsim), nrow=nrow(S2))
saveseedsS2 <- lapply(1:nrow(seedsS2), FUN= function(x) saveRDS(seedsS2[x,],file = paste0(filepath,"Siminput/S2/seedsS2_",x,".rds")))
set.seed(22)
bsseedsS2 <- array(replicate(nrow(S2)*nsim,sample(1:1000000000, size = B, replace=F)),dim=c(nrow=B,ncol=nsim,height=nrow(S2)))
savebsseedsS2 <- lapply(1:(dim(bsseedsS2)[3]), FUN= function(x) saveRDS(bsseedsS2[,,x],file = paste0(filepath,"Siminput/S2/bsseedsS2_",x,".rds")))


# S3, seeds
set.seed(3)
seedsS3 <- matrix(sample(1:10000000, size=nrow(S3)*nsim), nrow=nrow(S3))
saveseedsS3 <- lapply(1:nrow(seedsS3), FUN= function(x) saveRDS(seedsS3[x,],file = paste0(filepath,"Siminput/S3/seedsS3_",x,".rds")))
set.seed(33)
bsseedsS3 <- array(replicate(nrow(S3)*nsim,sample(1:1000000000, size = B, replace=F)),dim=c(nrow=B,ncol=nsim,height=nrow(S3)))
savebsseedsS3 <- lapply(1:(dim(bsseedsS3)[3]), FUN= function(x) saveRDS(bsseedsS3[,,x],file = paste0(filepath,"Siminput/S3/bsseedsS3_",x,".rds")))
