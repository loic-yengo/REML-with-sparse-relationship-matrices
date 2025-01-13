library(MASS)
## Simulate data for example...
hsq <- 0.6
csq <- 0.1
Np <- 100000 # number of sibling pairs
G <- rnorm(n=Np,mean=0.5,sd=0.038) ## IBD sharing
Y <- do.call("rbind",lapply(1:Np,function(i){
  mvrnorm(n=1,mu=c(0,0),Sigma = matrix(c(1,hsq * G[i] + csq,hsq * G[i] + csq,1),2,2) )
}))


## write indiv. ID file
iid <- paste0(rep(paste0("FAM",1:Np),each=2),rep(c("_SIB1","_SIB2"),Np))
write.table(cbind(iid,iid),"Example_data.grm.id",quote=F,row.names = F,col.names = F,sep="\t")

## write relationship file
N <- 2 * Np
grm <- rbind(
  cbind(0:(N-1),0:(N-1),1),
  cbind(seq(0,N-1,by=2),seq(1,N-1,by=2),G)
)
write.table(grm,"Example_data.grm.sp",quote=F,row.names = F,col.names = F,sep="\t")

## write covariates file
covars <- cbind.data.frame(iid,iid,matrix(rnorm(N*5),nrow=N,ncol=5))
write.table(covars,"Example_data.covar.txt",quote=F,row.names = F,col.names = F,sep="\t")

## write phenotype file
write.table(cbind(iid,iid,c(t(Y))),"Example_data.pheno.txt",quote=F,row.names = F,col.names = F,sep="\t")


## Run test
source("sparseREML_v0.6_toshare.R")
model_AE   <- sparseREML(prefixGRM="Example_data",
                         phenoFile="Example_data.pheno.txt",
                         covarFile="Example_data.covar.txt",
                         model="ACE",
                         verbose = TRUE,algorithm = "ML")

print(model_AE)

