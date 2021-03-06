
rm(list = ls())

setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

isSVD = 0

nIter = 100
nCores = 2

mVec = c(1,2,5,10)

dataName = "CPAC200"
# dataName = "desikan"
# dataName = "JHU"
# dataName = "slab907"
# dataName = "slab1068"
# dataName = "Talairach"

source("function_collection.R")
require(parallel)

tmpList = read_data(dataName, DA=F)
A_all = tmpList[[1]]
n = tmpList[[2]]
M = tmpList[[3]]
rm(tmpList)

P = add(A_all)/M
diag(P) = 0

dVec = 1:n
nD = length(dVec)

for (m in mVec) {
  print(c(m, isSVD))
  
  error_P_hat = matrix(0, nD, nIter)
  error_A_bar = matrix(0, nD, nIter)
  
  out <- mclapply(1:nIter, function(x) test_rank1(m, dVec, P, isSVD), 
                  mc.cores=nCores)
  out = array(unlist(out), dim = c(nD, nIter))
  
  error_P_hat = out
  
  if (isSVD) {
    fileName = paste("../../Result/result_", dataName, "_brute_rank1_",
                     "m_", m, "_svd.RData", sep="")
  } else {
    fileName = paste("../../Result/result_", dataName, "_brute_rank1_",
                     "m_", m, "_eig.RData", sep="")
  }
  
  save(error_P_hat, n, m, dVec, nIter, file=fileName)
  
}