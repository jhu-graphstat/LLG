rm(list = ls())
# setwd("E:/GitHub/LLG/Code/R")
setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

Diff_ZG_Mean = list()
Diff_ZG_UB = list()
Diff_ZG_LB = list()
Diff_USVT_Mean = list()
Diff_USVT_UB = list()
Diff_USVT_LB = list()

dataName = "JHU"
# dataName = "desikan"
# dataName = "CPAC200"

mVec = c(1, 5, 10, 50)

for (iM in 1:length(mVec)) {
  m = mVec[iM]
  
  # Eigen-decomposition
  fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_eig.RData", sep="")
  load(fileName)
  
  error_A_bar = error_A_bar^2*n*(n-1)
  error_P_hat = error_P_hat^2*n*(n-1)
  
  error_P_hat_ZG = error_P_hat_ZG^2*n*(n-1)
  error_P_hat_USVT = error_P_hat_USVT^2*n*(n-1)
  
  Diff_ZG_Mean[[dataName]][iM] = mean(error_A_bar - error_P_hat_ZG)
  Diff_ZG_LB[[dataName]][iM] = Diff_ZG_Mean[[dataName]][iM] -
    sqrt(var(error_A_bar - error_P_hat_ZG))/sqrt(length(error_A_bar))*1.96
  Diff_ZG_UB[[dataName]][iM] = Diff_ZG_Mean[[dataName]][iM] +
    sqrt(var(error_A_bar - error_P_hat_ZG))/sqrt(length(error_A_bar))*1.96
  
  Diff_USVT_Mean[[dataName]][iM] = mean(error_A_bar - error_P_hat_USVT)
  Diff_USVT_LB[[dataName]][iM] = Diff_USVT_Mean[[dataName]][iM] -
    sqrt(var(error_A_bar - error_P_hat_USVT))/sqrt(length(error_A_bar))*1.96
  Diff_USVT_UB[[dataName]][iM] = Diff_USVT_Mean[[dataName]][iM] +
    sqrt(var(error_A_bar - error_P_hat_USVT))/sqrt(length(error_A_bar))*1.96
  
  print(sqrt(var(error_P_hat_ZG/error_A_bar))/sqrt(length(error_A_bar))*1.96)
  print(sqrt(var(error_P_hat_USVT/error_A_bar))/sqrt(length(error_A_bar))*1.96)
  
  print(sum(error_A_bar > error_P_hat_ZG)/length(error_A_bar))
  print(sum(error_A_bar > error_P_hat_USVT)/length(error_A_bar))
}

