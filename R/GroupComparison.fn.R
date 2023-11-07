  #' Function to calculate the statistics for each Group
  #'
  #' @description  function to calculate the statistics for each Group: Number of Samples,
  #' mean, standard deviation (SD), standard error (SE), median, confident interval, p-value of ANOVA
  #'
  #' @param  Data vector for response values
  #' @param  GroupName vector for group names
  #' @param  na.rm whether to remove value for calculation
  #' @param  conf.level confidence level
  #'
  #' @return  a dataframe for the statistics for each Group. Number of Samples, mean
  #' @return  standard deviation (SD), median,  upper and lower bounds of CI, p-value of ANOVA
  #'
  #' @references Zhang T, Dong X, Chen C, Wang D, Zhang XD. RespirAnalyzer: an R package for continuous monitoring of respiratory signals.
  #' @examples data("HqData")
  #' PP_Hq <- HqData
  #' filenames <- row.names(PP_Hq)
  #' q <- -10:10
  #' ClassNames <- c(substr(filenames[1:19], start = 1, stop = 3),
  #'                 substr(filenames[20:38], start = 1, stop = 5))
  #' Class <- unique(ClassNames)
  #' Data <- GroupComparison.fn(PP_Hq[,1],ClassNames)
  #' Data
  #'

GroupComparison.fn <- function(Data,GroupName, na.rm=TRUE, conf.level = 0.95)
{
  N <- tapply(Data,GroupName,function(x){sum(!is.na(x))})
  Mean <- tapply(Data,GroupName, mean, na.rm = na.rm)
  SD <- tapply(Data,GroupName, sd, na.rm = na.rm)
  SE <- SD/(N^1/2)
  median <- tapply(Data,GroupName, median, na.rm = na.rm)

  Class <- unique(GroupName)
  CIlower_vec <- CIupper_vec <-  pval_vec <- rep(NA, length(Class))
  idx=0
  for (i in 1:length(Class)){
    idx=idx+1
    data <- Data[which(GroupName == Class[i])]
    Ttest_list <- t.test(data,conf.level = conf.level)
    CIlower_vec[idx] <- Ttest_list$conf.int[1]
    CIupper_vec[idx] <- Ttest_list$conf.int[2]
  }
  ANOVA <- summary(a.aov <- aov(Data~GroupName, data= data.frame(Data,GroupName)))
  pval_vec <- rep(ANOVA[[1]]["Pr(>F)"][[1]][1], length(Class))
  Results <- cbind(N,Mean,SD,SE,median,CIlower_vec, CIupper_vec,pval_vec)

  dimnames(Results) <- list(Group <- Class,SummaryStat=c("N","Mean","SD","SE","Median",
                                                         "CIlower","CIupper","P-value"))

  return(Results)
}
