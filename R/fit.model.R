#' Function to fit the MFDFA result with the extended binomial multifractal model.
#'
#' @description function to fit the result of Multifractal detrended fluctuation analysis (MFDFA) with the extended binomial multifractal model.
#' Return the results as a vector which contain the parameters of the model and the goodness of fit
#'
#' @param Hq a nurmeric vector for the generalized Hurst exponent.
#' @param q a vector of integers, q-order of the moment.
#'
#' @return  a vector for fitting parameters."a" and "b" is the coefficients of the extended binomial multifractal model.
#' "Goodness" is the goodness of fit
#' @references Zhang T, Dong X, Chen C, Wang D, Zhang XD. RespirAnalyzer: an R package for continuous monitoring of respiratory signals.
#'
#' @examples data("TestData") # load Data from TestData dataset
#' Fs=50 ## sampling frequency is 50Hz
#' Peaks=find.peaks(Data[,2],Fs)
#' PP_interval=diff(Peaks[,1])/Fs
#' exponents=seq(3, 8, by=1/4)
#' scale=2^exponents
#' q=-10:10
#' m=2
#' Result <- MFDFA(PP_interval, scale, m, q)
#' Coeff <- fit.model(Result$Hq,q)
#' Coeff



fit.model<-function(Hq,q){
  Hurst=data.frame(q,Hq);names(Hurst)<-c("q","Hq")
  if (0 %in% q)
    Hurst <- Hurst[-which(q==0),]

  a0=exp(-Hurst$Hq[1]/log(2));b0=exp(-Hurst$Hq[length(Hurst$q)]/log(2))
  m <- nls(Hurst$Hq ~ (1/Hurst$q)*(1-log(a^Hurst$q+b^Hurst$q)/log(2)),
           start=list(a=a0,b=b0),data=Hurst)

  Goodness <- cor(Hurst$Hq, predict(m))
  A <- summary(m)
  Parameter <- c(A$parameters[,1],Goodness)
  names(Parameter) <- c("a","b","Goodness")
  return(Parameter)
}
