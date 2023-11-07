#' Function to calculate the power spectral density (PSD)
#'
#' @description function to calculate the power spectral density (PSD) of a time series. Methods derived from A. Eke (2000) "Physiological time series: distinguishing fractal noises from motions" <doi:10.1007/s004249900135>.
#'
#' @param series a numeric vector, with data for a regularly spaced time series.
#' @param plot logical. whether to draw the plot of log power vs. log frequency
#' @param min,max the optional values. Frequency range of power spectral density. The default value is 1/2 and 1/8 and cannot be set to a negative number
#'
#' @return a value of spectral exponent(beta) which the the slope of the  of the fitting line on plot of log power vs. log frequency
#' @references Zhang T, Dong X, Chen C, Wang D, Zhang XD. RespirAnalyzer: an R package for continuous monitoring of respiratory signals.
#'
#' @examples data("TestData")
#' Fs <- 50
#' Peaks <- find.peaks(Data[,2],Fs,lowpass=TRUE,freq=1,MovingAv=FALSE,
#'                     W=FALSE,filter=TRUE,threshold=0.05)
#' PP_interval=diff(Peaks[,1])/Fs
#' LowPSD(series=PP_interval,plot=TRUE,min=1/64, max=1/2)
#'
LowPSD <- function(series,plot=TRUE,min=1/8, max=1/2){
  if (min < 0 || max < 0 || max < min)  stop(
    "'max' and 'min' must be positive value and max>min")
  A<- series-mean(series)
  n=length(series)
  x=1:n
  W<- 1-(2*x/(n+1)-1)^2
  B=A*W
  k=(B[n]-B[1])/(n-1)
  b=-(B[n]-B[1]*n)/(n-1)
  y=k*x+b
  C=B-y
  #####
  #####
  Y=fft(C)
  f=(2:ceiling(n/2))/n
  position=which(f>min&f<max)
  P=abs(Y*Conj(Y)/n)
  x=log10(f[position]);y=log10(P[position])
  fit<-lm(y~x)
  if(plot){
    plot(x,y,"l",
         xlab="log10(frequency)",ylab="log10(power)", axes=FALSE,cex.lab=1.6)
    axis(1); axis(2, las=2); box()
    cor(x,y)
    abline(fit)
  }
  slope <- coef(fit)
  return(slope)
}
