# RespirAnalyzer
 Analysis Functions of Respiratory Data
# Description
 Provides functions for the complete analysis of respiratory data. Consists of a set of functions that allow to preprocessing respiratory data,  calculate both regular statistics and nonlinear statistics, conduct group comparison and visualize the results. Especially, Power Spectral Density ('PSD') (A. Eke (2000) <doi:10.1007/s004249900135>), 'MultiScale Entropy(MSE)' ('Madalena Costa(2002)' <doi:10.1103/PhysRevLett.89.068102>) and 'MultiFractal Detrended Fluctuation Analysis(MFDFA)' ('Jan W.Kantelhardt' (2002) <doi:10.1016/S0378-4371(02)01383-3>) were applied for the analysis of respiratory data. 
# Install
```R
#install "remotes" if it was not installed before
#install.packages("remotes")
remotes::install_github("dongxinzheng/RespirAnalyzer")
```
# Example
```R
# load Data from TestData dataset. The TestData is the respiratory data obtained from a healthy people which can download from the Fantasia database (PhysioNet, https://physionet.org/content/fantasia/1.0.0/)
data("TestData")
Seriesplot.fn(Data[1:2000,1],Data[1:2000,2],points=FALSE,
              xlab="Time(s)",ylab="Respiratory")
Fs=50 ## sampling frequency of the TestData is 50Hz
#### find the peaks of the TestData
Peaks <- find.peaks(Data[,2],Fs,lowpass=TRUE,freq=1,MovingAv=FALSE,
                    W=FALSE,filter=TRUE,threshold=0.05)
points(Data[Peaks[2:13,1],1],Data[Peaks[2:13,1],2],col=2)
#### caculate the peak-to-peak intervals of the TestData
PP_interval <- diff(Peaks[,1])/Fs
Seriesplot.fn(1:length(PP_interval),PP_interval,points=FALSE,xlab="Count",
          ylab="Inter-breath Interval(s)")
#### Smoothing the series by Moving Average or Low pass filter can make the peaks to be found more precise
#### Moving Average and Low pass filter are the argument in the function of find.peaks
#### calculating Moving Average of a series
W <- FS <- 50
Data[,3] <- MovingAverage(Data[,2],W)
Seriesplot.fn(Data[1:2000,1],Data[1:2000,2],points=FALSE,
              xlab="Time(s)",ylab="Respiratory")
lines(Data[1:2000,1],Data[1:2000,3],col=2)
#### Low pass filter of a series
bf <- signal::butter(2, 2/Fs, type="low")
Data[,4] <- signal::filtfilt(bf,Data[,2])
Seriesplot.fn(Data[1:2000,1],Data[1:2000,2],points=FALSE,
              xlab="Time(s)",ylab="Respiratory")
lines(Data[1:2000,1],Data[1:2000,4],col=2)
#### calcualte the multiscale entropy of rawdata
scale_raw <- seq(1,90,2)
MSE <-  MSE(Data$V2[seq(1,100000,2)], tau=scale_raw, m=2, r=0.15, I=40000)
Seriesplot.fn(MSE$tau ,MSE$SampEn,points=TRUE,
              xlab="Scale",ylab="Sample entropy")
#### calcualate the multiscale entropy of peak-to-peak intervals
scale_PP <- 1:10 # setting the scale of entropy
MSE <-  MSE(PP_interval, tau=scale_PP, m=2, r=0.15, I=40000)
Seriesplot.fn(MSE$tau ,MSE$SampEn,points=TRUE,
              xlab="Scale",ylab="Sample entropy")

#### Power Spectral Density (PSD) analysis of the peak-to-peak intervals
LowPSD(PP_interval, plot=TRUE,min=1/64, max=1/2)
#### MultiFractal Detrended Fluctuation Analysis of the peak-to-peak intervals
exponents=seq(3, 9, by=1/4)
scale=2^exponents #Vector of scales
q=-10:10 #q-order of the moment
m=2 #An integer of the polynomial order for the detrending
Result <- MFDFA(PP_interval, scale, m, q)
MFDFAplot.fn(Result,scale,q,model = TRUE)
#### fit the MFDFA result with the extended binomial multifractal model
Coeff <- fit.model(Result$Hq,q)
Coeff
#### Caculate the parameter of multifractal model
Para<- -log(Coeff)/log(2);Para[3]=Para[1]-Para[2]
names(Para)<-c("Hmax","Hmin","i÷H")
Para

#### plot MFDFA results by individual
data("HqData")
PP_Hq <- HqData
filenames <- row.names(PP_Hq)
q=-10:10 # define the range of q-order of the moment
ClassNames <- c(substr(filenames[1:19], start = 1, stop = 3),
                substr(filenames[20:38], start = 1, stop = 5))
Class <- unique(ClassNames) # Obtain group name
col_vec <- rep(NA, nrow(PP_Hq) )
pch_vec <- rep(16, nrow(PP_Hq) )
for( i in 1:length(Class) ) { col_vec[ ClassNames == Class[i] ] <- i }
Individualplot.fn(q,PP_Hq,Name=Class,col=col_vec,pch=pch_vec, xlab="q",ylab="Hurst exponent")
legend("topright", legend=paste0(Class, "(N=", table( ClassNames ), ")"),
      col=1:4, cex=1, lty=1, pch=16)

#### plot the mean and error bar of by MFDFA results group
data("HqData")
PP_Hq <- HqData
filenames <- row.names(PP_Hq)
q <- -10:10 # define the range of q-order of the moment
ClassNames <- c(substr(filenames[1:19], start = 1, stop = 3),
                substr(filenames[20:38], start = 1, stop = 5))
Class <- unique(ClassNames)
for (i in 1:length(q)){
  Data <- GroupComparison.fn(PP_Hq[,i],ClassNames)
  Result_mean_vec <- Data[,"Mean"]
  Result_sd_vec <- Data[,"SE"]
  if( i == 1 ) {
    Result_mean_mat <- Result_mean_vec
    Result_sd_mat <- Result_sd_vec
  } else {
    Result_mean_mat <- rbind(Result_mean_mat, Result_mean_vec)
    Result_sd_mat <- rbind(Result_sd_mat, Result_sd_vec)
 }
}
#### Hurst exponent of q 0-10
Groupplot.fn (q[1:10],Result_mean_mat[1:10,],Class,errorbar = Result_sd_mat[1:10,],
              xRange = NA, yRange = NA, col = NA, pch = rep(16,4), Position = "topright",
              cex.legend = 1, xlab="q",ylab="Hurst exponent",main = "")
#### Hurst exponent of q -10-0
Groupplot.fn (q[11:21],Result_mean_mat[11:21,],Class,errorbar = Result_sd_mat[11:21,],
              xRange = NA, yRange = NA, col = NA, pch = rep(16,4), Position = "topright",
              cex.legend = 1, xlab="q",ylab="Hurst exponent",main = "")


``` 
