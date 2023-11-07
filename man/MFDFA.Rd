\name{MFDFA}
\alias{MFDFA}
\title{MultiFractal Detrended Fluctuation Analysis}
\usage{
MFDFA(tsx, scale, m, q)
}
\arguments{
  \item{tsx}{Univariate time series (must be a vector).}
  \item{scale}{Vector of scales.}
  \item{m}{An integer of the polynomial order for the detrending.}
  \item{q}{q-order of the moment.}
}
\value{
A list of the following elements:
\itemize{
  \item\code{Hq} q-order Hurst exponent.
  \item\code{tau_q} Mass exponent.
  \item\code{hq} Holder exponent.
  \item\code{Dq} singularity dimension.
  \item\code{Fqi} q-order fluctuation function.
  \item\code{line} linear fitting line of fluctuation function.
}
}
\description{
Applies the MultiFractal Detrended Fluctuation Analysis (MFDFA) to time series.
}
\references{
Zhang T, Dong X, Chen C, Wang D, Zhang XD. RespirAnalyzer: an R package for continuous monitoring of respiratory signals.
}
\examples{
data("TestData") # load Data from TestData dataset
Fs <- 50
Peaks <- find.peaks(Data[,2],Fs,lowpass=TRUE,freq=1,MovingAv=FALSE,
                    W=FALSE,filter=TRUE,threshold=0.05)
head(Peaks)
PP_interval <- diff(Peaks$PeakIndex)/Fs
## Computing Multifractal
exponents=seq(3, 9, by=1/4)
scale=2^exponents
q=-10:10
m=2
Result <- MFDFA(PP_interval, scale, m, q)
Coeff <- fit.model(Result$Hq,q)
print(Coeff)
Para<- -log(Coeff)/log(2)
Para[3]=Para[1]-Para[2]
names(Para)<-c("Hmax","Hmin","DeltaH")
Para

PP_Hq <- Result$Hq
PP_hq <- Result$hq
PP_Dq <- Result$Dq
PP_Para <-Para
}
