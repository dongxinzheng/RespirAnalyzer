

#### Functions

### Beginning of functions to be updated
Groupplot.fn <-
  function (x, Average, GroupName, errorbar = NA, xRange = NA,
            yRange = NA, col = NA, pch = NA, Position = "topright", cex.legend = 0.75,
            xlab = "", ylab = "", main = "")
  {
    stopifnot(is.vector(x, mode = "numeric") || length(is.na(x)) ==
                0)
    N_Group <- dim(Average)[2]
    if (is.na(xRange[1]))
      xRange <- range(x)
    if (is.na(col[1]))
      col <- 1:N_Group
    if (is.na(pch[1]))
      pch <- 1:N_Group
    if (is.na(yRange[1]))
      yRange <- range(Average)
    plot(xRange, yRange, type = "n", axes = FALSE, xlab = xlab,
         ylab = ylab, cex.lab=1.6)
    axis(1)
    axis(2, las=2)
    box()
    for (i in 1:N_Group) {
      points(x, Average[, i], col = col[i], pch = pch[i])
      lines(x, Average[, i], col = col[i])
      y1.vec <- Average[, i] - errorbar[, i]
      y2.vec <- Average[, i] + errorbar[, i]
      segments(x, y1.vec, x, y2.vec, col = col[i])
    }
    legend(Position, legend = GroupName, col = col, cex = cex.legend,
           lty = 1, pch = pch)
    title(main = main)
  }

Individualplot.fn <-
  function (x, y, Name = NA, xRange = NA, yRange = NA, col = NA,
            pch = NA, Position = "topright", cex.legend = 0.75, xlab = "",
            ylab = "", main = "")
  {
    stopifnot(is.vector(x, mode = "numeric") || length(is.na(x)) ==
                0)
    N <- dim(y)[1]
    if (is.na(xRange[1]))
      xRange <- range(x)
    if (is.na(col[1]))
      col <- 1:N
    if (is.na(pch[1]))
      pch <- 1:N
    if (is.na(yRange[1]))
      yRange <- range(y)
    plot(xRange, yRange, type = "n", axes = FALSE, xlab = xlab,
         ylab = ylab, cex.lab=1.6)
    axis(1)
    axis(2, las=2)
    box()
    for (i in 1:N) {
      points(x, y[i, ], col = col[i], pch = pch[i])
      lines(x, y[i, ], col = col[i])
    }
    legend(Position, legend = Name, col = col, cex = cex.legend,
           lty = 1, pch = pch)
    title(main = main)
  }


LowPSD <-
  function (series, plot = TRUE, min = 1/8, max = 1/2)
  {
    if (min < 0 || max < 0 || max < min)
      stop("'max' and 'min' must be positive value and max>min")
    A <- series - mean(series)
    n = length(series)
    x = 1:n
    W <- 1 - (2 * x/(n + 1) - 1)^2
    B = A * W
    k = (B[n] - B[1])/(n - 1)
    b = -(B[n] - B[1] * n)/(n - 1)
    y = k * x + b
    C = B - y
    Y = fft(C)
    f = (2:ceiling(n/2))/n
    position = which(f > min & f < max)
    P = abs(Y * Conj(Y)/n)
    x = log10(f[position])
    y = log10(P[position])
    fit <- lm(y ~ x)
    if (plot) {
      plot(x, y, "l", xlab = "log10(frequency)", ylab = "log10(power)", axes=FALSE,cex.lab=1.6)
      axis(1); axis(2, las=2); box()
      cor(x, y)
      abline(fit)
    }
    slope <- coef(fit)
    return(slope)
  }

MFDFAplot.fn <-
  function (Result, scale, q, cex.lab = 1.6, cex.axis = 1.6, col.points = 1,
            col.line = 1, lty = 1, pch = 16, lwd = 2, model = TRUE, cex.legend = 1)
  {
    if (model) {
      Coeff <- fit.model(Result$Hq, q)
      Model <- Result[1:4]
      Model$Hq <- (1/q) * (1 - log(Coeff["a"]^q + Coeff["b"]^q)/log(2))
      Model$Hq[which(q == 0)] <- -log(Coeff["a"] * Coeff["b"])/(2 *
                                                                  log(2))
      Model$tau_q <- -log(Coeff["a"]^q + Coeff["b"]^q)/log(2)
      Model$hq <- diff(Model$tau_q)
      Model$Dq <- q[1:(length(q) - 1)] * Model$hq - Model$tau_q[1:(length(q) -
                                                                     1)]
    }
    #    layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE), heights = c(4, 4))
    #    oldpar <- par(no.readonly = TRUE)
    #    on.exit(par(oldpar))
    #    par(mai = c(0.8, 1, 0.8, 0.4))
    xRange <- range(log2(scale))
    yRange <- range(log2(Result$Fqi))
    plot(xRange, yRange, type = "n", axes = FALSE, xlab = expression("log"[2] *
                                                                       "(Scale)"), ylab = expression("log"[2] * "(F"[q] * ")"),
         cex.lab = cex.lab, cex.axis = cex.axis, main = "A: q-order Fluctuation function")
    Index <- c(1, which(q == 0), which(q == q[length(q)]))
    axis(2, las=2)
    axis(1)
    box()
    for (i in 1:3) {
      k <- Index[i]
      points(log2(scale), log2(Result$Fqi[, k]), col = col.points +
               i - 1, pch = pch)
      lines(log2(scale), Result$line[, k], type = "l", col = col.points +
              i - 1, lwd = lwd)
    }
    legend("bottomright", c(paste("q", "=", q[Index], sep = " ")),
           cex = cex.legend, lwd = c(lwd, lwd, lwd), pch = c(pch,
                                                             pch, pch), bty = "n", col = col.points:(col.points +
                                                                                                       2))
    #    par(mai = c(0.8, 1, 0.8, 0.4))
    plot(q, Result$Hq, col = col.points, axes = F, ylab = expression("h"[q]),
         pch = pch, cex.lab = cex.lab, cex.axis = cex.axis, main = "B: Hurst exponent",
         ylim = range(Result$Hq))
    if (model) {
      lines(q, Model$Hq, col = col.line, lwd = lwd, lty = lty)
      legend("topright", c("Data", "Model"), cex = cex.legend,
             lwd = c(-1, lwd), pch = c(pch, -1), bty = "n", col = c(col.points,
                                                                    col.line))
    }
    axis(1, cex = 4)
    axis(2, cex = 4, las=2)
    box()
    #    par(mai = c(0.8, 1, 0.8, 0.4))
    plot(q, Result$tau_q, col = col.points, axes = F, cex.lab = cex.lab,
         cex.axis = cex.axis, main = "C: Mass exponent", pch = 16,
         ylab = expression(tau[q]))
    if (model) {
      lines(q, Model$tau_q, col = col.line, lwd = lwd, lty = lty)
      legend("bottom", c("Data", "Model"), cex = cex.legend,
             lwd = c(-1, lwd), pch = c(pch, -1), bty = "n", col = c(col.points,
                                                                    col.line))
    }
    axis(1, cex = 4)
    axis(2, cex = 4, las=2)
    box()
    #    par(mai = c(0.8, 1, 0.8, 0.4))
    plot(Result$hq, Result$Dq, col = col.points, axes = F, pch = 16,
         main = "D: Multifractal spectrum", ylab = bquote("f (" ~
                                                            alpha ~ ")"), cex.lab = cex.lab, cex.axis = cex.axis,
         xlab = bquote(~alpha))
    axis(1, cex = 4)
    axis(2, cex = 4, las=2)
    box()
    if (model) {
      lines(Model$hq, Model$Dq, col = col.line, lwd = lwd,
            lty = lty)
      legend("bottom", c("Data", "Model"), cex = cex.legend,
             lwd = c(-1, lwd), pch = c(pch, -1), bty = "n", col = c(col.points,
                                                                    col.line))
    }
  }

fit.model <-
  function (Hq, q)
  {
    Hurst = data.frame(q, Hq)
    names(Hurst) <- c("q", "Hq")
    if (0 %in% q)
      Hurst <- Hurst[-which(q == 0), ]
    a0 = exp(-Hurst$Hq[1]/log(2))
    b0 = exp(-Hurst$Hq[length(Hurst$q)]/log(2))
    m <- nls(Hurst$Hq ~ (1/Hurst$q) * (1 - log(a^Hurst$q + b^Hurst$q)/log(2)),
             start = list(a = a0, b = b0), data = Hurst)
    Goodness <- cor(Hurst$Hq, predict(m))
    A <- summary(m)
    Parameter <- c(A$parameters[, 1], Goodness)
    names(Parameter) <- c("a", "b", "Goodness")
    return(Parameter)
  }

####### End of Functions



if( !require("RespirAnalyzer") ) install.packages( "RespirAnalyzer" )
library(RespirAnalyzer)

# load Data from TestData dataset
data("TestData")

########################### Figure 1 ############################
par( mfrow=c(2,1), mar=c(5.1, 5.1, 3.1, 2.1) )

##### Figure 1A: Raw data of air flow
# Moving Average
W <- Fs <- 50
fit.MA <- MovingAverage(Data[,2],W)

# Low pass filter
bf <- signal::butter(2, 2/Fs, type="low")
fit.LPF <- signal::filtfilt(bf, Data[,2])

nF = 480
#Seriesplot.fn(Data[1:nF,1]-2000, Data[1:nF,2], points=FALSE,
#              xlab="Time in Seconds", ylab="Air Flow")

plot( Data[1:nF,1]-2000, Data[1:nF,2], axes=FALSE, type="n",
      xlab="Time in Seconds",ylab="Air Flow",
      main="A", cex.lab=1.6, cex.main=2)
axis(1); axis(2, las=2); box()
#points(Data[1:nF,1]-2000, Data[1:nF,2], cex=0.5)
lines(Data[1:nF,1]-2000, Data[1:nF,2], cex=2)
lines(Data[1:nF,1]-2000, fit.MA[1:nF], col=2)
lines(Data[1:nF,1]-2000, fit.LPF[1:nF], col=3)
legend( "bottomright", legend=c("Raw data", "Moving average", "Low pass filter"),
        col=1:3, lty=rep(1,3) )

##### Figure 1B: IBI data
Fs=50 ## sampling frequency is 50Hz
Peaks <- find.peaks(Data[,2], Fs, lowpass=TRUE, freq=1, MovingAv=FALSE,
                    W=FALSE, filter=TRUE, threshold=0.05)
#points(Data[Peaks[2:13,1],1],Data[Peaks[2:13,1],2],col=2)
PP_interval <- diff(Peaks[,1])/Fs
plot( 1:length(PP_interval), PP_interval, axes=FALSE, type="n",
      xlab="Air flow cycle series", ylab="IBI in seconds",
      main="B", cex.lab=1.6, cex.main=2)
axis(1); axis(2, las=2); box()
lines( 1:length(PP_interval), PP_interval)


par( mfrow=c(1,2), mar=c(5.1, 5.1, 3.1, 2.1) )
##### Figure 1C: MSE
scale_raw <- seq(1,33)  #scale_raw <- seq(1, 90, 2)
scale_PP <- 1:10
mse.raw <-  MSE(Data$V2[seq(1,100000,2)], tau=scale_raw, m=2, r=0.15, I=40000)
mse.IBI <-  MSE(PP_interval, tau=scale_PP, m=2, r=0.15, I=40000)

plot( mse.raw$tau ,mse.raw$SampEn, axes=FALSE, xlab="Scale",ylab="Sample entropy",
      main="C", cex.lab=1.6)#, cex.main=1)
axis(1, cex=2); axis(2, las=2); box()
lines(mse.raw$tau ,mse.raw$SampEn)

points(mse.IBI$tau ,mse.IBI$SampEn, col=2)
lines(mse.IBI$tau ,mse.IBI$SampEn, col=2)
legend( "bottomright", legend=c("Raw data of air flow", "IBI of air flow"),
        col=1:2, lty=1:2, pch=rep(1,1) )


####  Figure 1D: PSD analysis
LowPSD(PP_interval, plot=TRUE,min=1/64, max=1/2)
title(main="D")

########################### Figure 2 ############################
par( mfrow=c(3,2), mar=c(5.1, 5.1, 3.1, 2.1) )

#par(mfrow=c(3,2))
#### MFDFA: Figure 2A-D
exponents=seq(3, 9, by=1/4)
scale=2^exponents
q=-10:10
m=2
Result <- MFDFA(PP_interval, scale, m, q)
MFDFAplot.fn(Result, scale, q, model = TRUE)

#### fit.model
Coeff <- fit.model(Result$Hq,q)
Coeff
#        a         b  Goodness
#0.4538859 0.7752902 0.9965735
Para<- -log(Coeff)/log(2);Para[3]=Para[1]-Para[2]
names(Para)<-c("Hmax","Hmin","¡÷H")
Para
#     Hmax      Hmin       ¡÷H
#1.1395984 0.3671916 0.7724067

# par(mfrow=c(1,2))
#### Individualplot
data("HqData")
PP_Hq <- HqData
filenames <- row.names(PP_Hq)
q=-10:10
ClassNames <- c(substr(filenames[1:19], start = 1, stop = 3),
                substr(filenames[20:38], start = 1, stop = 5))
Class <- unique(ClassNames)
col_vec <- rep(NA, nrow(PP_Hq) )
pch_vec <- rep(16, nrow(PP_Hq) )
for( i in 1:length(Class) ) { col_vec[ ClassNames == Class[i] ] <- i }
Individualplot.fn(q,PP_Hq,Name=Class,col=col_vec,pch=pch_vec, xlab="q",ylab="Hurst exponent")
legend("topright", legend=paste0(Class, "(N=", table( ClassNames ), ")"),
       col=1:4, cex=1, lty=1, pch=16, bg = "white")
title(main="E: individuals")

#### Groupplot
data("HqData")
PP_Hq <- HqData
filenames <- row.names(PP_Hq)
q <- -10:10
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

Groupplot.fn(q, Result_mean_mat, Class, errorbar = Result_sd_mat,
             xRange = NA, yRange = NA, col = NA, pch = rep(16,4), Position = "topright",
             cex.legend = 1, xlab="q",ylab="Hurst exponent",main = "")
title(main="E: Groups")


