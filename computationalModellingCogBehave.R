#Random walk model
nreps <- 1e4
nsamples <- 2000

drift <- 0.0 #noninformative stimulus
sdrw <- 0.3
criterion <- 3
t2tsd <- c(0.8,0.0)

latencies <- rep(0,nreps)
responses <- rep(0,nreps)
evidence <- matrix(0,nreps,nsamples+1)

for (i in c(1:nreps)){
  sp = rnorm(1,0,t2tsd[1])
  dr = rnorm(1,0,t2tsd[2])
  evidence[i,] <-
    cumsum(c(sp,rnorm(nsamples,dr,sdrw))) #adds all points in a vector
p<- which(abs(evidence[i,])>criterion)[1] #returns position of values in a vector that meet a criterion
responses[i] <- sign(evidence[i,p]) #Sign returns the sign of a set of numbers
latencies[i] <- p
          
}

#Plot up to 5 random walk paths
tbpn <- min(nreps,5)
plot(1:max(latencies[1:tbpn])+10,type="n",las=1,
     ylim=c(-criterion-.5,criterion+.5),
     ylab = "Evidence", xlab = "decision time")
for (i in c(1:tbpn)){
  lines(evidence[i,1:(latencies[i]-1)])
}
abline(h=c(criterion,-criterion), lty="dashed")

#Plot histogram of latencies

par(mfrow=c(2,1))

toprt <- latencies[responses>0]
topprop <- length(toprt)/nreps
hist(toprt,col="gray",
     xlab="Decision time", xlim = c(0,max(latencies)),
     main = paste("Top responses
                  ( ", as.numeric(topprop),") m=",as.character(signif(mean(toprt),4)),
                  sep =""),las=1)

botrt <- latencies[responses<0]
botprop <- length(botrt)/nreps
hist(botrt, col = "gray",
     xlab = "Decision time",xlim = c(0,max(latencies)),
     main = paste("Bottom responses
                  (",as.numeric(botprop),
                  ") m=", as.character(signif(mean(botrt),4)),
                  sep=""),las = 1)


#Chapter 3 - Least Squar Estimation

#set functions for rmsd and getregpred. These are called by optim later on

#plot data and current predictions
getregpred <- function(parms, data){
  getregpred <- parms["b0"] + parms["b1"] * data[,2]

  #wait with drawing a graph until key is pressed
    par(ask=FALSE)
  
  plot(data[,2], type = "n", las = 1, ylim = c(-2,2),
       xlim = c(-2,2),xlab="X",ylab="Y")
  par(ask = FALSE)
  points (data[,2],data[,1], pch=21, bg="gray")
  lines(data[,2], getregpred, lty = "solid")
  
  return(getregpred)
  }

#obtain current predictions and compute discrepancy

rmsd <- function(parms,data1){
  preds <- getregpred(parms,data1)
  rmsd <- sqrt(sum((preds-data1[,1])^2)/length(preds))
}
#Create a linear model to demonstrate parameter estimation techniques using RMLS


rho <- 0.8
intercept <- .0
nDataPts <- 20

#Generate synthetic data

data <- matrix(0,nDataPts,2)
data[,2] <- rnorm(nDataPts)
data[,1] <- rnorm(nDataPts) * sqrt(1.0-rho^2) + data[
  ,2]*rho + intercept

#Do conventional regression analysis

lm(data[,1] ~ data[,2])

#Assign starting values
startParms <- c(-1.0,.2)
names(startParms) <- c("b1","b0")

#obtain parameter estimates

xout <- optim(startParms, rmsd,data1=data, method = 'SANN')

#Bootstrapping example

# Discrepancy for power forgetting function

powerdiscrep <- function(params,rec,ri){
  
  if (any (params<0)|any(params>1)) return(1e6)
  pow_pred <- params["a"] * (params["b"] * ri + 1)^(-params["c"])
                              return(sqrt( sum((pow_pred-rec)^2)/length(ri) ))
    
}

#Carpenter et al experiment

rec <- c(.93,.88,.86,.66,.47,.34)
ri <- c(0.0035,1,2,7,14,42)

#initialise starting values

sparams <- c(1,0.05,0.7)
names(sparams) <- c('a','b','c')

#obtain best fitting estimate
pout <- optim(sparams,powerdiscrep,rec=rec,ri=ri)
pow_pred <- pout$par["a"] * (pout$par["b"]*c(0:max(ri))+1)^(-pout$par["c"])

#plot data and predictions

x11()
par(cex.axis=1.2,cex.lab=1.4)
par(mar=(c(5,5,3,2)+0.1),las=1)
plot(ri,rec,
     xlab="Retention Interval (Days)",
     ylab = "Proportion Items Retained",
     ylim = c(0.3,1), xlim=c(0,43), xaxt = "n", type = "n")
lines(c(0:max(ri)), pow_pred,lwd=2)
points(ri,rec,pch=21,bg="dark grey", cex=2)
dev <- pow_pred[ri+1]
for(x in c(1:length(ri))){
  lines(c(ri[x],ri[x]), c(dev[x],rec[x]), lwd=1)
}

axis(1,at=c(0:43))

#perform boostrapping analysis

ns <- 55
nbs <- 1000
bsparams <- matrix(NA,nbs,length(sparams))
bspow_pred <- pout$par["a"] * (pout$par["b"]*ri + 1)^(-pout$par["c"])

for ( i in c(1:nbs)){
  recsynth <- vapply(bspow_pred,
                     FUN = function(x) mean(rbinom(ns,1,x)), numeric(1))
  bsparams[i,] <- unlist(optim(pout$par,powerdiscrep,rec=recsynth,ri=ri)$par)
}

#function to plot a histogram

histoplot <- function (x,l4x) {
  hist(x,xlab=l4x,main="",xlim=c(0,1),cex.lab=1.5,cex.axis=1.5)
  lq <- quantile(x,0.025)
  abline(v=lq,lty="dashed",lwd=2)
  uq <- quantile(x,0.975)
  abline(v=uq,lty = "dashed",lwd=2)
  return(c(lq,uq))
}

x11(5,2)
par(mfcol=c(1,3),las=1)
for (i in c(1:dim(bsparams)[2])){
  print(histoplot(bsparams[,i],names(sparams)[i]))
}


#Maximum likelihood parameter estimation (Ch4)

#Using MLE on the GCM model.

source("/Users/oliverclark/Documents/Dropbox/PhD/Modelling/GCMpred.R")

N<- 2*80 #numbr of trials in original study i.e. 2 resp per face

N_A <- round(N*.968) #N_B is implicitly N-N_A

c <- 4
w <- c(0.19, 0.12, 0.25,0.45)

stim <- as.matrix(read.table("/Users/oliverclark/Documents/Dropbox/PhD/Modelling/faceStim.csv", sep = ','))

exemplars <- list(a=stim[1:5], b = stim[6:10,])

preds <- GCMpred(stim[1,], exemplars, c, w)

likelihood <- dbinom(N_A, size = N, prob = preds[1])


#C 4.4 GCM and Maximum Likelihood

GCMprednoisy <- function (probe, exemplars, c, w, sigma, b){
  #Calculate the likelihood of N_A 'A' responses out of N, given parameter c
  #Stim is a single vector representing the stimulus to be categorised
  #exemplars is a list of exemplars in memory and the second list item
  #is the B exemplars in memory
  # c is the scaling parameter and w is a weight vector for each stim dimension
  # for a large number of categories you could use lapply
  
  dist <- list()
  for (ex in exemplars){
    dist[[length(dist)+1]] <- apply(as.array(ex), 1, 
                                    function(x)
                                      sqrt(sum(w*(x-probe)^2)))
  }
  
  sumsim <- unlist(lapply(dist, function(a)
    sum(exp(-c*a))))
  
  r_prob <- c(0,0)
  r_prob[1] <- pnorm(sumsim[1] - sumsim[2] - b, sd =sigma)
  rprob[2] <- 1 - r_prob[1]
return(r_prob)
  
}

library (dfoptim)

GCMutil <- function(theta, stim, exemplars, data, N, retpreds){
  nDat <- length(data)
  dev <- rep(NA,nDat)
  preds <- dev
  
  c <- theta[1]
  w <- theta[2]
  w[2] <- (1-w[1]) * theta[2]
  w[3] <- (1-sum(w[1:2]))
  w[4] <- (1-sum(w[1:3]))
  sigma <- theta[5]
  b <- theta[6]
  for (i in 1:nDat){
    p <- GCMprednoisy(stim[i,], exemplars, c, w, sigma, b)
    dev[i] <- -2*log(dbinom(data[i], size = N, prob = p[1]))
    preds[i] <- p[1]
  }
  
  if (retpreds){
    return(preds)
} else {
  return(sum(dev))
}
}

N <- 2*40

stim <- as.matrix(read.table("/Users/oliverclark/Documents/Dropbox/PhD/Modelling/computational-modelling/codeFromBook/Chapter4/faceStim.csv", sep = ","))

exemplars <- list(a=stim[1:5,], b = stim[6:10,])

data <- scan(file = "/Users/oliverclark/Documents/Dropbox/PhD/Modelling/computational-modelling/codeFromBook/Chapter4/facesDataLearners.txt")
data <- ceiling(data * N)

bestfit <- 1e4

for (w1 in c (0.25,0.5,0.75)){
  for (w2 in c(0.25,0.5,0.75)){
    for (w3 in c(0.25,0.5,0.75)){}
    print(c(w1,w2,w3))
    fitres <- nmkb(par= c(1,w1,w2,w3,1,0.2),
                   fn = function(theta)
                     GCMutil(theta, stim, exemplars, data, N, FALSE),
                   lower = c(0,0,0,0,0,-5),
                   upper = c(10,1,1,1,10,5),
                   control = list(trace=0))
    print(fitres)
    if (fitres$value>bestfit){
      bestres <- fitres
      bestfit <-fitres$value
    }
  }
}


preds <- GCMutil(bestres$par, stim, exemplars, data, N, TRUE)

pdf(file = "GCMfits.pdf", width = 5, height = 5)
plot(preds, data/N, xlab = "Data", ylab = "Predictions")
dev.off()

bestres
theta <- bestres$par
w <- theta[2]
w[2] <- (1-w[1]) * theta[3]
w[3] <- (1-sum(w[1:2]))*theta[4]
w[4] <- (1-sum(w[1:3]))
print(w)

# Chapter 5.4 Aggregating using Vincent Averaging

nsubj <- 30
nobs <- 20
q_p <- c(.1,.3,.5,.7,.9)

shift <- rnorm(nsubj,250,.50)
scale <- rnorm(nsubj,200,50)
shape <- rnorm(nsubj,2,0.25)

params <- rbind(shift,scale,shape)

print (rowMeans(params))

# rows are participants, columns are observations
dat <- apply(params, 2, function(x) rweibull(nobs, shape = x[3], 
                                              scale = x[2]) + x[1])

#Calculate sample quantiles for each participant

kk <- apply(dat,2,function(x) quantile(x,probs=q_p))

## Fitting Quantile averaging

vinq <- rowMeans(kk)

# Fit the shifted Weibull to averaged quantiles

weib_qdev <- function(x,q_emp,q_p){
  if ( any(x <= 0)){
    return(1e6)
  }
  q_pred <- qweibull(q_p, shape = x[3], scale = x[2])+x[1]
                     dev <- sqrt(mean((q_pred-q_emp)^2))
  }

res <- optim(c(225,225,1),
             function(x) weib_qdev(x,vinq,q_p))

print (res)