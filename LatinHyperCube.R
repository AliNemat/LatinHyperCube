## Latin hypercube method.

## Libraries
library(scatterplot3d)
#attach(hyperCubeSampledata)
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

# input parameters
Name_vars=c( 'TensionR' , 'Diameter', 'Contract')
binNum=10
inputNum=3 
min_vars=c(1,0,0)
max_vars=c(10,3,12)

## define variables
hyperCubeSample=matrix(1:(binNum*inputNum), nrow = binNum, ncol = inputNum)
colnames(hyperCubeSample) <- Name_vars
bin_sample=matrix(1:(binNum*inputNum), nrow = binNum, ncol = inputNum)

## Find one random sample inside each bin for each input parameter
bin_size=(max_vars-min_vars)/binNum
for (i in 1:binNum) {
  for(j in 1:inputNum) {
    rnd<-runif(1, 0, 1) 
    bin_sample[i,j]=min_vars[j]+(i-1)*bin_size[j]+rnd*bin_size[j]
  }
}

## Generate a matrix which contains in each coloumn, corresponding to one input parameters, 
#  random non-repeated numbers from 1 to 10 
tensionRand=sample(c(1:binNum), size=binNum, replace=FALSE)
diamterRand=sample(c(1:binNum), size=binNum, replace=FALSE)
contractRand=sample(c(1:binNum), size=binNum, replace=FALSE)
sampleExp=cbind (tensionRand, diamterRand, contractRand)

## use random matrix, generated above, to select values already chosen inside each bin. 
for (j in 1:inputNum) {
  for (i in 1:binNum ) {
    a=sampleExp[i,j]
    hyperCubeSample[i,j]=bin_sample[a,j]
  }
}

## plot the resutls
hyperCubeSampledata=as.data.frame(hyperCubeSample)
scatterplot3d(TensionR,Diameter,Contract, pch=16, highlight.3d=TRUE,
              type="h", main="Latin Hypercube Sampling",grid=TRUE) 
addgrids3d(TensionR,Diameter,Contract, grid = c("xy", "xz", "yz"))
