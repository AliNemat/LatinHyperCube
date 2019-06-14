########################### sensitivity analysis ####################################################

## Libraries
setwd("C:/Users/Shixin Xu/Ali_Files/UCR/Research/Cell height/R")
library(scatterplot3d)
library(ggplot2)
library(gridExtra)

library(ggpubr)
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
source('pcor.R')
source ('Purify_Input.R')
source ('Purify_Output.R')

#######################Latin hypercube samping method ############################################
# input parameters
Name_vars=c( 'TensionR' , 'Diameter', 'Contract')
binNum=10
inputNum=3
min_vars=c(1,0,0)
max_vars=c(10,2,12)
createNewSample=FALSE
## define variables
hyperCubeSample=matrix(1:(binNum*inputNum), nrow = binNum, ncol = inputNum)
bin_sample     =matrix(1:(binNum*inputNum), nrow = binNum, ncol = inputNum)
colnames(hyperCubeSample) <- Name_vars


if (createNewSample) {
  ## Find one random sample inside each bin for each input parameter
  bin_size=(max_vars-min_vars)/binNum
  for(j in 1:inputNum) {
    for (i in 1:binNum) {
      rnd<-runif(1, 0, 1) 
      bin_sample[i,j]=min_vars[j]+(i-1)*bin_size[j]+rnd*bin_size[j]
    }
  }
  
  ## Put randomly beans for each input parameter next to each other to create one set of beans for simulation
  tensionRand=sample(c(1:binNum), size=binNum, replace=FALSE)
  diamterRand=sample(c(1:binNum), size=binNum, replace=FALSE)
  contractRand=sample(c(1:binNum), size=binNum, replace=FALSE)
  sampleExp=cbind (tensionRand, diamterRand, contractRand)
  
  ## replace beans index with actual value of that paramerer inside the bean 
  for (j in 1:inputNum) {
    for (i in 1:binNum ) {
      hyperCubeSample[i,j]=bin_sample[sampleExp[i,j],j]
    }
  }
}else {
  
  ############################# load the hypercube sample that we used for simulation ###################
  attach (HyperCube_parameters_output_March4th2019)
  inputs<-cbind ( TensionR[2:11],D_perecent_mapped[2:11],Contract[2:11]/6)
  for (j in 1:inputNum) {
    for (i in 1:binNum ) {
      hyperCubeSample[i,j]=inputs [i,j]
    }
  }
  
}

## plot and write the latin hypercube sampling results 
pdf("HypercubeMarch19th2019.pdf",width=12, height=8, useDingbats=FALSE, bg="white")

hyperCubeSampledata=as.data.frame(hyperCubeSample)

scatterplot3d(hyperCubeSampledata$TensionR,
              hyperCubeSampledata$Diameter,
              hyperCubeSampledata$Contract, 
              pch=16, highlight.3d=TRUE,
              type="h", main="Latin Hypercube Sampling",grid=TRUE) 
addgrids3d(hyperCubeSampledata$TensionR,hyperCubeSampledata$Diameter,hyperCubeSampledata$Contract
           , grid = c("xy", "xz", "yz"))

dev.off()

write.table(hyperCubeSampledata, file = "HyperCubeResultMarch11th2019.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)


####################Load output results from external code ############################################
outputNum=3
Name_Outs=c( 'curvature', 'height','nucLoc')
attach (HyperCube_parameters_output_March4th2019)
outputs<-cbind ( Curvature[2:11],height[2:11],Nucleus[2:11])

################################ Merge input and output data #########################################
inputsAndOutputs<-cbind ( hyperCubeSampledata,outputs)
Names<-c ( Name_vars,Name_Outs )
colnames(inputsAndOutputs) <-Names
inputsAndOutputsData=as.data.frame(inputsAndOutputs)

## just give a short name to make earier wrking with it.
dataT=inputsAndOutputsData

################## Partial correlation analysis of input and output data  #############################

# Removing linear correlation of input parameters on the other input parameters
tensionRPure<-vector()
diameterPure<-vector()
contractPure<-vector()

#tensionRPure<-dataT$TensionR #Purify_Input (dataT$TensionR, dataT$Diameter, dataT$Contract )
#diameterPure<-dataT$Diameter #Purify_Input (dataT$Diameter, dataT$TensionR, dataT$Contract )
#contractPure<-dataT$Contract #Purify_Input (dataT$Contract, dataT$TensionR, dataT$Diameter )

tensionRPure<-Purify_Input (dataT$TensionR, dataT$Diameter, dataT$Contract )
diameterPure<-Purify_Input (dataT$Diameter, dataT$TensionR, dataT$Contract )
contractPure<-Purify_Input (dataT$Contract, dataT$TensionR, dataT$Diameter )

inputPure=cbind (tensionRPure, diameterPure, contractPure)

# Removing linear correlations of output parameter, Curvature, with respect to different inputs 

curvatureT<-vector()
curvatureD<-vector()
curvatureC<-vector()
#curvatureT<-Purify_Output (dataT$curvature, dataT$Diameter, dataT$Contract)
#curvatureD<-Purify_Output (dataT$curvature, dataT$TensionR, dataT$Contract)
#curvatureC<-Purify_Output (dataT$curvature, dataT$TensionR, dataT$Diameter)

curvatureT<-Purify_Output (dataT$curvature, dataT$Diameter, dataT$Contract)
curvatureD<-Purify_Output (dataT$curvature, dataT$TensionR, dataT$Contract)
curvatureC<-Purify_Output (dataT$curvature, dataT$TensionR, dataT$Diameter)

## Removing linear correlations of output parameter, nucLoc, with respect to different inputs 

nucLocT<-vector()
nucLocD<-vector()
nucLocC<-vector()
nucLocT<-Purify_Output (dataT$nucLoc, dataT$Diameter, dataT$Contract)
nucLocD<-Purify_Output (dataT$nucLoc, dataT$TensionR, dataT$Contract)
nucLocC<-Purify_Output (dataT$nucLoc, dataT$TensionR, dataT$Diameter)


## Removing linear correlations of output parameter, height, with respect to different inputs 
heightT<-vector()
heightD<-vector()
heightC<-vector()
heightT<-Purify_Output (dataT$height, dataT$Diameter, dataT$Contract)
heightD<-Purify_Output (dataT$height, dataT$TensionR, dataT$Contract)
heightC<-Purify_Output (dataT$height, dataT$TensionR, dataT$Diameter)


############## merge purified input and output file as a data frame ################################

InputAndOutputPure=cbind (tensionRPure, diameterPure, contractPure, 
                          curvatureT,   curvatureD,   curvatureC,
                          nucLocT,      nucLocD,      nucLocC,      
                          heightT,      heightD,      heightC)

dataPure=as.data.frame(InputAndOutputPure)

pdf("plotSensitivity_June12th.pdf",width=12, height=8, useDingbats=FALSE, bg="white")
par(mfrow=c(3,3))

plot1<-ggscatter(dataPure, x ="tensionRPure", y ="curvatureT", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "tensionRPure", ylab = "curvatureT")

plot2<-ggscatter(dataPure, x ="diameterPure", y ="curvatureD", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "diameterPure", ylab = "curvatureD")


plot3<-ggscatter(dataPure, x ="contractPure", y ="curvatureC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "contractPure", ylab = "curvatureC")


################################################################
plot4<-ggscatter(dataPure, x ="tensionRPure", y ="nucLocT", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "tensionRPure", ylab = "nucLocT")

plot5<-ggscatter(dataPure, x ="diameterPure", y ="nucLocD", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "diameterPure", ylab = "nucLocD")


plot6<-ggscatter(dataPure, x ="contractPure", y ="nucLocC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "contractPure", ylab = "nucLocC")
#########################################################################

plot7<-ggscatter(dataPure, x ="tensionRPure", y ="heightT", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "tensionRPure", ylab = "heightT")

plot8<-ggscatter(dataPure, x ="diameterPure", y ="heightD", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "diameterPure", ylab = "heightD")


plot9<-ggscatter(dataPure, x ="contractPure", y ="heightC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "contractPure", ylab = "heightC")
###################################################################################



grid.arrange(plot1, plot2,plot3,plot4, plot5, plot6, plot7, plot8,plot9,ncol=3, nrow = 3)


dev.off()


fitCT<-coefficients(lm (dataPure$curvatureT ~dataPure$tensionRPure))
fitCD<-coefficients(lm (dataPure$curvatureD ~dataPure$diameterPure))
fitCC<-coefficients(lm (dataPure$curvatureC ~dataPure$contractPure))
curvatureCoef=cbind (fitCT[[2]] ,fitCD[[2]],fitCC[[2]])

fitHT<-coefficients(lm (dataPure$heightT ~dataPure$tensionRPure))
fitHD<-coefficients(lm (dataPure$heightD ~dataPure$diameterPure))
fitHC<-coefficients(lm (dataPure$heightC ~dataPure$contractPure))
heightCoef=cbind (fitHT[[2]] ,fitHD[[2]],fitHC[[2]])

fitNT<-coefficients(lm (dataPure$nucLocT ~dataPure$tensionRPure))
fitND<-coefficients(lm (dataPure$nucLocD ~dataPure$diameterPure))
fitNC<-coefficients(lm (dataPure$nucLocC ~dataPure$contractPure))
nucLocCoef=cbind (fitNT[[2]] ,fitND[[2]],fitNC[[2]])


par(mfrow=c(1,3))
barplot (curvatureCoef, names.arg=c("ECM Tension", "N Diameter", "Contraction"))
barplot (heightCoef, names.arg=c("ECM Tension", "N Diameter", "Contract"))
barplot (nucLocCoef, names.arg=c("ECM Tension", "N Diameter", "Contract"))