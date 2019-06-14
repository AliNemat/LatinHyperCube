Purify_Output <- function(out,x,y) {
  
  fit1<-coefficients(lm (out~x))
  fit2<-coefficients(lm (out~y))
  
  outPure<-vector()
  for (i in 1:10) {
    outPure[i]= out[i]  - ( fit1[[2]]*x[i]+fit1[[1]]) - ( fit2[[2]]*y[i]+fit2[[1]])
  }
  
  return (outPure)
}
