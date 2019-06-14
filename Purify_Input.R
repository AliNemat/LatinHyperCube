Purify_Input <- function(x,y,z) {
  fit1<-coefficients(lm (x~y))
  fit2<-coefficients(lm (x~z))
  xPure<-c(0, 0, 0, 0, 0, 0,0, 0,0,0)
  for (i in 1:10) {
    
    xPure[i]= x[i]  -  ( fit1[[2]]*y[i]+fit1[[1]]) - ( fit2[[2]]*z[i]+fit2[[1]])
  }
  
  return (xPure)
  
}
