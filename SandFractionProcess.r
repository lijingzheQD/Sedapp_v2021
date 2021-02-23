'This function is used to optimize the sandstone content.
  There are 2 arguments:
    total.add: 
      total thickness of the formation, vector;
    sand.fraction: 
      the original calculated sandstone content, vector.
  The output is a new vector, that is an updated sand.fraction!'
#####################################################################################################################
SandFractionProcess = function(total.add,sand.fraction){
  sand.fraction[sand.fraction<0] = 0
  sand.fraction[sand.fraction>1] = 1
  sand.fraction[total.add<max(total.add)*0.001] = NA
  
  here.x = 1:length(sand.fraction)
  here.matrix = cbind(here.x,sand.fraction)
  here.matrix = here.matrix[which(!is.na(sand.fraction)),]
  
  library(smoothr)
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  source('ChazhiFunction.r')
  smooth1 = smooth_ksmooth(here.matrix,smoothness = 5)
  sand.fraction.new = ChazhiFunction(table = smooth1
                                     ,nx = length(sand.fraction)
                                     ,start.v = min(cbind(here.x,sand.fraction)[,1])
                                     ,end.v = max(cbind(here.x,sand.fraction)[,1]))
  sand.fraction.new[sand.fraction.new<0] = 0
  sand.fraction.new[sand.fraction.new>1] =1
  return(sand.fraction.new[,2])
}


#example
if(0){
  aaa = c(1,2,3,6,6,6,NA,7,100,9)/10
  plot(aaa,pch=19)
  aab = SandFractionProcess(total.add = 1:10,sand.fraction = aaa)
  points(aab,col=21)
}

