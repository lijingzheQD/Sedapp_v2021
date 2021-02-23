'This function is to make linear interpolation.
    L:
      the length of x; 
    nx:
      mesh number; 
    table:
      sample topography'
ChazhiFunction <- function(L, nx, table){
  hangshu=dim(table)[1]
  
  samplex=table[,1]
  sampley=table[,2]
  x=seq(0,L,length.out = nx)
  y=seq(0,0,length.out = nx)
  
  for (i in (hangshu-1):1){
        y[(x<= samplex[i+1])] = (sampley[i+1]-sampley[i])/(samplex[i+1]-samplex[i]) * (x[x<= samplex[i+1] ]-samplex[i]) + sampley[i]
  }
  y[x < samplex[1] ] =sampley[1]
  y[x > samplex[hangshu] ] =sampley[hangshu]
         return(y)
}
