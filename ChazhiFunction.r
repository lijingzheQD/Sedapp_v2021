'This function is used to make linear interpolation.
  There are 4 arguments:
      nx: 
        mesh number; 
      table: 
        sample topography.
      start.v: 
        the start value of the first column in the new output table;
      end.v: 
        the end value of the first column in the new output table;
  The output is a new table.
  Note: table.max is the maximum value of the first column of the input table. 
    Similar is the table.min, and they both should satisfy:
    table.max >= end.v > start.v >= table.min. Otherwise, warnings!'
ChazhiFunction <- function( table, nx, start.v = min(table[,1]), end.v=max(table[,1])){
 
  if(start.v < min(table[,1])){warning('start.v < min(table[,1])')}
  if(start.v > end.v){warning('start.v > end.v')}
  if(end.v > max(table[,1])){warning('end.v > max(table[,1])')}
  L = end.v - start.v
  hangshu=dim(table)[1]
  
  samplex=table[,1]
  sampley=table[,2]
  x=seq(start.v,end.v,length.out = nx)
  y=seq(0,0,length.out = nx)
  
  for (i in (hangshu-1):1){
        y[(x<= samplex[i+1])] = (sampley[i+1]-sampley[i])/(samplex[i+1]-samplex[i]) * (x[x<= samplex[i+1] ]-samplex[i]) + sampley[i]
  }
  #These two sentences are used for interpolation
  y[x < samplex[1] ] =sampley[1]
  y[x > samplex[hangshu] ] =sampley[hangshu]
  
  newtable = cbind(x,y)
         
  return(newtable)
}


#example===========================================================================================================================================
if(0){
  topo.samples= cbind(
  c(    0  ,50000,  170000 ),
  c(   105 ,102,  101))
  aaa = ChazhiFunction(table = topo.samples,nx = 100)
  plot(aaa[,2])
}