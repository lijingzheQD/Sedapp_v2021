'this function is used to compute the horizontal sum of a matrix; 
  The input is a matrix;
  The output is a vector.
'
heng.sum = function(mat1){
  if(is.matrix(mat1))
  {hangshu.mat1 = dim(mat1)[1]
    lieshu.mat1 = dim(mat1)[2]
    hengsum = rep(0,hangshu.mat1)
    for (i in 1:hangshu.mat1) {
                                hengsum[i] = sum(mat1[i,])
                              }
  } else  {hengsum = mat1}
  hengsum
}
#example======================================================================================
if(0){
  aaa = matrix(1:9,3,3);aaa
  aab = heng.sum(aaa);aab
}