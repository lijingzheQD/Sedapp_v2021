'This function is used to normalized a vector.
 for example: c(1,2,3,4,5) ¡ú c(0.00, 0.25, 0.50, 0.75, 1.00)
  The input is a vector, and the output is also a vector.
'
guiyihua = function(variable){
  new = (variable - min(variable))/(max(variable)-min(variable))
  # print(new)
  return(new)
}