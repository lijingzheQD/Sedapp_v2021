'This function is used to optimize the color setting
  There are four arguments:
    input: 
      Input vector or matrix;
    upperLimit: 
      Maximum. Generally, it is the maximum value of the input vector;
    lowerLimit: 
      Minimum. Generally, it is the minimum value of the input vector;
    color.use: 
      Color plate used. The default is colorRampPalette(heat.colors(3))(101).'
BestColorFunctionV2.1 = function(input,
                             upperLimit = max(input),
                             lowerLimit = min(input),
                             color.use = jet.colors(101)
                            ){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    source('UniformizationFunction.r')
    source('jet_colors.r')
 
  
  input.guiyihua = (input - lowerLimit)/(upperLimit - lowerLimit)
  input.guiyihua = guiyihua(input.guiyihua)
  zcolor = color.use[input.guiyihua*100+1]
  if(is.vector(input)){return(zcolor)}else{
    nrowinput = nrow(input)
    return(matrix(zcolor,nrow = nrowinput))
  }
}

'example===========================
input = 1:100
zcolor = BestColorFunction(input)
plot(input,col=zcolor)
'