'Function to plot a color bar
  There are 6 arguments:
    lut: 
      Color vector to be made, e.g., heat.colors(10);
    min: 
      The minimum value indicated by the above vector;
    max: 
      The maximum value indicated by the above vector;
    nticks: 
      The number of scale numbers in the legend;
    ticks: 
      Scale number in legend;
    title: 
      Title.
  The output is to plot.'

color.bar <- function(lut, min, max=-min, nticks=6, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  # dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# example=====================================================================================================
if(0){
  aaa = color.bar(lut = heat.colors(10),min = -1,max = 1)
}
