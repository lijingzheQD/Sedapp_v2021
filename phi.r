#the porosity of sandstone is a function of the depth (meters):
phi.sandstone = function(depth){
  phi = 0.49*exp(depth/(-3700))
}
# the porosity of shale is a function of the depth (meters):
phi.shale = function(depth){
  phi = 0.573*exp(-1.89*depth/1000)
}
# the porosity of silt is a function of the depth (meters):
phi.silt = function(depth){
  phi = 0.659*exp(-1.55*depth/1000)
}  
# the porosity of a mixed layer is a combination of the above functions:
phi = function(sand.fraction,depth){
  phi = sand.fraction*phi.sandstone(depth) + (1-sand.fraction)*phi.shale(depth)
}
#example============================================================================================
if(0){
  depth = c(0:6000)
  plot(depth,phi.sandstone(depth)
       # ,xlim = range(0:6000)
       ,ylim = range(0,0.7)
       ,col = 'red')
  points(depth,phi.shale(depth),col = 'green')
  points(depth,phi.silt(depth),col = 'orange')
  points(depth,phi(0.5,depth),col = 'blue')
}