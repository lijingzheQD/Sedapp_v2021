#install rstudioapi/smoothr first
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#preparation################################################################################################
load("rfvmtools.Rdata")
library(Matrix)
source('ChazhiFunction.r')
source('color.bar.R')
source('phi.r')
source('horizontal sum.r')
source('SandFractionProcess.r')
source('DualLithoChenjiFunctionV2.3.R')
source('BestColorFunctionV2.1.R')

#Process the input##################################################################### 
total.subsidence.samples = estimatedd.subsidence.samples#Copy the Subsidence rate
total.subsidence.samples[,2] = estimatedd.subsidence.samples[,2]*0.68*4/3
DS = ChazhiFunction(table = total.subsidence.samples,nx = nx)[,2]/t1  #subsidence rate. The differential subsidence is also available! Unit[L/T].
m <-CreateMesh1D(nx, L)
t=0 #the start time
k = ceiling((t1-t)/dt) #computed times
t1p = t1
kp = ceiling((t1p-t)/dt)#
tperiod=seq(0,(k-1)*dt,length.out = k) #time series
#boundary conditions######################################
bc <- CreateBC1D(m)
bc$left$a = bc_left_a
bc$left$b = bc_left_b
bc$left$c = bc_left_c
bc$right$a = bc_right_a
bc$right$b = bc_right_b
bc$right$c = bc_right_c
#display settings##################################################
x.index=seq(1,nx,length.out = nx) 
x.index.sim=x.index[x.index%%10==0]
huatuxlim = c(0,L)
huatuylim = range(min(topo.samples[,2])- 2000,max(topo.samples[,2]+5))

#initial topography################################################
h0=ChazhiFunction(topo.samples,nx)[,2]
z= BuildCellVariable(h0, bc)
plot(seq(0,L,length.out = nx),h0,ylim = huatuylim)
#sea level settings#################################################
sl=sl.series
maxdepth=max(sl-min(h0))#user-defined, for a better demonstration color.

#Start, Prepare the loop##################################################################################################################
i=1
z.new = z
z.value.hold=h0#for no-erosion
DL = 0; DW = 0
#Create some matrices for storage####################
x.matrix = as.vector(m$cell.centers$x) %o% rep(1,kp)
wd.matrix = matrix(0,nr=nx,nc=kp)
z.matrix = matrix(0,nr=nx,nc=kp)
layer.thickness.matrix = matrix(0,nr=nx,nc=kp)
sand.fraction.matrix = matrix(0,nr=nx,nc=kp)
porosity.matrix = matrix(phi(0.5,0),nr=nx,nc=kp)
binxianx.tvec = rep(0,kp)
binxiany.tvec = rep(0,kp)
slope.fluvial.tvec = rep(0,kp)
time.series.t = rep(0,kp)
DL.series = list()
DW.series = list()
#End, Prepare the loop##################################################################################################################

# start the loop##################################################################################################################
for (i in 1:kp) {
  sl.now = sl.series[i]
  #record the previous z value#####################################################
  if(i>=2){z.value.hold=z.new$value[2 : (nx+1)]  }#for no-erosion
  #compute the water depth value#############################################################
  wd=pmax(sl[i] -z.new$value[2:(nx+1)],0)#wd is updated here!
  #compute the distance from the river mouth point#############################################################
  x.value=m$cell.centers$x
  dx=x.value[2]-x.value[1]
  hekoux=x.value[which.min( abs(   z.value.hold-sl[i]  )   )]
  jlyuanshi=as.vector(abs(x.value-hekoux))
  jlyuanshi[wd==0]=0
  
  #compute start##################################################################################################
  z.newlist = DualLithoChenjiFunction(  m = m
                                        ,bc = bc
                                        ,z = z.new
                                        ,sl.now = sl.now
                                        ,expected.sedimentation.sand = expected.sedimentation.sand
                                        ,expected.sedimentation.mud = expected.sedimentation.mud
                                        ,dt = dt
                                        ,dco = df.dco
                                        ,zhishu = df.zhishu
                                        ,wdscale.sand = df.wdscale.sand
                                        ,wdscale.mud = df.wdscale.mud #the distinction of sand and mud.
                                        ,env.energy = env.energy
                                        ,df.EA = df.EA
                                        ,dep.ero.ratio = dep.ero.ratio
  )
  z.new = z.newlist$z.new
  
  #Start, Compute the sand fraction#############################################################
  layer.thickness.sha.feiling = z.newlist$layer.thickness.sha.feiling[2:(nx+1)]
  # layer.thickness.sha.feiling [layer.thickness.sha.feiling<=0] = 0
  layer.thickness.ni.feiling = z.newlist$layer.thickness.ni.feiling[2:(nx+1)]
  # layer.thickness.ni.feiling [layer.thickness.ni.feiling<=0] = 0
  sand.fraction = layer.thickness.sha.feiling/ (layer.thickness.sha.feiling + layer.thickness.ni.feiling)
  sand.fraction = SandFractionProcess(total.add = layer.thickness.sha.feiling + layer.thickness.ni.feiling
                                      ,sand.fraction = sand.fraction)
    #End, Compute the sand fraction#############################################################
  
  ###########################start, record some data###########################################
  wd.matrix[,i] = wd
  z.matrix[,i] = z.newlist$z.new$value[2 : (nx+1)]
  if(i>=2){layer.thickness.matrix[,i] = z.matrix[,i] - z.matrix[,i-1]}else{
    layer.thickness.matrix[,i] = z.matrix[,i] - h0} #a column is the x series; a row is the layers
  # layer.thickness.matrix[layer.thickness.matrix<0] = 0
  sand.fraction.matrix[,i] = sand.fraction
  time.series.t[i] = t
  #########################end, record some data###########################################
  
  #####################################Start, Isostasy####################################
  layer.thickness.matrix.feifu = layer.thickness.matrix
  layer.thickness.matrix.feifu[layer.thickness.matrix<0]=0
  if(i>=2){DL = (2300-1000)/(3300-1000)*layer.thickness.matrix.feifu[,i]*(1-phi(sand.fraction.matrix[,i],0))
  Dw = 1000/(3300-1000)*(sl[i]-sl[i-1])}else{ DL = 0; DW = 0 }
  DW[wd==0] = 0 ; DW[is.na(DW)] = 0
  #Start, Update some data!
  z.new$value[2 : (nx+1)] = z.new$value[2 : (nx+1)] - DL -DW#The z value has been updated here!
  DL.series[[i]] = DL
  DW.series[[i]] = DW
  if(isostasy){  for(ii in 1:i){z.matrix[,ii] = z.matrix[,ii] - DL -DW }###################zheyaokan##########################
    h0 = h0 - DL -DW}
  #End, Update some data!
  #####################################End, Isostasy####################################
  
  #####################################Start, Tectonic####################################
  #Start, Update some data!
  if(is.array(DS)){
    z.new$value[2 : (nx+1)] = z.new$value[2 : (nx+1)] - DS[,i]*dt#The z value has been updated here!
    for(ii in 1:i){z.matrix[,ii] = z.matrix[,ii] - DS[,i]*dt }
    h0 = h0 - DS[,i]*dt
  }else{
    z.new$value[2 : (nx+1)] = z.new$value[2 : (nx+1)] - DS*dt#The z value has been updated here!
    for(ii in 1:i){z.matrix[,ii] = z.matrix[,ii] - DS*dt }
    h0 = h0 - DS*dt}
  #End, Update some data!
  #####################################End, Tectonic####################################
  
  #################################Start, Compaction process###########################################################
  if(i>=2){for(j in 1:i){
    depth.old = heng.sum(layer.thickness.matrix.feifu[ ,j:i]) - layer.thickness.matrix.feifu[,i]#The distance from layer j to the top layer i-The top layer i
    depth.new = heng.sum(layer.thickness.matrix.feifu[ ,j:i])#The distance from layer j to the top layer i
    phi.old = phi(sand.fraction.matrix[,j] ,depth.old)#Subtract the latest layer
    phi.new = phi(sand.fraction.matrix[,j] ,depth.new)
    compaction.degree = (1 - phi.old)/(1 - phi.new)
    #Start, Update some data!
    layer.thickness.matrix[,j] = compaction.degree * layer.thickness.matrix[,j]#update the thickness!
    z.matrix[,j] = h0 + heng.sum(layer.thickness.matrix[,1:j])
    #End, Update some data!
  }}
  #Start, Update some data!
  z.new$value[2 : (nx+1)] = z.matrix[,i]#The z value has been updated here!
  # print(layer.thickness.matrix[100,1])
  #End, Update some data!
  ###############################End, Compaction process###########################################################
  
  #################################Start, Porosity###########################################################
  for(j in 1:i){
    porosity.matrix[,j] = phi(sand.fraction.matrix[,j],max(z.matrix[,i]) - z.matrix[,j] + overlying.layer.thickness)
  }
  #################################End, Porosity#########################################################
  
  #Start, Compute the shoreline points#############################################################
  # compute the binxianx
  options(warn=-1)#mute the warnings£¡
  binxianx.tvec[i] = max(x.value[wd.matrix[,i]<=0])+(L/nx/2)
  options(warn=0)
  if(binxianx.tvec[i]==-Inf){binxianx.tvec[i]=L}
  # compute the binxiany
  for(j in 1:i){
    tryCatch(#for keeping going on the loop
      expr = { binxiany.tvec[j] = z.matrix[,j][which.min(abs(x.value-binxianx.tvec[j]))]      },
      error = function(e){ print(e) }
    )
    
  }
  #Modify the latest shoreline point
  binxiany.tvec[i] = sl[i]
  myloessmodel1 = loess(x.value~z.matrix[,i],degree = 1 )
  binxianx.tvec[i] = predict(myloessmodel1,binxiany.tvec[i])
  #End, Compute the shoreline points#############################################################
  
  ###########################################Start, Drawing#########################################
  if(i!=-1){
    #draw the original landscape
    plot(x.value
         ,h0
         ,xlim = range(huatuxlim)
         ,ylim = range(huatuylim)
         ,xlab = "x"
         ,ylab = "z"
         ,cex = 1
         ,sub=paste("t = ",t)
    )
    #draw the layers
    if(color_facies1_or_litho2_or_porosity3==2){color.use.matrix = sand.fraction.matrix}
    else if(color_facies1_or_litho2_or_porosity3==3){color.use.matrix = porosity.matrix}
    else{color.use.matrix = wd.matrix}
    for(ii in 1:i){
      points( x.value 
              ,h0 + heng.sum(layer.thickness.matrix[,1:ii])
              ,col = BestColorFunctionV2.1(color.use.matrix,upperLimit = 109,lowerLimit = 0,color.use = color.use)[,ii]
              ,pch = 19
              ,cex = 0.561
      )
      if(ii%%5==0){lines(x.value 
                         ,h0 + heng.sum(layer.thickness.matrix[,1:ii]))}
      abline( h = sl[i],col='blue') #draw the sea level
      lines(template.basin.samples[,1],-estimatedd.subsidence.samples[,2])
    }
  }
  ###########################################End, Drawing##########################################
  
  #¼ÆËãslope#############################################################
  slope.fluvial=(z.new$value[2 : (nx+1)][1]-binxiany.tvec[i])/binxianx.tvec[i]
  print(paste('Fluvial slope is ',slope.fluvial,sep = ''))
  print(paste('Slope angle is ',atan(slope.fluvial)/pi*180,"¡ã",sep = ''))
  #end¼ÆËãslope#############################################################
  
  #output start``````````````````````````````````````````````````````````````````````````````````````````````````````
  if(i%%10==0){
    my_list <- list(  "z.matrix" = z.matrix[,1:i]
                      , "shoreline" = cbind(binxianx.tvec,binxiany.tvec)
                      , "sand.fraction.matrix" = sand.fraction.matrix[,1:i]
                      , "porosity.matrix" = porosity.matrix
                      , "sand.fraction.matrix" = sand.fraction.matrix[,1:i]
                      , "myloessmodel1" = myloessmodel1
                      , "layer.thickness.matrix" = layer.thickness.matrix[,1:i]
                      , "wd.matrix" = wd.matrix
                      , "layer.thickness.matrix" = layer.thickness.matrix[,1:i]
                      , "h0_k" = h0
                      , "DL.series" = DL.series
                      , "DW.series" = DW.series)
    #``````````````````````````````````````````````````````````````````````````````````````````````````````
    t1p=i*2
    wdnow = paste('~/gaobei',t1p,sep = '')
    dir.create(wdnow)
    setwd(wdnow)
    write.csv(my_list$wd.matrix ,file = 'wd.matrix.csv',row.names = FALSE)
    saveRDS(my_list,paste('my_list',t1p,'.rds',sep = ''))
    # write.csv(riverchannel ,file = riverchannelname,row.names = FALSE)
    listnames = names(my_list)
    for (ijk in 1:length(listnames) ){
      listnames[ijk] = paste(listnames[ijk],'.csv',sep = '')
      if(grepl('matrix',listnames[ijk])){write.csv(my_list[[ijk]] ,file = listnames[ijk],row.names = FALSE)}
    }
  }
  
  
  #output end````````````````````````````````````````````````````````````````````````````````````````````````````````
  t=t+dt
  print(paste('time is',t))#;if(t>60){break}
  i=i+1
}
