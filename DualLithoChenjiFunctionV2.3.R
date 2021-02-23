'This chenji function for dual lithologies 
  The arguments: 
    m: 
      mesh, see ( m <-CreateMesh1D(nx, L)) £»
    bc: 
      Boundary condition. see( CreateBC1D(m)) , a(du/dx) + bu = c £»
    z: 
      Elevation, the h in the equation, see ( BuildCellVariable(h0, bc), h0 is the topo series after interpolation) £»
    sl.now:
      the current sea level;
    expected.sedimentation.sand:
      the expected sedimentation of sand;
    expected.sedimentation.mud:
      the expected sedimentation of mud;
    dt: 
      Step size£»
    dco: 
      The pre exponential factor in the transport coefficient, see the expression of K;
    zhishu: 
      Exponent, see the expression of K;
    wdscale.sand: 
      The coefficient of jl for sand, see the expression of K;
    wdscale.mud: 
      The coefficient of jl for mud, see the expression of K;
    env.energy: 
      The environment energy factor, see the expression of K;
    df.EA: 
      The activation engery, see the expression of K;
    dep.ero.ratio: 
      The potential ratio of the deposition rate to the erosion rate.
'
DualLithoChenjiFunction = function(    m 
                                      ,bc 
                                      ,z 
                                      ,sl.now
                                      ,expected.sedimentation.sand 
                                      ,expected.sedimentation.mud 
                                      ,dt 
                                      ,dco 
                                      ,zhishu 
                                      ,wdscale.sand 
                                      ,wdscale.mud 
                                      ,env.energy 
                                      ,df.EA 
                                      ,dep.ero.ratio = 100
                                      ){
  #########################################################################################################################################################################
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    source('SandFractionProcess.r')
    source('ChenjiFunctionWithDepEroRatioAndExpectedSedimentationV1.r')
  
  #compute the next water depth and jl according to the topography and sea level==========================================================================================
    wd=pmax(sl.now -z.new$value[2:(nx+1)],0)#wd is a series greater than 0. update here!!!
    #compute the distance from the river mouth point---------------------------
      x.value=m$cell.centers$x
      hekoux=x.value[which.min( abs(   z.new$value[2:(nx+1)]-sl[i]  )   )]
      jlyuanshi=as.vector(abs(x.value-hekoux))
      jl = jlyuanshi;  jl[wd==0]=0 
  #compute the sand first==================================
    z.new.sand = ChenjiFunctionWithDepEroRatioAndExpectedSedimentation(   m = m
                                    ,jl = jl
                                    ,bc = bc
                                    ,z = z
                                    ,expected.sedimentation = expected.sedimentation.sand
                                    ,dt = dt
                                    ,dco = dco
                                    ,zhishu = zhishu
                                    ,wdscale = wdscale.sand
                                    ,env.energy = env.energy
                                    ,df.EA = df.EA
                                    ,dep.ero.ratio = dep.ero.ratio)
    # points(m$cell.centers$x
    #        ,z.new.sand$value[2:(nx+1)],col='yellow')
  #then compute the mud==================================
    z.new.mud = ChenjiFunctionWithDepEroRatioAndExpectedSedimentation(   m = m
                                 ,jl = jl
                                 ,bc = bc
                                 ,z = z.new.sand
                                 ,expected.sedimentation = expected.sedimentation.mud
                                 ,dt = dt
                                 ,dco = dco
                                 ,zhishu = zhishu
                                 ,wdscale = wdscale.mud #distinction of sand and mud
                                 ,env.energy = env.energy
                                 ,df.EA = df.EA
                                 ,dep.ero.ratio = dep.ero.ratio)
    # points(m$cell.centers$x
    #        ,z.new.mud$value[2:(nx+1)],col='gray')
  #combine the results of sand and mud
    #amount of sand
      sand.add = z.new.sand$value - z$value
      # sand.add[sand.add<0] = 0
    #amout of mud
      mud.add = z.new.mud$value - z.new.sand$value
      # mud.add[mud.add<0] = 0
    #compute the sand fraction
      total.add = sand.add + mud.add
      sand.fraction.original = sand.add/(sand.add + mud.add)
      sand.fraction.new = SandFractionProcess(total.add = total.add
                                              ,sand.fraction = sand.fraction.original)
    #update for output
      z.new = z
      z.new$value = z$value + sand.add + mud.add
      # points(m$cell.centers$x
      #        ,z.new$value[2:(nx+1)],col='black',cex=0.5)
      layer.thickness.sha.feiling = sand.add; layer.thickness.sha.feiling = pmax(layer.thickness.sha.feiling,0)
      layer.thickness.ni.feiling = mud.add; layer.thickness.ni.feiling = pmax(layer.thickness.ni.feiling,0)
      output1 = list('z.new' = z.new
                     ,'sand.fraction' = sand.fraction.new
                     ,'wd' = wd
                     ,'layer.thickness.sha.feiling' = layer.thickness.sha.feiling
                     ,'layer.thickness.ni.feiling' = layer.thickness.ni.feiling)
    return(output1)
}
#example==============================================================================================================
if(0){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  #################################################################################################
  load("rfvmtools.Rdata")
  library(Matrix)
  source('ChazhiFunction.r')
  source('color.bar.R')
  source('phi.r')
  source('horizontal sum.r')
  source('jet_colors.r')
  source('SandFractionProcess.r')
  source('DualLithoChenjiFunctionV2.3.R')
  source('BestColorFunctionV2.1.R')
  # library(export)
  #inputs and settings#################################################################################################
  #mesh settings
  L <- 15000 #long
  dt = 2
  t1 = 400
  nx = 300
  sed.supply = 0.0003795
  topo.samples= cbind(
    c(    0  ,1890,3900,  L ),
    c(  9 ,0, -86  ,-6.8))
  estimatedd.subsidence.samples = cbind(
    c(    0,1890.03,2044.67,2216.49,2405.5,2731.96,3250,3402.06,3436.43,4673.5,5910.65,7285.22,8556.7,11271.5,13505.2,14003.4,15000 ),
    c(  0,0,28,441.7,722,1107,1409,1530,1537,1530,1373,1239,1000,530,77,0,0))
  template.basin.samples = cbind(
    c(    0,1890.03,2044.67,2216.49,2405.5,2731.96,3250,3402.06,3436.43,4673.5,5910.65,7285.22,8556.7,11271.5,13505.2,14003.4,15000 ),
    c(  0,0,28,441.7,722,1107,1409,1479,1486,1479,1373,1239,1135,568,77,0,0))
  total.subsidence.samples = estimatedd.subsidence.samples; total.subsidence.samples[,2] = estimatedd.subsidence.samples[,2]*0.68*4/3
  sl.series = rep(0,ceiling(t1/dt)) #+ sin(1:(ceiling(t1/dt)))*0.5
  DS = ChazhiFunction(table = total.subsidence.samples,nx = nx)[,2]/t1  #
  df.dco = 168e6 #
  df.EA = 968e2 #
  color_facies1_or_litho2_or_porosity3 = 1
  dep.ero.ratio = 1e6
  isostasy = FALSE
  overlying.layer.thickness = 0
  m <-CreateMesh1D(nx, L)
  #time frame
  t=0 #the start time
  k = ceiling((t1-t)/dt) #computed times
  
  
  t1p = t1##################################################################################################
  kp = ceiling((t1p-t)/dt)##################################################################################################
  
  tperiod=seq(0,(k-1)*dt,length.out = k) #time series
  #boundary conditions######################################
  bc <- CreateBC1D(m)
  bc$left$a = 0
  bc$left$b = 1
  bc$left$c = 9
  bc$right$a = 1
  bc$right$b = 0
  bc$right$c = 0
  #display settings##################################################
  ps=2.0 #set pause time
  colorpiandu=1 #bias of color, normal=1.
  x.index=seq(1,nx,length.out = nx) 
  x.index.sim=x.index[x.index%%10==0]
  huatuxlim = c(0,L)
  huatuylim = range(min(topo.samples[,2])- 2000,max(topo.samples[,2]+5))
  #df parameters######################################################
  df.zhishu=2 #
  env.energy = 100 #
  expected.sedimentation.sand = 0.135 #
  expected.sedimentation.mud = expected.sedimentation.sand * 5 #
  df.wdscale.sand = 0.1 #
  df.wdscale.mud = 0.75*df.wdscale.sand
  #initial topography################################################
  h0=ChazhiFunction(topo.samples,nx)[,2]
  z= BuildCellVariable(h0, bc)
  plot(seq(0,L,length.out = nx),h0,ylim = huatuylim)
  #sea level settings#################################################
  sl=sl.series
  maxdepth=max(sl-min(h0))#
  
  
  #prepare the loop##################################################################################################################
  i=1
  z.new = z
  z.value.hold=h0#for no-erosion
  DL = 0; DW = 0
  #create some matrices for storage purpose####################
  x.matrix = as.vector(m$cell.centers$x) %o% rep(1,kp)
  wd.matrix = matrix(0,nr=nx,nc=kp)
  z.matrix = matrix(0,nr=nx,nc=kp)
  layer.thickness.matrix = matrix(0,nr=nx,nc=kp)
  sand.fraction.matrix = matrix(0,nr=nx,nc=kp)
  porosity.matrix = matrix(phi(0.5,0),nr=nx,nc=kp)
  binxianx.tvec = rep(0,kp)
  binxiany.tvec = rep(0,kp)
  slope.fluvial.tvec = rep(0,kp)
  zcolor.matrix = matrix(0,nr=nx,nc=kp)
  zcolor.facies = matrix(0,nr=nx,nc=kp)
  zcolor.litho = matrix(0,nr=nx,nc=kp)
  zcolor.porosity = matrix(0,nr = nx,nc=kp)
  time.series.t = rep(0,kp)
  DL.series = list()
  DW.series = list()
  # thickness.matrix = matrix(0,nr=nx,nc=k)
  #end the preparation of the loop##################################################################################################################
  
  # start the loop##################################################################################################################
    i = 1  
    sl.now = sl.series[i]
    #record the previous z value#####################################################
    if(i>=2){z.value.hold=z.new$value[2 : (nx+1)]  }#for no-erosion
    #compute the water depth value#############################################################
    wd=pmax(sl[i] -z.new$value[2:(nx+1)],0)#
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
                                          ,wdscale.mud = df.wdscale.mud #
                                          ,env.energy = env.energy
                                          ,df.EA = df.EA
                                          ,dep.ero.ratio = dep.ero.ratio
    )
}
