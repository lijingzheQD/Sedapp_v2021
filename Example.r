#inputs and settings#################################################################################################
  #mesh settings################################
    #The length of the study area
      L <- 15000 
    #The time step
      dt = 2 
    #The geologic time span
      t1 = 400 
    #Number of Steps
      nx = 300 
  #basin parameters################################
    #The initial topography
      topo.samples= cbind(
        c(    0  ,1890,3900,  L ),
        c(  9 ,0, -86  ,-6.8)) 
    #The Subsidence rate
      estimatedd.subsidence.samples = cbind(
        c(    0,1890.03,2044.67,2216.49,2405.5,2731.96,3250,3402.06,3436.43,4673.5,5910.65,7285.22,8556.7,11271.5,13505.2,14003.4,15000 ),
        c(  0,0,28,441.7,722,1107,1409,1530,1537,1530,1373,1239,1000,530,77,0,0))
    #The final topography
      template.basin.samples = cbind(
        c(    0,1890.03,2044.67,2216.49,2405.5,2731.96,3250,3402.06,3436.43,4673.5,5910.65,7285.22,8556.7,11271.5,13505.2,14003.4,15000 ),
        c(  0,0,28,441.7,722,1107,1409,1479,1486,1479,1373,1239,1135,568,77,0,0))
    #The sea level curve (vector, length = ceiling(t1/dt) )
      sl.series = rep(0,ceiling(t1/dt)) 
    #The parameter Der
      dep.ero.ratio = 1e6 
    #If the isostasy is considered (logical)
      isostasy = FALSE 
    #The overlying layer thickness
      overlying.layer.thickness = 0 
    #The preexponential factor (Alpha in Eq. (9) in https://gmd.copernicus.org/preprints/gmd-2020-256/)
      df.dco = 168e6 
    #The characteristic energy (E in Eq. (9))
      df.EA = 968e2 
    #The exponent scale (Eta in Eq. (4))
      df.zhishu=2 
    #The adjust factor (Epsilon in Eq. (9))
      env.energy = 100 
    #The sediment supply rate of sand
      expected.sedimentation.sand = 0.135 
    #The sediment supply rate of mud
      expected.sedimentation.mud = expected.sedimentation.sand * 5 
    #The distance scale for sand (c in Eq. (9))
      df.wdscale.sand = 0.1 
    #The distance scale for mud (c in Eq. (9))
      df.wdscale.mud = 0.75*df.wdscale.sand 
    #The property to be displayed (choose from 1,2 and 3).
      color_facies1_or_litho2_or_porosity3 = 1 
    #Boundary conditions
      bc_left_a = 0
      bc_left_b = 1
      bc_left_c = 9
      bc_right_a = 1
      bc_right_b = 0
      bc_right_c = 0
    #The color palette (The default is rev(topo.colors(101)))
      color.use = rev(topo.colors(101))

# Source the code
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
  source('Sedapp to be called.R')
  