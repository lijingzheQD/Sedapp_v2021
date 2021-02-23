'This function is the solution function of diffusion equation 
  The form of the equation is as follows:
    partial h/partial t = K*partial 2h / partial x2
      K is transport coefficient: K = dco*exp((-(wdscale*jl)**zhishu)/df.EA)+env.energy
  The arguments: 
    m: 
      mesh, see ( m <-CreateMesh1D(nx, L)) £»
    jl: 
      the distance from the eustary£»
    bc: 
      Boundary condition. see( CreateBC1D(m)) , a(du/dx) + bu = c £»
    z: 
      Elevation, the h in the equation, see ( BuildCellVariable(h0, bc), h0 is the topo series after interpolation) £»
    dt: 
      Step size£»
    dco: 
      The pre exponential factor in the transport coefficient, see the expression of K;
    zhishu: 
      Exponent, see the expression of K;
    wdscale: 
      The coefficient of jl, see the expression of K;
    alfac.val: 
      An internal parameter of the diffusion equation. i.e. alfac*partial h/partial t = K*partial 2h / partial x2£»
    env.energy: 
      The environment energy factor, see the expression of K;
    df.EA: 
      The activation engery, see the expression of K;
    dep.ero.ratio: 
      The potential ratio of the deposition rate to the erosion rate£»
    expected.sedimentation: 
      The sediment supply rate. 
'
#Pure chenji function#########################################################################################################################################
  PureComputeChenji = function( m
                                  ,jl
                                  ,bc
                                  ,z
                                  ,env.energy
                                  ,dco
                                  ,wdscale
                                  ,zhishu
                                  ,df.EA
                                  ,dt
                                  ,alfac.val=1
                                  ){
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
      load("rfvmtools.Rdata")
      source('linearinterpolation.r')
   
    
    D = BuildCellProperty(dco*exp((-(wdscale*jl)**zhishu)/df.EA)+env.energy, bc)###for computation use
    Dave = HarmonicMean(D)###for computation use
    Mdiff = DiscretizeDiffusionTerm(Dave)###for computation use
    sbc=DiscretizeBoundaryCondition(bc)###for computation use
    alfac = BuildCellProperty(alfac.val, bc)#Repeat from this step!!!
    S.tran= DiscretizeTransientTerm1D(z, dt, alfac)###for computation use
    M =-Mdiff+sbc$coef.matrix + S.tran$coef.matrix###for computation use
    RHS = sbc$rhs.vector+S.tran$rhs.vector###for computation use
    z.new = SolveEquations(m,M, RHS)
    return(z.new)
  }
#Differentiated chenji function for deposition and erosion######################################################################################################################################################
ChenjiFunctionWithDepEroRatioAndExpectedSedimentation <- function(
  m
  ,jl
  ,bc
  ,z
  ,expected.sedimentation
  ,dt=1
  ,dco= 1e7
  ,zhishu= 2
  ,wdscale=1
  # ,alfac.val = 1
  ,env.energy = 1000
  ,df.EA = 6000000
  ,dep.ero.ratio = 10
) 
{
  alfac.val = 1 
   
  #Start, an initial step for differential deposition and erosion=====================
  z.new = PureComputeChenji(m = m
                              ,jl = jl
                              ,bc = bc
                              ,z = z
                              ,env.energy = env.energy
                              ,dco = dco
                              ,wdscale = wdscale
                              ,zhishu = zhishu
                              ,df.EA = df.EA
                              ,dt = dt
                              ,alfac.val = alfac.val
                              )
  #End, an initial step for differential deposition and erosion=====================
  #Start the differential deposition and erosion==================================
  z.new.value = z.new$value[2:(nx+1)]
  z.original.value = z$value[2:(nx+1)]
  alfac.val=rep(1,nx)
  alfac.val[ z.new.value < z.original.value ]=dep.ero.ratio#better for lands
  #Start to repeat the computation======================================================================
  z.newnew = PureComputeChenji(m = m
                              ,jl = jl
                              ,bc = bc
                              ,z = z
                              ,env.energy = env.energy
                              ,dco = dco
                              ,wdscale = wdscale
                              ,zhishu = zhishu
                              ,df.EA = df.EA
                              ,dt = dt
                              ,alfac.val = alfac.val
                              )
  z.newnew.value = z.newnew$value[2:(nx+1)]
  z.value.creation = (z.newnew.value - z.original.value)
  z.value.creation.sum = sum(z.value.creation)/length(z.newnew$value[2:(nx+1)])
  #get the ratio
    shiji.lilun.ratio = z.value.creation.sum / (expected.sedimentation*dt)
    shiji.lilun.ratio[shiji.lilun.ratio == 0] = 1e-10
  #update the values
     
    z.newnewnew = PureComputeChenji(m = m
                                ,jl = jl
                                ,bc = bc
                                ,z = z
                                ,env.energy = env.energy
                                ,dco = dco/shiji.lilun.ratio
                                ,wdscale = wdscale
                                ,zhishu = zhishu
                                ,df.EA = df.EA
                                ,dt = dt
                                ,alfac.val = alfac.val
                                )
  z.value.creation.sum = sum( z.newnewnew$value[2:(nx+1)] - z.original.value )/length(z.newnew$value[2:(nx+1)])
  
  return(z.newnewnew)
}





#example====================================================================
if(0){
  source('~/R/RFVMToolStartUp.yong.r');source('~/R/Sedapp/success whole function/linearinterpolation.r')
  L <- 170000 #long
  nx <- 1500 #divide
  m <-CreateMesh1D(nx, L)######################
  jl = 0######################
  dt = 0.25
  bc <- CreateBC1D(m)######################
  bc$left$a = 1
  bc$left$b = 0
  bc$left$c = -0.000035
  
  topo.samples= cbind(    c(    0  ,50000,  L ),    c(   105 ,102,  101))
  h0=ChazhiFunction(L,nx,topo.samples)
  z= BuildCellVariable(h0, bc)######################
  
  plot(z$value[2:(nx+1)])
  
  jl = 102 - h0; jl[jl<0] = 0
  
  for (i in 1:10){
    chenji.output = ChenjiFunctionWithDepEroRatioAndExpectedSedimentation(m,jl,bc,z,dt = dt,expected.sedimentation = 0.1,dep.ero.ratio = 168.0000)
    # new.z = chenji.output$z
    z$value = chenji.output$value
    points(z$value[2:(nx+1)],col=2)
    # print(sedimentation.sum)
  }
}




