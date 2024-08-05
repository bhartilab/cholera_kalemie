library(deSolve) #ode
library(Matrix)
library(msm) #rtnorm
library(lubridate)
library(dplyr)

# loading functions needed to run the sampler and providing values for environmental variables
# replace the path toward the script "other_functions_migrN.R" on your machine
source("XXXX/R code/other_functions_migrN.R")

# load weekly aggregared reported cholera cases living in Kalemie "cases.rda" in data
load("XXXX/data/cases.rda") # This file is NOT on the repository, see th  README file.

# load the vacci_step data corresponding to the assumption on vaccine effectiveness (see README.md)
# load(file.path(path.data, "vacci_step3.rda"))
# load(file.path(path.data, "vacci_step3_s1.rda"))
# load(file.path(path.data, "vacci_step3_s2.rda"))
load(file.path(path.data, "vacci_step3_s3.rda"))

# load the functions providing the values of the environmental variables and the data they rely on
source("XXXX/R code/rad_lag2.R")
source("XXXX/R code/rain_lag2.R")
source("XXXX/R code/sst_lag.R")
source("XXXX/R code/chlor_lag.R")

load("XXXX/data/rad_kal.rda")
load("XXXX/data/prcp_sat.rda")
load("XXXX/data/m_rad.rda")
load("XXXX/data/m_prcp.rda")
load("XXXX/data/temp.rda")
load("XXXX/data/chlor.rda")
load("XXXX/data/chlor_30.rda")


################################################
### Preparing the environmental variables    ###
################################################

rad_t <- rad_lag2(model=TRUE)
rain_t <- rain_lag2(lag=0, normalize=FALSE, regular=FALSE, raw=TRUE)
rain_amp_t <- rain_lag2(lag=0, normalize=TRUE, regular=FALSE, raw=TRUE)
sst_t <- sst_lag(raw=TRUE)
chlor_t <- chlor_lag(expand=TRUE)


###################################
### Setting up some parameters  ###
###################################

total.time = 118 #total weeks 
t.obs.seq  = seq(0,total.time,by=1)
N          = 262963
y          = cases

# Change which vacci_step to use based on the assumption on vaccine effectiveness (see README.md)
vacc.ag    = vacci_step3_s1$n

### VACCINATION RATE
eta = stepfun( x = t.obs.seq, 
               y = c(0,vacc.ag/N,0))

suffix = "SUFFIX_FOR THE_NAME_OF_YOUR_FILES"

#################
## ODE setup   ##
#################

R1 = c(-1, 1, 0, 0, 0) # contact s to i 
R2 = c(-1, 0, 1, 0, 0) # vax s to r
R3 = c(0,-1, 1, 0, 0) # recover/remove i to r
R4 = c(1, 0, -1, 0, 0) # immunity waning
R5 = c(0, 0, 0, 1, 0) # fecal contamination
R6 = c(0, 0, 0, -1, 0) # environmental drivers
R7 = c(0, 0, 0, 0, 1) # contact s to i 

R  = cbind(R1, R2, R3, R4, R5, R6, R7)

######################
### MCMC Set-up    ###
######################
iterations   = 300000

burn.in=200000
# Defining the size of the chunks used to regularly save thee samples
write.states = 5000

store.iter   = 1

# Starting values for the sampler
theta.1.curr=c(-0.3430312, 0.3318278, 0.008952414, 0.270325, 0.3623217, 0.05620103, 0.003411574, 1.495769, 1.468861, 0.5324734, 0.4909023)
theta.2.curr=c(4285.66867,75.35729,1.00540)
theta.3.curr=c(0.01736722, -2.466238)

# Adjust based on the assumption on the penalty used to reflect the spatially targeted vaccination in areas with high attack rates
s_capt = 0.8

clock=0

############################
### tuning parameters    ###
############################
c0                    = 1
c1                    = 0.8

######################################
## tuning paramters for theta.1     ##
######################################

theta.1.len             = length(theta.1.curr)
theta.1.batch           = 100
theta.1.tune            = 2.3^2 / theta.1.len
theta.1.jumps           = 0
t.theta.1               = 1
r.theta.1.opt           = 0.23

Sig.theta.1            = rep(1,theta.1.len) 

theta.1.block.matrix    = matrix( 0, nrow = theta.1.batch, ncol = theta.1.len )
theta.1.batch.iter      = 1 
theta.1.block.matrix[theta.1.batch.iter,] = theta.1.curr


theta.1.states = matrix( NA, nrow = write.states, ncol = length(theta.1.curr))


######################################
## tuning paramters for theta.2     ##
######################################

theta.2.len             = length(theta.2.curr) 
theta.2.batch           = 100
theta.2.tune            = 2.3^2 / theta.2.len
theta.2.jumps           = 0
t.theta.2               = 1
r.theta.2.opt           = 0.23

Sig.theta.2             = rep(1,theta.2.len) 
# Sig.theta.2             = 1 

theta.2.block.matrix    = matrix( 0, nrow = theta.2.batch, ncol = theta.2.len )
theta.2.batch.iter      = 1 
theta.2.block.matrix[theta.2.batch.iter,] = theta.2.curr

theta.2.states = matrix( NA, nrow = write.states, ncol = theta.2.len )


# ######################################
# ## tuning paramters for theta.3     ##
# ######################################

theta.3.len             = length(theta.3.curr) 
theta.3.batch           = 100
theta.3.tune            = 2.3^2 / theta.3.len
theta.3.jumps           = 0
t.theta.3               = 1
r.theta.3.opt           = 0.23

theta.3.batch.iter      = 1 

Sig.theta.3            = rep(1,theta.3.len) 

theta.3.block.matrix    = matrix( 0, nrow = theta.3.batch, ncol = theta.3.len )
theta.3.batch.iter      = 1 
theta.3.block.matrix[theta.3.batch.iter,] = theta.3.curr


theta.3.states = matrix( NA, nrow = write.states, ncol = length(theta.3.curr) )



###########################
### prior parameters &  ###
### parameters spaces   ###
###########################

sigma.theta.1 = rep(0.5, length(theta.1.curr))
sigma.theta.2 = c(0.01*N, 0.005*N, 0.5)
sigma.theta.3 = c(0.10, 1)

theta.1.prior.mean = c(0,0,0,0,0,0,1/(5*52),0,1.4,0,0)

theta.1.lower = c( -Inf,-Inf,0,0,0,0,1/(8*52),-Inf,1.15,0,0)
theta.1.upper = c( Inf,Inf,Inf,Inf,Inf,2,1/52,Inf,1.65,Inf,Inf)

theta.2.lower = c(0.01*N,0, 0) 
theta.2.upper = c(0.80*N,10000, 20)

theta.3.prior.mean = c(0.01,5)
theta.3.lower = c(0.0001, -Inf)
theta.3.upper = c(0.3, Inf)

# #################
# ### Solve ode ###
# #################
dt            = 0.01 
X.mesh        = round(seq(0, total.time, by = dt),2)
X.obs.index   = which(X.mesh %in% t.obs.seq)
X.mesh.len    = length(X.mesh)

init.curr     = c(theta.2.curr[1:2], N-theta.2.curr[1]-theta.2.curr[2], 10^theta.2.curr[3], 0)
X.ode.curr    = ode(init.curr,X.mesh,dXdt,theta.1.curr,"rk4")[,-1]
# the function dXdt is defined in "other_functions_migrN.R"
I.ode.curr    = X.ode.curr[-1,5]
I.int.curr    = colSums( matrix(c(0,diff(I.ode.curr)),nrow=1/dt)) 


for(iters in 1:iterations){
  ##############################
  ### Propose theta         ####
  ##############################
  
  theta.1.star = rtnorm( n     = theta.1.len, 
                         mean  = theta.1.curr,
                         sd    = sqrt(theta.1.tune)*sqrt(Sig.theta.1),
                         lower = theta.1.lower,
                         upper = theta.1.upper )

  # solve ODE system for theta.star
  init.curr     = c(theta.2.curr[1:2], N-theta.2.curr[1]-theta.2.curr[2], 10^theta.2.curr[3], 0)
  X.ode.star  = ode(init.curr,X.mesh,dXdt,theta.1.star,"rk4")[,-1]
  I.ode.star     = X.ode.star[-1,5]
  I.int.star     = colSums( matrix(c(0, diff(I.ode.star)),nrow=1/dt)) 
  
  if( all((I.int.star>0) & !is.na(I.int.star)) ){
    #####################################
    ### Evaluate log posterior theta  ###
    ### and proposal densities        ###
    #####################################

    log.post.star = log_post_theta.1( theta.1.input = theta.1.star,
                                      theta.3.input = theta.3.curr,
                                      I.int.input   = I.int.star)
    
    log.post.curr = log_post_theta.1( theta.1.input = theta.1.curr,
                                      theta.3.input = theta.3.curr,
                                      I.int.input   = I.int.curr)
    
    ## Truncated normal proposal density evaluations for MH ratio 
    prop.curr = dtnorm( theta.1.curr, mean = theta.1.star, 
                        sd = sqrt(theta.1.tune)*sqrt(Sig.theta.1), 
                        lower=theta.1.lower, 
                        upper = theta.1.upper,
                        log=TRUE )
    
    prop.star = dtnorm( theta.1.star, mean = theta.1.curr, 
                        sd = sqrt(theta.1.tune)*sqrt(Sig.theta.1), 
                        lower = theta.1.lower, 
                        upper = theta.1.upper, 
                        log   = TRUE ) 
    
    #####################################
    ### Evaluate log MH ratio A.theta  ##
    ### and accept reject step         ##
    #####################################
    A.theta = log.post.star - log.post.curr + sum(prop.curr) - sum(prop.star)
    
    if( log( runif(1) ) < A.theta ){
      
      theta.1.curr  = theta.1.star 
      theta.1.jumps = theta.1.jumps + 1
      
      I.int.curr    = I.int.star 
      X.ode.curr    = X.ode.star 
      
    }
    
  }
  
  # ##############################
  # ### Propose theta.2       ####
  # ##############################
  
  theta.2.star = rtnorm(
                         n    = theta.2.len, 
                         mean  = theta.2.curr, 
                         sd    = sqrt(theta.2.tune)*sqrt(Sig.theta.2),
                         lower = theta.2.lower,
                         upper = theta.2.upper
                    )

  init.star     = c(theta.2.star[1:2], N-theta.2.star[1]-theta.2.star[2], 10^theta.2.curr[3], 0)
  X.ode.star  = ode(init.star,X.mesh,dXdt,theta.1.curr,"rk4")[,-1]
  I.ode.star     = X.ode.star[-1,5]
  I.int.star     = colSums(matrix(c(0, diff(I.ode.star)),nrow=1/dt))
  
  if( all((I.int.star>0) & !is.na(I.int.star)) ){
    #####################################
    ### Evaluate log posterior theta  ###
    ### and proposal densities        ###
    #####################################
    log.post.star = log_post_theta.2( theta.2.input = theta.2.star,
                                      theta.3.input = theta.3.curr,
                                      I.int.input   = I.int.star)
    
    log.post.curr = log_post_theta.2( theta.2.input = theta.2.curr,
                                      theta.3.input = theta.3.curr,
                                      I.int.input   = I.int.curr)
    
    prop.curr = dtnorm( theta.2.curr,
                        mean = theta.2.star, 
                        sd = sqrt(theta.2.tune)*sqrt(Sig.theta.2), 
                        lower=theta.2.lower, 
                        upper = theta.2.upper,
                        log=T )
    
    prop.star = dtnorm( theta.2.star,
                        mean = theta.2.curr, 
                        sd = sqrt(theta.2.tune)*sqrt(Sig.theta.2), 
                        lower = theta.2.lower, 
                        upper = theta.2.upper, 
                        log   = T )

    #####################################
    ### Evaluate log MH ratio A.theta  ##
    ### and accept reject step         ##
    #####################################
    A.theta = log.post.star - log.post.curr + sum(prop.curr) - sum(prop.star)
    
    if( log( runif(1) ) < A.theta ){
      
      theta.2.curr  = theta.2.star 
      theta.2.jumps = theta.2.jumps + 1
      
      X.ode.curr    = X.ode.star 
      I.int.curr    = I.int.star 
      
    }
  }
  
  ##############################
  ### Propose theta.3       ####
  ##############################
  
  theta.3.star = rtnorm( n     = theta.3.len, 
                         mean  = theta.3.curr, 
                         sd    = sqrt(theta.3.tune)*sqrt(Sig.theta.3),
                         lower = theta.3.lower,
                         upper = theta.3.upper )

  #####################################
  ### Evaluate log posterior theta  ###
  ### and proposal densities        ###
  #####################################
  log.post.star = log_post_theta.3( theta.3.input = theta.3.star,
                                    I.int.input   = I.int.curr )
  
  log.post.curr = log_post_theta.3( theta.3.input = theta.3.curr,
                                    I.int.input   = I.int.curr )

  ## Truncated normal proposal density evaluations for MH ratio 
  prop.curr = dtnorm( theta.3.curr,
                      mean  = theta.3.star, 
                      sd    = sqrt(theta.3.tune)*sqrt(Sig.theta.3), 
                      lower = theta.3.lower, 
                      upper = theta.3.upper,
                      log   = T )
    
  prop.star = dtnorm( theta.3.star, 
                      mean  = theta.3.curr, 
                      sd    = sqrt(theta.3.tune)*sqrt(Sig.theta.3), 
                      lower = theta.3.lower, 
                      upper = theta.3.upper,
                      log   = T )
    
    #####################################
    ### Evaluate log MH ratio A.theta  ##
    ### and accept reject step         ##
    #####################################
    A.theta = log.post.star - log.post.curr + sum(prop.curr) - sum(prop.star)
    
    if( log( runif(1) ) < A.theta ){
      
      theta.3.curr  = theta.3.star 
      theta.3.jumps = theta.3.jumps + 1
      
    }
  
  if(iters<=burn.in){
    #########################
    ## Auto tune theta.1  ###
    #########################
    if( theta.1.batch.iter == theta.1.batch ){
      
      theta.1.block.matrix[theta.1.batch.iter,] = as.vector( theta.1.curr) 
      
      t.theta.1        = t.theta.1 + 1
      r.theta.1        = theta.1.jumps / theta.1.batch
      gam.1          = 1/t.theta.1^(c1)
      gam.2          = gam.1*c0
      theta.1.tune     = exp( log(theta.1.tune) + gam.2*(r.theta.1-r.theta.1.opt) )
      
      Sig.theta.1      = Sig.theta.1 + gam.2*( apply(theta.1.block.matrix,2,var) - Sig.theta.1  )
      
      
      theta.1.jumps      = 0
      theta.1.batch.iter = 1 
      
    }else{
      
      theta.1.block.matrix[theta.1.batch.iter,] = as.vector(theta.1.curr) 
      theta.1.batch.iter                      = theta.1.batch.iter + 1
      
    }
    
    # #########################
    # ## Auto tune theta.2  ###
    # #########################
    if( theta.2.batch.iter == theta.2.batch ){
      
      theta.2.block.matrix[theta.2.batch.iter,] = as.vector( theta.2.curr ) 
      
      t.theta.2        = t.theta.2 + 1
      r.theta.2        = theta.2.jumps / theta.2.batch
      gam.1          = 1/t.theta.2^(c1)
      gam.2          = gam.1*c0
      theta.2.tune     = exp( log(theta.2.tune) + gam.2*(r.theta.2-r.theta.2.opt) )
      
      Sig.theta.2      = Sig.theta.2 + gam.2*( apply(theta.2.block.matrix,2,var) - Sig.theta.2  )
      
      
      theta.2.jumps      = 0
      theta.2.batch.iter = 1 
      
    }else{
      
      theta.2.block.matrix[theta.2.batch.iter,] = theta.2.curr 
      theta.2.batch.iter                        = theta.2.batch.iter + 1
      
    }

    #########################
    ## Auto tune theta.3  ###
    #########################
    if( theta.3.batch.iter == theta.3.batch ){
      
      theta.3.block.matrix[theta.3.batch.iter,] = as.vector( theta.3.curr ) 
      
      t.theta.3        = t.theta.3 + 1
      r.theta.3        = theta.3.jumps / theta.3.batch
      gam.1          = 1/t.theta.3^(c1)
      gam.3          = gam.1*c0
      theta.3.tune     = exp( log(theta.3.tune) + gam.3*(r.theta.3-r.theta.3.opt) )
      
      Sig.theta.3      = Sig.theta.3 + gam.3*( apply(theta.3.block.matrix,2,var) - Sig.theta.3  )
      
      
      theta.3.jumps      = 0
      theta.3.batch.iter = 1 

    }else{
      
      theta.3.block.matrix[theta.3.batch.iter,] = theta.3.curr 
      theta.3.batch.iter                        = theta.3.batch.iter + 1
      
    }
  }

  ###########################################################
  ###### Store states if burnin period has passed ###########
  ###########################################################
  # I save the samples 
  if(iters > burn.in){
    if(store.iter == write.states){
      
      theta.1.states[store.iter, ]   = theta.1.curr
      theta.2.states[store.iter, ]   = theta.2.curr
      theta.3.states[store.iter,]   = theta.3.curr

      store.iter                     = 1
      clock=clock+1

      save(
        theta.1.states,  
        file=file.path(path.root, paste0("inst/vacci-model/sample/theta1_", suffix, clock, ".rda")))

      save(
        theta.2.states, 
        file=file.path(path.root, paste0("inst/vacci-model/sample/theta2_", suffix, clock, ".rda")))

      save(
        theta.3.states, 
        file=file.path(path.root, paste0("inst/vacci-model/sample/theta3_", suffix, clock, ".rda")))

      theta.1.states <- matrix(NA, nrow=write.states, ncol = theta.1.len)
      theta.2.states <- matrix( NA, nrow = write.states, ncol = theta.2.len )
      theta.3.states <- matrix( NA, nrow = write.states, ncol = theta.3.len )

    }else{
      theta.1.states[store.iter,]   = theta.1.curr
      theta.2.states[store.iter,]   = theta.2.curr
      theta.3.states[store.iter,]   = theta.3.curr

      store.iter                     = store.iter+1
    }  
  }
  print(iters)
}
