## Script to produce figure 4 bifurcation diagrams that show how cover of Coral and Macroalgae at equilibrium change as a function of fishing effort, f, and fisher selectivity, s.

# Assign colors -----
# Variable colors
Ccol <- '3182bd' # blue
Mcol <- '#8c510a' # dark brown


# Set parameter values (using default values) -----
phi_C <- 0.001    # coral open recruitment rate (per year)
phi_I <- 0.0001   # immature mac open recruitment rate (per year)
g_C <- 0.1        # rate at which corals overgrow turf (per year)
g_I <-0.6    # rate at which immature macroalgae overgrow turf (per year)
g_M <- 0.6   # rate at which mature macroalgae overgrow turf (per year)
d_C <- 0.05      # coral natural mortality rate (per year)
d_I <- 0.5    # immature mac natural mortality rate (per year)
d_M <- 0.3    # mature mac natural mortality rate (per year)
r <- 0.5      # local recruitment of immature macroalgae onto turf, produced by the mature stage (per year)
gamma <- 0.5  # scaling constant to slow growth of mature mac over coral vs. free space
omega <- 2     # maturation rate of immature mac into the mature stage
alpha_P <- 10  # herbivory rate for parrotfish 
alpha_U <- 10  # herbivory rate for unicornfish
K_P <- 0.2   # baseline parrotfish carrying capacity 
K_U <- 0.2   # baseline unicornfish carrying capacity 
e <- 0.02  # conversion efficiency of algae by fish
zeta <- 0.2 # scaling constant that reduces herbivory on mature macroalgae
rho <- 4 # parrotfish scalar multiplier



# Run simulation: initially high coral/low macroalgae cover + unicornfish preference (s=0.25) -----

# Set time series to run over
tset <- seq(from = 0, to = 1000, length.out = 10000)

# Set parrotfish and unicornfish growth rates and parrotfish preference
rho <- 4
s <- 0.25 

# create a vector 'f' of values from 0 to 1 
fset <- seq(from = 0, to = 0.15, length.out = 100)

# create holding vectors for equilibrium values of cover for all state variables 
C.star.fbif <- NaN*fset
I.star.fbif <- NaN*fset
M.star.fbif <- NaN*fset
P.star.fbif <- NaN*fset
U.star.fbif <- NaN*fset

# create outer for() loop
# for each value j from 1 to the end of f (so for each value of f)
for(j in 1:length(fset)){
  # assign value for 'f' for each iteration to plug into equations
  f.bif <- fset[j] 
  
  # here we're creating holding vectors for the inner loop
  # create a holding vector for the covers we're calculating at each timestep (for the length of tset) under a disturbed reef scenario 
  C.star.f <- NaN*tset 
  C.star.f[1] <- 0.9  
  I.star.f <- NaN*tset 
  I.star.f[1] <- 0.05
  M.star.f <- NaN*tset 
  M.star.f[1] <- 0.05
  P.star.f <- NaN*tset 
  P.star.f[1] <- 0.1
  U.star.f <- NaN*tset 
  U.star.f[1] <- 0.1
  
  # calculating cover of each state variable for each value of 'f'
  # for each value of 'f'
  for(i in 2:length(tset)){
    # change in time
    dt <- tset[i] - tset[i-1]
    
    # use variable assignment to create state variable objects to plug into calculations 
    C <- C.star.f[i-1]
    I <- I.star.f[i-1]
    M <- M.star.f[i-1]
    P <- P.star.f[i-1]
    U <- U.star.f[i-1]
    
    # calculate change in cover/fish at each timestep
    dC <- dt*( phi_C*(1-M-I-C) + g_C*C*(1-M-I-C) - d_C*C - gamma*g_M*M*C)
    dI <- dt*( phi_I*(1-M-I-C) + r*M*(1-M-I-C) + g_I*I*(1-M-I-C) - d_I*I - omega*I - alpha_U*U*I - alpha_P*P*I)
    dM <- dt*( omega*I + g_M*M*(1-M-I-C) + gamma*g_M*M*C - d_M*M - (zeta*alpha_U)*U*M)
    dP <- dt*( (rho*(e*alpha_P*(1-M-I-C) + e*alpha_P*I))*P *(1- P/K_P) - f.bif*s*P)
    dU <- dt*( (e*alpha_U*(1-M-I-C) + e*alpha_U*I + e*(zeta*alpha_U)*M)*U *(1- U/K_U) - f.bif*(1-s)*U)
    
    
    # calculate cover/fish at that timestep
    C.star.f[i] <- C + dC
    I.star.f[i] <- I + dI	
    M.star.f[i] <- M + dM
    P.star.f[i] <- P + dP
    U.star.f[i] <- U + dU
  }
  
  # final outputs are cover/fish at equilibrium
  C.star.fbif[j] <- C.star.f[length(tset)]
  I.star.fbif[j] <- I.star.f[length(tset)]
  M.star.fbif[j] <- M.star.f[length(tset)]
  P.star.fbif[j] <- P.star.f[length(tset)]
  U.star.fbif[j] <- U.star.f[length(tset)]
}

# Save simulations
C.star.f.HighC <- C.star.fbif
I.star.f.HighC <- I.star.fbif
M.star.f.HighC <- M.star.fbif
P.star.f.HighC <- P.star.fbif
U.star.f.HighC <- U.star.fbif


# Run simulation: initially low coral/high macroalgae cover + unicornfish preference (s=0.25) -----

# Set time series to run over
tset <- seq(from = 0, to = 1000, length.out = 10000)

# Set parrotfish and unicornfish growth rates and parrotfish preference
rho <- 4
s <- 0.25

# create a vector 'f' of values from 0 to 1 
fset <- seq(from = 0, to = 0.15, length.out = 100)

# create holding vectors for equilibrium values of cover for all state variables 
C.star.fbif <- NaN*fset
I.star.fbif <- NaN*fset
M.star.fbif <- NaN*fset
P.star.fbif <- NaN*fset
U.star.fbif <- NaN*fset

# create outer for() loop
# for each value j from 1 to the end of f (so for each value of f)
for(j in 1:length(fset)){
  # assign value for 'f' for each iteration to plug into equations
  f.bif <- fset[j] 
  
  # here we're creating holding vectors for the inner loop
  # create a holding vector for the covers we're calculating at each timestep (for the length of tset) under a disturbed reef scenario 
  C.star.f <- NaN*tset 
  C.star.f[1] <- 0.05  
  I.star.f <- NaN*tset 
  I.star.f[1] <- 0.4
  M.star.f <- NaN*tset 
  M.star.f[1] <- 0.3
  P.star.f <- NaN*tset 
  P.star.f[1] <- 0.1
  U.star.f <- NaN*tset 
  U.star.f[1] <- 0.1
  
  # calculating cover of each state variable for each value of 'f'
  # for each value of 'f'
  for(i in 2:length(tset)){
    # change in time
    dt <- tset[i] - tset[i-1]
    
    # use variable assignment to create state variable objects to plug into calculations 
    C <- C.star.f[i-1]
    I <- I.star.f[i-1]
    M <- M.star.f[i-1]
    P <- P.star.f[i-1]
    U <- U.star.f[i-1]
    
    # calculate change in cover/fish at each timestep
    dC <- dt*( phi_C*(1-M-I-C) + g_C*C*(1-M-I-C) - d_C*C - gamma*g_M*M*C)
    dI <- dt*( phi_I*(1-M-I-C) + r*M*(1-M-I-C) + g_I*I*(1-M-I-C) - d_I*I - omega*I - alpha_U*U*I - alpha_P*P*I)
    dM <- dt*( omega*I + g_M*M*(1-M-I-C) + gamma*g_M*M*C - d_M*M - (zeta*alpha_U)*U*M)
    dP <- dt*( (rho*(e*alpha_P*(1-M-I-C) + e*alpha_P*I))*P *(1- P/K_P) - f.bif*s*P)
    dU <- dt*( (e*alpha_U*(1-M-I-C) + e*alpha_U*I + e*(zeta*alpha_U)*M)*U *(1- U/K_U) - f.bif*(1-s)*U)
    
    # calculate cover/fish at that timestep
    C.star.f[i] <- C + dC
    I.star.f[i] <- I + dI	
    M.star.f[i] <- M + dM
    P.star.f[i] <- P + dP
    U.star.f[i] <- U + dU
  }
  
  # final outputs are cover/fish at equilibrium
  C.star.fbif[j] <- C.star.f[length(tset)]
  I.star.fbif[j] <- I.star.f[length(tset)]
  M.star.fbif[j] <- M.star.f[length(tset)]
  P.star.fbif[j] <- P.star.f[length(tset)]
  U.star.fbif[j] <- U.star.f[length(tset)]
}

# Save simulations
C.star.f.LowC <- C.star.fbif
I.star.f.LowC <- I.star.fbif
M.star.f.LowC <- M.star.fbif
P.star.f.LowC <- P.star.fbif
U.star.f.LowC <- U.star.fbif


# Plot bifurcations for unicornfish selectivity (s=0.25) -----


## Plot bifurcation diagrams (C* as a function of fishing effort at various s_P values)
# Colors
# High i.c.c = darkest
# Low i.c.c. = intermediate
# No i.c.c. = lightest
Ccol.H <- '#08519c'
Ccol.L <- '#6baed6'
Ccol.No <- ''
Icol.H <- '#238b45'
Icol.L <- '#74c476'
Icol.No <- ''
Mcol.H <- '#8c510a'
Mcol.L <- '#cc4c02'
Mcol.No <- '#'

Pcol.H <- '#fd8d3c'
Pcol.L <- '#fed976'
Pcol.No <- ''
Ucol.H <- '#99d8c9'
Ucol.L <- '##edf8b1'
Ucol.No <- ''


# Bifurcation of Proportional Coral Cover at Equilibrium by Total Fishing Effort (f) with two initial conditions: high and low initial coral cover
plot(fset, C.star.f.HighC, 
     type = 'l', 
     lwd = 4, 
     col = Ccol.H, 
     xlab = "", 
     ylab = "", 
     las = 1, 
     ylim=c(0,1), 
     main= "",
     cex.axis=1.5)
lines(fset, C.star.f.LowC, 
      col=Ccol.L, 
      lyt=1, 
      lwd=4)

;
legend("topleft",
       c('High initial coral cover','Low initial coral cover'),
       lty=1,
       lwd=4,
       col=c(Ccol.H, Ccol.L),
       cex=1.5,
       bty="n")


# saved as image with dimensions: 700 width x 550 height

# Bifurcation of Proportional Mature Macroalgae at Equilibrium by Total Fishing Effort (f) with two initial conditions: high and low initial coral cover
plot(fset, M.star.f.HighC, 
     type = 'l', 
     lwd = 4, 
     col = Mcol.H, 
     xlab = "", 
     ylab = "", 
     las = 1, 
     ylim=c(0,1), 
     main= "",
     cex.axis=1.5)
lines(fset, M.star.f.LowC, 
      col=Mcol.L, 
      lyt=1, 
      lwd=4)

;
legend("topleft",
       c('High initial coral cover','Low initial coral cover'),
       lty=1,
       lwd=4,
       col=c(Mcol.H, Mcol.L),
       cex=1.5,
       bty="n")

# Bifurcation of Proportional Immature Macroalgae Cover at Equilibrium by Total Fishing Effort (f) with two initial conditions: high and low initial coral cover
plot(fset, I.star.f.HighC, 
     type = 'l', 
     lwd = 4, 
     col = Icol.H, 
     xlab = "", 
     ylab = "", 
     las = 1, 
     ylim=c(0,1), 
     main= "",
     cex.axis=1.5)
lines(fset, I.star.f.LowC, 
      col=Icol.L, 
      lyt=1, 
      lwd=4)

;
legend("topleft",
       c('High initial coral cover','Low initial coral cover'),
       lty=1,
       lwd=4,
       col=c(Icol.H, Icol.L),
       cex=1.5,
       bty="n")
