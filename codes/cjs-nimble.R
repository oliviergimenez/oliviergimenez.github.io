#---------- SIMULATE DATA FROM CJS MODEL (from Kery and Schaub book)

# Define function to simulate a capture-history (CH) matrix 
simul.cjs <- function(PHI, P, marked){
	n.occasions <- dim(PHI)[2] + 1
	CH <- matrix(0, ncol = n.occasions, nrow = sum(marked)) 
	# Define a vector with the occasion of marking 
	mark.occ <- rep(1:length(marked), marked[1:length(marked)]) 
	# Fill the CH matrix 
	for (i in 1:sum(marked)){ 
		CH[i, mark.occ[i]] <- 1 # Write an 1 at the release occasion 
		if (mark.occ[i]==n.occasions) next
		for (t in (mark.occ[i]+1):n.occasions){ 
			# Bernoulli trial: has individual survived occasion? 
			sur <- rbinom(1, 1, PHI[i,t-1]) 
			if (sur==0) break # If dead, move to next individual 
			# Bernoulli trial: has individual been recaptured? 
			rp <- rbinom(1, 1, P[i,t-1]) 
			if (rp==1) CH[i,t] <- 1 
			} # t loop 
		} # i loop 
return(CH) 
} 

# parameter values 
k <- 6 # number of capture occasions 
marked <- rep(50, k-1) # annual number of newly marked individuals 
nind <- sum(marked)
phi <- rep(0.65, k-1) # survival 
p <- rep(0.4, k-1) # detection

# matrices with survival and recapture probabilities 
PHI <- matrix(phi, ncol = k-1, nrow = sum(marked)) 
P <- matrix(p, ncol = k-1, nrow = sum(marked)) 

# simulate 
y <- simul.cjs(PHI, P, marked) 

#---------- FIT CJS MODEL TO SIMULATED DATA USING NIMBLE

# load Nimble (see Nimble manual for installation)
library(nimble)

# input code (identical to WinBUGS/JAGS code)
dipper <- modelCode({
phi ~ dunif(0,1)
p ~ dunif(0,1)
for(i in 1:nind) {
	x[i, first[i]] <- 1
	for(t in (first[i]+1):k) {
		mu_x[i,t] <- phi * x[i,t-1]
		mu_y[i,t] <- p * x[i,t]
		x[i,t] ~ dbern(mu_x[i,t])
		y[i,t] ~ dbern(mu_y[i,t])
							} # t loop
				} # i loop
})

# constant values (see Nimble manual)
get.first <- function(x) min(which(x!=0)) # get occasion of marking 
first <- apply(y, 1, get.first) 
constants <- list(k=k,nind=nind,first=first)

# data
data <- list(y = y)

# initial values
x_init <- array(1, dim=c(nind,k))
for(i in 1:nind) {
    x_init[i, 1:first[i]] <- NA
}
inits <- list(phi = 0.6, p = 0.4, x = x_init)

# compile model (see Nimble manual)
cjs <- nimbleModel(dipper, name='cjs', constants=constants, data=data, inits=inits)
# cjs$getNodeNames()
Ccjs <- compileNimble(cjs)
cjsSpec <- MCMCspec(cjs, print = TRUE)
cjsMCMC <- buildMCMC(cjsSpec)
CcjsMCMC <- compileNimble(cjsMCMC, project = cjs) 

# run code for 10000 iterations
niter <- 10000
set.seed(0)
CcjsMCMC(niter)

# values from posterior distributions
samples <- as.matrix(nfVar(CcjsMCMC, 'mvSamples'))

# traceplots
par(mfrow = c(1, 2))
plot(samples[ , 'phi'], type = 'l', xlab = 'iteration',ylab = expression(phi))
plot(samples[ , 'p'], type = 'l', xlab = 'iteration',ylab = expression(p))

# mean values
mean(samples[, 'phi'])
mean(samples[, 'p'])


