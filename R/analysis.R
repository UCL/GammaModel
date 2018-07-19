#-----------------------------------------------------------------------------------------------------------
# Full analysis to replicate results in paper
#-----------------------------------------------------------------------------------------------------------
# Overheads
#-----------------------------------------------------------------------------------------------------------
source('functions.R')
data <- read.csv('sheep.goat.NISP.csv',row.names=1)
models.redding <- read.csv('models.redding.csv',row.names=1)
models.payne <- read.csv('models.payne.csv',row.names=1)
#-----------------------------------------------------------------------------------------------------------
# Goodness of fit: Calculating p-values under each model for table 4
#-----------------------------------------------------------------------------------------------------------
# specify one dataset (site) and one model. 

# For example, how well the 'FON1' data fits the payne meat model 
GOF(counts = data['FON1',] , model = models.payne['meat',])

# Or 'TRA1' fits the redding energy model
GOF(counts = data['TRA1',] , model = models.redding['energy',])

#-----------------------------------------------------------------------------------------------------------
# Maximum Likelihood (ML) Gamma parameters for Fig 1
#-----------------------------------------------------------------------------------------------------------
# specify one dataset, and if necessary specify the class ages 

# for example find the Maximum Likelihood Gamma parameters of the 'FON1' data, using the default Payne class ages, showing the DEoptimR trace
gamma.ML.parameters(counts = data['FON1',], trace=T)

# find the Maximum Likelihood Gamma Parameters of the 'TRA1' data, and plot the Gamma disttribution defined by these parameters
par <- gamma.ML.parameters(counts = data['TRA1',])
x <- seq(0,10,length.out=10000)
shape <- par[1]
rate <- par[1]/par[2]
y <- dgamma(x, shape, rate)
plot(x,y,type='l',xlab='age', ylab='probability density')
#-----------------------------------------------------------------------------------------------------------
# MCMC joint parameter distributions for Fig 2
#-----------------------------------------------------------------------------------------------------------
# specify one dataset, and if necessary specify the class ages 

# for example the full joint distribution of the Gamma parameters for the 'TRA1' data, using the default Payne class ages, plotting the MCMC chain to check behaviour
pars <- mcmc(counts = data['TRA1',], plot.chain=T)
#-----------------------------------------------------------------------------------------------------------
# AIC and BIC for table 5
#-----------------------------------------------------------------------------------------------------------
# first, find the log Maximum Likelihood Estimate (MLE) under the best fitting Gamma model, and the best fitting Age Class model

# for example using the 'TRA1' data
gamma.ll <- gamma.log.MLE(counts = data['TRA1',])
ageclass.ll <- multinomial.log.MLE(counts = data['TRA1',])

# next calculate the AIC under both models (Gamma has 2 parameters, Multinomial has 8)
gamma.AIC <- (2*2) - (2*gamma.ll)
ageclass.AIC <- (2*8) - (2*ageclass.ll)

# similarly calculate the BIC
gamma.BIC <- (log(47)*2) - (2*gamma.ll)
ageclass.BIC <- (log(47)*8) - (2*ageclass.ll)
#-----------------------------------------------------------------------------------------------------------
# Ensure all Plot functions are included
#-----------------------------------------------------------------------------------------------------------

