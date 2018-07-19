#-----------------------------------------------------------------------------------------------------------
# All functions required for analysis
#-----------------------------------------------------------------------------------------------------------
require(gtools)
require(dplyr)
require(combinat)
require(DEoptimR)
require(scales)
require(LaplacesDemon)
#-----------------------------------------------------------------------------------------------------------
check.counts.format <- function(counts){
	# a few basic checks to ensure the format of the raw counts data is in the correct format
	# Basic required format is a one row data.frame comprising column names of age classes (or multi-classes) in capital letters from 'A'
	# For example, Payne's 1973 caprine age classes are 'A','B','C','D','E','F','G','H',and 'I'
	# More likely, the raw counts will include some multi-class assignments, for example 'ABCD','A','AB','B','BC' etc. Order does not matter.
	# row name specifies the name of the data, e.g. the site/phase or any other ID
	
	# counts: a one-row data.frame comprising integer counts, with column names of age classes (or multi-classes) in capital letters from 'A'
	
	# check is a data.frame
	if(!is.data.frame(counts))stop('Raw count data must be a dataframe')

	# Check just one row
	if(nrow(counts)!=1)stop('Raw count data must comprise a single row dataframe')

	# check integers
	if(!is.integer(as.matrix(counts)))stop('Raw count data must comprise only integers')

	# check col names are unique
	if(sum(duplicated(names(counts)))!=0)stop('Column names of raw count data must be unique')

	# check col names are capital letters
	letters <- unique(strsplit(paste(names(counts),collapse=''),split='')[[1]])
	if(sum(!letters%in%LETTERS)!=0) stop("Column names of raw count data must comprise only capital letters, eg 'EFG' or 'E' ")
	}
#-----------------------------------------------------------------------------------------------------------
unique.class.maker <- function(x){
	# memory demands can be reduced hugely by retaining a multiclass that ALWAYS (every row of 'counts') have counts combined 

	# x: string of class names, comprising capital letters from A. Eg c('ABCD','A','AB','B',...)

	classes.individual <- sort(unique(strsplit(paste(x,collapse=''),split='')[[1]]))
	keep <- c()
	for(n in 1:length(classes.individual)){
		member <- grepl(classes.individual[n],x) * nchar(x)
		member[member==0] <- NA
		keep <- c(keep, which(member==min(member,na.rm=T)))
		}
return(x[unique(keep)])}
#-----------------------------------------------------------------------------------------------------------
all.arrangements <- function(counts){
	# generate a matrix of all possible arrangements of the observed counts
 	# Memory demands from the total number of arrangements rapidly becomes huge if there are many counts across many multiclasses
	# If data have large numbers, calculating likelihoods exactly therefore becomes inpractical using this approach.

	# counts: data.frame satisfying the requirements of check.counts.format()


	check.counts.format(counts)

	classes <- unique.class.maker(names(counts))

	m1 <- matrix(0,1,length(classes)); colnames(m1) <- classes

	for(n in 1:length(counts)){

		# overheads
		class.counts <- counts[n]

		# separate into the classes found in 'classes'
		member <- c()
		for(c in 1:length(classes))member[c] <- (grepl(classes[c],names(class.counts)))
		class.names <- classes[member]

		# all possible arrangements of class counts, for the n th (multi)-class
		x <- t(xsimplex(length(class.names),as.integer(class.counts)))
		colnames(x) <- class.names
		m2 <- matrix(0,nrow(x),ncol(m1)); colnames(m2) <- classes
		m2[,colnames(x)] <- x

		# save memory, store as integers
		storage.mode(m2) <- 'integer'

		# all possible combinations of m1 and m2
		multi.1 <- do.call("rbind", rep(list(m1), nrow(m2)))
		multi.2 <- m2[rep(1:nrow(m2),each=nrow(m1)),]
		both <- multi.1 + multi.2

		# keep unique rows, and iterate
		m1 <- distinct(as.data.frame(both))
		}
return(m1)}
#-----------------------------------------------------------------------------------------------------------
random.arrangement <- function(counts){
	# Helper function required by GOF()
	# Generates one possible random realisation (conversion of multi-class counts into single class counts), under the assumption that the 
	# archaeologist's original multiclass assumption represents an equal probability of belonging to each of the individual classes. 
	# I.e, a count belonging to multiclass 'AB' was believed to belong to either class 'A' or 'B' with an equal probability. 

	# counts: a one row data.frame satisfying the requirements of check.counts.format(). I.e, handles just one dataset at a time.
	
	# ensure only one row
	if(nrow(counts)>1)warning('random.arrangement() must be handed a one-row data.frame. Only the first row has been used.')
	counts <- counts[1,]

	# scrap any with a count of zero
	counts <- counts[,counts>0,drop=F]

	# separate into uncertain and certain calls, only uncertains need sampling
	cert.counts <- counts[nchar(names(counts))==1]
	uncert.counts <- counts[nchar(names(counts))>1]

	# blank df of counts
	letters <- LETTERS[1:9]
	all.counts <- as.data.frame(matrix(0,1,length(letters))); names(all.counts) <- letters

	# add summary of random assignment to each uncertain class
	for(n in 1:length(uncert.counts)){
		rand <- sample(strsplit(names(uncert.counts[n]),'')[[1]],as.numeric(uncert.counts[n]),replace=T)
		tb <- table(rand)
		all.counts[,names(tb)] <- all.counts[,names(tb)] + as.integer(tb)
		}

	# add the certain counts
	all.counts[,names(cert.counts)] <- all.counts[,names(cert.counts)] + as.integer(cert.counts)
return(all.counts)}
#-----------------------------------------------------------------------------------------------------------
combine.classes <- function(comb,model){
	# Helper function required by GOF()
	# Aggregates individual classes to match the model classes. 
	# Column names of 'comb' must be single capital letters that occur within the column names of 'model'

	# comb: data.frame of counts, in single classes only. Eg, output by random.arrangement()
	# model: data.frame of model age class probabilities, with headers

	class.agg <- names(model)
	comb.agg <- matrix(,nrow(comb),length(class.agg))
	for(n in 1:length(class.agg)){
		v <- strsplit(class.agg[n],'')[[1]]
		comb.agg[,n] <- rowSums(comb[,v,drop=F])
		}

	# keep only the unique combinations
	if(nrow(comb.agg)>1)comb.agg <- comb.agg[!duplicated(comb.agg),,drop=F]
	comb.agg <- as.data.frame(comb.agg); colnames(comb.agg) <- class.agg
return(comb.agg)}
#--------------------------------------------------------------------------------------------------------
GOF <- function(counts, model, N=10000){
	# Goodness of Fit
	# The general approach to calculating a p-value for n observed counts given the model probabilities is to calculate the probability of 
	# all possible arrangements of the n counts, then sum the probabilities that are smaller or equal to the probability of the observed arrangement.
	# This approach can be easily achieved if each count was assigned to only one class, but becomes computationally expensive when dealing with the additional problem of multi-class assignments.

	# Instead, the p-value can be estimated by randomly sampling arrangements of the observed data (this automatically takes care of their frequencies),
	# calculating a p-value for each using a chi-squared test, then calculating the mean p-value. This automatically becomes the probability of the observed data being as or more extreme
	# (in contrast to Fisher's method of combining p-values which answers a different question of how to combine p-values from independent trials).

	# counts: a one row data.frame satisfying the requirements of check.counts.format(). I.e, handles just one dataset at a time.
	# model: data.frame of model age class probabilities, with headers

	check.counts.format(counts)

	all.p <- numeric(N)
	for(n in 1:N){

		# random assignment of multi-class counts into single classes
		r <- random.arrangement(counts)

		# aggregate single classes to match the model classes
		r <- combine.classes(r,model)

		# chi-squared p-value
		x <- as.integer(r)
		p <- as.numeric(model)
		
		# remove model zero probability classes (meaningless for Chisq)
		i <- abs(p)>0.00000001
		x <- x[i]
		p <- p[i]
		
		all.p[n] <- suppressWarnings(chisq.test(x=x,p=p)$p.value)

		if(n%in%seq(0,N,by=1000))print(paste(n,'of',N,'samples completed'))
		}
	res <- formatC(mean(all.p), format = "e", digits = 1)

return(res)}
#-------------------------------------------------------------------------------------------------------
convert.class.ages <- function(class.ages,multi.classes){
	# helper function to assist in computational efficiency, since class.ages can be combined where observed data permit
	# Checks various requirements first.
	
	# ensure multi.classes only comprise sequential classes
	for(n in 1:length(multi.classes)){
		char <- (strsplit(multi.classes[n],split='')[[1]])
		i <- which(LETTERS%in%char)
		if(length(i)>1){
			problem <- sum(diff(i)!=1)
			if(problem!=0)stop(paste('Multiclass',multi.classes[n],'is not valid, must comprise sequential letters'))
			}
		}

	# ensure multi.classes only comprises capital letters
	mc <- strsplit(paste(multi.classes,collapse=''),split='')[[1]]
	if(sum(!mc%in%LETTERS)!=0)stop('multi.classes must comprise only capital letters')

	# ensure class ages only comprises capital letters
	ca <- strsplit(paste(names(class.ages),collapse=''),split='')[[1]]
	if(sum(!ca%in%LETTERS)!=0)stop('names of class.ages must comprise only capital letters')

	# ensure class.ages and multi.class ages comprise the same letters
	match <- identical(sort(unique(mc)),  sort(unique(ca)))
	if(!match)stop('names of class.ages and multi.classes dont match')

	# convert the class.ages to match the multi.classes
	N <- length(multi.classes)
	new <- as.data.frame(matrix(,1,N)); names(new) <- multi.classes
	for(n in 1:N){
		i <- strsplit(multi.classes[n],split='')[[1]][1] == names(class.ages)
		new[n] <- class.ages[i]
		}
return(new)}
#-------------------------------------------------------------------------------------------------------
gamma.loglik <- function(x,class.ages,shape,mean){
	# log likelihood under any Gamma distribution 
	# converts the continuous PDF of any Gamma distribution into discrete probabilities (a multinomial probability distribution) for each age class
	# column names of 'x' and 'class.ages' must match, ie the classes (or multi-classes) must be the same

	# x: data.frame of all possible counts in each age class for which likelihoods are to be calculated and summed. Ie, the output of all.arrangements()
	# class.ages: a one row data.frame of the starting age of each age class or multiclass
	# These must be sequential, such that the start of the (n+1)th class equals the end of the nth class.
	# The final age class will automatically include all counts above this. The first starting age therefore must be zero. 
	# shape, mean: gamma shape parameters

	if(!identical(names(x),names(class.ages)))stop("The class names of 'x' and 'class.names' must match")
	if(class.ages[1]!=0)stop('The first value must be zero to ensure all possible ages are encompassed')
	rate <- shape/mean
	p <- diff(c(pgamma(q=as.numeric(class.ages),shape=shape,rate=rate),1))
	loglik <- log(sum(apply(x, 1, dmultinom, prob=p, log=FALSE)))

	# in the extreme, the probabilities can shrink to almost zero (loglik becomes -Inf). This causes trouble downstream in MCMC acceptance ratios
	if(is.infinite(loglik))loglik <- -99999
return(loglik)}
#-----------------------------------------------------------------------------------------------------------
multinomial.loglik <- function(x,p){
	# calculates the log likelihood under any multinomial distribution
	# column names of 'x' and 'p' must match, ie the classes (or multi-classes) must be the same

	# x: data.frame of all possible arrnagements of counts in each age class for which likelihoods are to be calculated and summed. Ie, the output of all.arrangements()
	# p: a one row data.frame of model probabilites in the (multi-class) age classes. 

	if(!identical(names(x),names(p)))warning("The class names of 'x' and 'p' must match")
	loglik <- log(sum(apply(x, 1, dmultinom, prob=as.numeric(p), log=FALSE)))
return(loglik)}
#-----------------------------------------------------------------------------------------------------------
proposal.function <- function(params,size=0.4){
	# helper function for MCMC
	# params: vector or two values representing the 'shape' and 'mean' parameters of the gamma distribution
	# size: hyperparameter controlling the size of the next random jump. May require tuning: smaller jumps result in a higher acceptance ratio (AR). Aim for AR between 0.2 and 0.6
	
	N <- length(params)
	jumps <- rep(size,N) 
	moves <- rnorm(N,0,jumps)
	new.params <- abs(params + moves)
return(new.params)}
#-----------------------------------------------------------------------------------------------------------
mcmc <- function(counts, class.ages=data.frame(A=0,B=1/6,C=1/2,D=1,E=2,F=3,G=4,H=6,I=8), N=30000, burn=2000, thin=5, plot.chain=T){
	# Markov Chain Monte Carlo. Generates a single chain of Gamma parameters (shape and mean) for some observed multi-class assignments
	# Returns the matrix of Gamma parameters after thinning and removing burn-in

	# counts: a one row data.frame satisfying the requirements of check.counts.format(). I.e, handles just one dataset at a time.
	# class.ages: a one row data.frame of the starting age of each age class. See gamma.loglik()
	# N: number of samples in the chain
	# burn: number of initial burn in samples to discard
	# thin: proportion of samples to discard. I.e. 5 = every 5th sample in the chain is retained.

	check.counts.format(counts)
	aa <- all.arrangements(counts)
	ca <- convert.class.ages(class.ages,names(aa))

	# initiate random Gamma parameters, between 0 and 5
	params <- runif(2,0,5)
	all.params <- matrix(,N,2)

	# mcmc loop
	accepted <- rep(0,N)
	for(n in 1:N){
		all.params[n,] <- params
		llik <- gamma.loglik(aa,ca,params[1],params[2])
		prop.params <- proposal.function(params)
		prop.llik <- gamma.loglik(aa,ca,prop.params[1],prop.params[2])
		ratio <- min(exp(prop.llik-llik),1)
		move <- sample(c(T,F),size=1,prob=c(ratio,1-ratio))
		if(move){
			params <- prop.params
			accepted[n] <- 1
			}
		if(n%in%seq(0,N,by=1000))print(paste(n,'of',N,'samples completed'))
		}

	# remove burnin
	N.sub <- all.params[burn:N,]

	# thinning
	i <- seq(burn,N,by=thin)
	res <- all.params[i,]

	# look at chain
	if(plot.chain){
		par(mfrow=c(3,2))
		plot(all.params[,1],type='l',xlab='Samples in chain',ylab='Gamma shape parameter',main='Entire chain')
		plot(all.params[,2],type='l',xlab='Samples in chain',ylab='Gamma mean parameter',main='Entire chain')
		hist(res[,1],breaks=100,xlab='Gamma shape parameter',main='Burn in removed and thinned',ylab='count')
		hist(res[,2],breaks=100,xlab='Gamma mean parameter',main='Burn in removed and thinned',ylab='count')
		plot(res,pch='.',xlab='Gamma shape parameter',ylab='Gamma mean parameter',main='Burn in removed and thinned')		
		}

	# print acceptance ratio
	print(paste('acceptance ratio = ',round(sum(accepted[burn:N])/nrow(N.sub),3)))

return(res)}
#-----------------------------------------------------------------------------------------------------------
gamma.ML.parameters <- function(counts,class.ages=data.frame(A=0,B=1/6,C=1/2,D=1,E=2,F=3,G=4,H=6,I=8),trace=F){
	# Search for the Maximum Likelihood Gamma parameters that fit any multi-class counts, using DEoptimR

	# counts: a one row data.frame satisfying the requirements of check.counts.format(). I.e, handles just one dataset at a time.
	# class.ages: a one row data.frame of the starting age of each age class. See gamma.loglik()

	check.counts.format(counts)
	aa <- all.arrangements(counts)
	class.ages <- convert.class.ages(class.ages,names(aa))

	fn <- function(x,aa,class.ages){-gamma.loglik(aa,class.ages,x[1],x[2])}
	res <- JDEoptim(lower=c(0,0),upper=c(20,20),fn=fn,aa=aa,class.ages=class.ages,tol=1e-7,NP=10,trace=trace)
return(round(res$par,3))}
#-----------------------------------------------------------------------------------------------------------
gamma.log.MLE <- function(counts,class.ages=data.frame(A=0,B=1/6,C=1/2,D=1,E=2,F=3,G=4,H=6,I=8),trace=F){
	# Search for the log Maximum Likelihood Estimate under the best fitting Gamma
	# Uses exactly the same algorithm as gamma.ML.parameters()

	# counts: a one row data.frame satisfying the requirements of check.counts.format(). I.e, handles just one dataset at a time.
	# class.ages: a one row data.frame of the starting age of each age class. See gamma.loglik()

	check.counts.format(counts)
	aa <- all.arrangements(counts)
	class.ages <- convert.class.ages(class.ages,names(aa))

	fn <- function(x,aa,class.ages){-gamma.loglik(aa,class.ages,x[1],x[2])}
	res <- JDEoptim(lower=c(0,0),upper=c(20,20),fn=fn,aa=aa,class.ages=class.ages,tol=1e-7,NP=10,trace=trace)
return(-res$value)}
#-----------------------------------------------------------------------------------------------------------
multinomial.ML.parameters <- function(counts,N=100,I=30,C=2){
	# Search for the Maximum Likelihood multinomial parameters (probabilities)
	# Search performed using Simulated Annealing to sample the Dirichlet simplex, since a required constraint is that the parameters sum to 1
	# Although column sums of 'aa' (normalised for unity) gives a very rough approximation of the MLE, it is incorrect since 'aa' provides all possible arrangements, but not the frequencies of each arrangement.
	# Instead, this function properly sums the likelihoods across all possible arrangements.

	# counts: a one row data.frame satisfying the requirements of check.counts.format(). I.e, handles just one dataset at a time.
	# N: pop size
	# I: iterations
	# C: cooling: annealing schedule

	check.counts.format(counts)
	aa <- all.arrangements(counts)
	R <- ncol(aa)

	# first iteration
	liks <- numeric(N)
	pars <- rdirichlet(N,rep(1,R))
	for(n in 1:N)liks[n] <- log(sum(apply(aa, 1, dmultinom, prob=as.numeric(pars[n,]), log=FALSE)))
	keep.pars <- pars[liks==max(liks),,drop=F]

	# subsequent iterations
	for(i in 1:I){
		liks <- numeric(N+1)
		pars <- rbind(keep.pars,rdirichlet(N,rep(0.1,R)+keep.pars*(C^i)))
		for(n in 1:(N+1))liks[n] <- log(sum(apply(aa, 1, dmultinom, prob=as.numeric(pars[n,]), log=FALSE)))
		keep.pars <- pars[liks==max(liks),,drop=F]	
		}

	res <- as.data.frame(matrix(round(keep.pars,4),1,R));names(res) <- names(counts) 
return(res)}
#-----------------------------------------------------------------------------------------------------------
multinomial.log.MLE <- function(counts,N=100,I=30,C=2){
	# Search for the log Maximum Likelihood Estimate under the best fitting Multinomial
	# Uses exactly the same algorithm as multinomial.ML.parameters()

	# counts: a one row data.frame satisfying the requirements of check.counts.format(). I.e, handles just one dataset at a time.
	# N: pop size
	# I: iterations
	# C: cooling: annealing schedule

	check.counts.format(counts)
	aa <- all.arrangements(counts)
	R <- ncol(aa)

	# first iteration
	liks <- numeric(N)
	pars <- rdirichlet(N,rep(1,R))
	for(n in 1:N)liks[n] <- log(sum(apply(aa, 1, dmultinom, prob=as.numeric(pars[n,]), log=FALSE)))
	keep.pars <- pars[liks==max(liks),,drop=F]

	# subsequent iterations
	for(i in 1:I){
		liks <- numeric(N+1)
		pars <- rbind(keep.pars,rdirichlet(N,rep(0.1,R)+keep.pars*(C^i)))
		for(n in 1:(N+1))liks[n] <- log(sum(apply(aa, 1, dmultinom, prob=as.numeric(pars[n,]), log=FALSE)))
		keep.pars <- pars[liks==max(liks),,drop=F]	
		}

	res <- max(liks)
return(res)}
#-----------------------------------------------------------------------------------------------------------

