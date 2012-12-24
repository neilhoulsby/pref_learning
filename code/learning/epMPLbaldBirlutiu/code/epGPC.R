#
# R script that implements the EP method for multitask preference learning.
#
# Author: Jose Miguel Hernandez Lobato
#
# Date: 28 Jan 2011
#

epGPCbirlutiu <- function(ratings, indexRatings, fullPriorMean, fullPriorCov) {

	# We use only the first two observations

	nRatings <- length(ratings)

	# We compute the local covariance and the local mean

	localPriorMean <- fullPriorMean[ indexRatings ]
	localPriorCov <- fullPriorCov[ indexRatings, indexRatings ]

	# We initialize the posterior approximation and the factor approximations

	# The first approximate factor

	f1Hat <- list(mu = rep(0, nRatings), vu = rep(Inf, nRatings))

	# The second approximate factor

	f2Hat <- list(mu = rep(0, nRatings), vu = rep(Inf, nRatings))

	# The posterior approximation with all the interesting information

	pa <- list(mu = rep(0, nRatings), vu = rep(Inf, nRatings),
		f1Hat = f1Hat, f2Hat = f2Hat, 
		ratings = ratings, indexRatings = indexRatings,
		fullPriorMean = fullPriorMean, fullPriorCov = fullPriorCov,
		fullPriorCovInv = chol2inv(chol(fullPriorCov)),
		localPriorMean = localPriorMean, localPriorCov = localPriorCov,
		localPriorCovInv = chol2inv(chol(localPriorCov)),
		indexRatings = indexRatings, Sigma = NULL)

	##
	## We start refinint the factors for the first iteration
	## 

	# We refine the second approximate factor

	pa$f2Hat$vu <- diag(pa$localPriorCov)

	##
	## Main loop of EP
	##

	damping <- 0.5
	convergence <- F
	iteration <- 1
	while ((!convergence && iteration <= 1000) || iteration < 2) {

		paOld <- pa

		##
		## We refine the first approximate factor
		##

		mOld <- pa$f2Hat$mu
		vOld <- pa$f2Hat$vu
		Y <- pa$ratings

		logZ <- pnorm(Y * mOld / sqrt(vOld + 1), log.p = T)
		ratio <- exp(-logZ + dnorm(mOld / sqrt(vOld + 1), log = T))
		alpha <- ratio * Y / sqrt(vOld + 1)
		beta <- -ratio * (Y * mOld / sqrt(vOld + 1) + ratio) / (1 + vOld)
                eta2HatNew <- -beta / (1 + beta * vOld)
                eta1HatNew <- (alpha - mOld * beta) / (1 + beta * vOld)

		vuHatNew <- eta2HatNew^-1
		muHatNew <- eta1HatNew / eta2HatNew

		index <- which(vuHatNew < 0)
		muHatNew[ index ] <- pa$f1Hat$mu[ index ]
		vuHatNew[ index ] <- pa$f1Hat$vu[ index ]

		vuHatNew[ is.infinite(vuHatNew) ] <- 1e300

			# We do damping

		pa$f1Hat$mu <- damping * muHatNew / vuHatNew + (1 - damping) * pa$f1Hat$mu / pa$f1Hat$vu
		pa$f1Hat$vu <- (damping * vuHatNew^-1 + (1 - damping) * pa$f1Hat$vu^-1)^-1
		pa$f1Hat$mu <- pa$f1Hat$vu * pa$f1Hat$mu

		##
		## We refine the second approximate factor
		##

		pa$Sigma <- chol2inv(chol(pa$localPriorCovInv + diag(as.double(pa$f1Hat$vu^-1), nRatings)))
		vuNew <- diag(pa$Sigma)
		muNew <- as.double(pa$Sigma %*% (pa$f1Hat$mu * pa$f1Hat$vu^-1 + pa$localPriorCovInv %*% pa$localPriorMean))

		# We update the fourth approximate factor
		
		vuHatNew <- 1 / (1 / vuNew - 1 / pa$f1Hat$vu)
		muHatNew <- vuHatNew * (muNew / vuNew - pa$f1Hat$mu / pa$f1Hat$vu)

		# We only update those terms with positive variances

		index <- which(vuHatNew < 0)
		vuHatNew[ index ] <- pa$f2Hat$vu[ index ]
		muHatNew[ index ] <- pa$f2Hat$mu[ index ]

		# We do damping

		pa$f2Hat$mu <- damping * muHatNew / vuHatNew + (1 - damping) * pa$f2Hat$mu / pa$f2Hat$vu
		pa$f2Hat$vu <- (damping * vuHatNew^-1 + (1 - damping) * pa$f2Hat$vu^-1)^-1
		pa$f2Hat$mu <- pa$f2Hat$vu * pa$f2Hat$mu

			# We update the posterior approximation

		pa$vu <- 1 / (1 / pa$f1Hat$vu + 1 / pa$f2Hat$vu)
		pa$mu <- pa$vu * (pa$f1Hat$mu / pa$f1Hat$vu + pa$f2Hat$mu / pa$f2Hat$vu)

		# We check for convergence

		convergence <- checkConvergence(pa, paOld, iteration)

		iteration <- iteration + 1
	}

	# We compute the full posterior mean and the full posterior covariance

	auxInvVar <- rep(0, nrow(pa$fullPriorCovInv))
	auxInvVar[ pa$indexRatings ] <- pa$f1Hat$vu^-1
	auxMean <- rep(0, length(auxInvVar))
	auxMean[ pa$indexRatings ] <- pa$f1Hat$mu
	pa$fullPosteriorCov <- chol2inv(chol(pa$fullPriorCovInv + diag(auxInvVar, nrow(pa$fullPriorCovInv))))
	pa$fullPosteriorMean <- pa$fullPosteriorCov %*% (pa$fullPriorCovInv %*% pa$fullPriorMean + auxInvVar * auxMean)

	# We are done

	pa
}

##
# Function that checks for convergence of the EP method
#

checkConvergence <- function(pa, paOld, iteration) {

	change <- max(abs(pa$mu - paOld$mu))

#	cat(iteration, change, "\n")
	
	if (change < 1e-3)
		T
	else 
		F
}

##
# Function that maps item ids to their corresponding indexes
#

mapRatingToIndex <- function(itemIdA, itemIdB, nTotalItems) {

	itemFeatures <- matrix(seq(1, nTotalItems), nTotalItems, 1)
	n <- nrow(itemFeatures)
	X <- cbind(apply(itemFeatures, 2, function(x) sapply(x, function(x)
		kronecker(matrix(1, n, 1), x))), kronecker(matrix(1, n, 1), itemFeatures))
	
	apply(cbind(itemIdA, itemIdB), 1, function(x) which(X[ ,1 ] == x[ 1 ] & X[ , 2 ] == x[ 2 ]))
}

##
# Function that generates a prediction for the preference of a user on several pairs of items
# in a vectorized form.
#

predictGPCbirlutiu <- function(pa, indexRatings) {

	list(m = pa$fullPosteriorMean[ indexRatings ], v = diag(pa$fullPosteriorCov)[ indexRatings ])
}

##
# Function that computes the preference learning kernel matrix.
#

computeKernel <- function(itemFeatures, lengthScale) {

	n <- nrow(itemFeatures)
	X <- cbind(apply(itemFeatures, 2, function(x) sapply(x, function(x)
		kronecker(matrix(1, n, 1), x))), kronecker(matrix(1, n, 1), itemFeatures))
	
	# We obtain the features of the different items

	m <- nrow(X)
	d <- ncol(X) / 2
	itemsA <- matrix(X[ , 1 : d ], m, d)
	itemsB <- matrix(X[ , (d + 1) : (2 * d) ], m, d)

	# We compute the kernel matrix using the preference kernel

	QA <- matrix(apply(itemsA^2, 1, sum), m, m)
	QB <- matrix(apply(itemsB^2, 1, sum), m, m)

	distanceA <- QA + t(QA) - 2 * itemsA %*% t(itemsA)
	KA <- exp(-0.5 * distanceA)

	distanceB <- QB + t(QB) - 2 * itemsB %*% t(itemsB)
	KB <- exp(-0.5 * distanceB)

	distanceAB <- QA + t(QB) - 2 * itemsA %*% t(itemsB)
	KAB <- exp(-0.5 * distanceAB)

	distanceBA <- QB + t(QA) - 2 * itemsB %*% t(itemsA)
	KBA <- exp(-0.5 * distanceBA)

	K <- KA + KB - KAB - KBA

	K + diag(rep(1e-5), m)
}

##
# Function that implements the EM algorithm for a particular training set
#

emBirlutiu <- function(itemFeatures, train) {

	# We initialize the mean and covariance for the priors

	fullPriorMean <- rep(0, nrow(itemFeatures)^2)
	fullPriorCov <- computeKernel(itemFeatures, 1)

	# We extract the userIds and the rating values

	userIdsTrain <- train[ , 3 ]
	ratingsTrain <- train[ , 4 ]

	# We start the EM algorithm

	fullPriorMeanOld <- Inf
	fullPriorCovOld <- Inf
	convergence <- F
	iteration <- 1
	while (!convergence && iteration <= 20) {
		
		# We fit the different Gaussian processes

		gpc <- list()
		for (k in unique(userIdsTrain)) {
			index <- which(userIdsTrain == k)
			indexRatings <- mapRatingToIndex(train[ index, 1 ], train[ index, 2 ], nrow(itemFeatures))
			gpc[[ k ]] <- epGPCbirlutiu(ratingsTrain[ index ], indexRatings, fullPriorMean, fullPriorCov)
		}

		# We update the prior mean 

		fullPriorMean <- 0
		for (k in unique(userIdsTrain))
			fullPriorMean <- fullPriorMean + gpc[[ k ]]$fullPosteriorMean
		fullPriorMean <- fullPriorMean / length(unique(userIdsTrain))

		# We update the prior covariance

		fullPriorCov <- 0
		for (k in unique(userIdsTrain)) {
			meanDiff <- gpc[[ k ]]$fullPosteriorMean - fullPriorMean
			fullPriorCov <- fullPriorCov + gpc[[ k ]]$fullPosteriorCov + meanDiff %*% t(meanDiff)
		}
		fullPriorCov <- fullPriorCov / length(unique(userIdsTrain))

		# We check for convergence

		change <- max(abs(fullPriorMean - fullPriorMeanOld))
		change <- max(change, max(abs(fullPriorCov - fullPriorCovOld)))
		if (change < 1e-2)
			convergence <- T
		print(change)

		fullPriorMeanOld <- fullPriorMean
		fullPriorCovOld <- fullPriorCov

		iteration <- iteration + 1
	}

	# We are done

	gpc
}
