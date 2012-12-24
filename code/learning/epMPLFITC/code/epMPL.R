#
# R script that implements the EP method for multitask preference learning.
#
# Author: Jose Miguel Hernandez Lobato
#
# Date: 23 Dec 2011
#

source("vbpcad.R")

source("epFITC.R")

##
# Function that implements the EP-VB algorithm for the multitask model
#
# @param	problemInfo	A list with the following fields:
#				
#				itemFeaturesA		-> n x q matrix with the features of the first rated item.
#				itemFeaturesB		-> n x q matrix with the features of the second rated item.
#				itemIdsA		-> n-dimensional vector with the id of the first item compared by the users.
#				itemIdsB		-> n-dimensional vector with the id of the first item compared by the users.
#				userIds			-> n-dimensional vector with the ids of the users.
#				ratings			-> n-dimensional vector with the ratings given by the users.
#				d			-> Number of latent functions.
#				lengthScaleItems	-> Length-scale of the cov. function for the shared latent functions.
#				lengthScaleUsers	-> Length-scale of the cov. function for user weights.
#				nTotalItems		-> Number of total items.
#				userFeatures		-> U x p matrix with the user features.
#

epMPLinternal <- function(problemInfo) {

	# We add the additional structure for the FITC approximation

	problemInfo <- addFITCauxiliarVariables(problemInfo)

	# Total number of different pairs of items rated and the number of users

	nRatings <- length(unique(problemInfo$ratingIds))
	nUsers <- length(unique(problemInfo$userIds))
	d <- problemInfo$d

	# We initialize the posterior approximation and the factor approximations

	# The first approximate factor

	f1Hat <- list(mu = matrix(0, nUsers, nRatings), vu = matrix(Inf, nUsers, nRatings))

	# The second approximate factor

	f2Hat <- list(mu = matrix(0, nUsers, nRatings), vu = matrix(Inf, nUsers, nRatings),
		mw = matrix(0, nUsers, d), vw = matrix(Inf, nUsers, d),
		mh = matrix(0, d, nRatings), vh = matrix(Inf, d, nRatings))

	# The third approximate factor

	f3Hat <- list(mw = matrix(0, nUsers, d), vw = matrix(Inf, nUsers, d))

	# The fourth approximate factor

	f4Hat <- list(mh = matrix(0, d, nRatings), vh = matrix(Inf, d, nRatings))

	# The posterior approximation

	a <- list(mu = matrix(0, nUsers, nRatings), vu = matrix(Inf, nUsers, nRatings),
		mw = matrix(0, nUsers, d), vw = matrix(1, nUsers, d),
		mh = matrix(0, d, nRatings), vh = matrix(Inf, d, nRatings))

	##
	## We start refinint the factors for the first iteration
	## 

	# We refine the fourth approximate factor

	a$vh <- f4Hat$vh <- matrix(problemInfo$fitcItems$diagKn, d, nRatings, byrow = T)

	# We refine the third approximate factor

	a$vw <- f3Hat$vw <- matrix(problemInfo$fitcUsers$diagKn, nUsers, d)

	# We refine the second approximate factor

	a$vu <- f2Hat$vu <- f3Hat$vw %*% f4Hat$vh

	# The variational solution is initialized to NULL

	variationalSolution <- NULL

	##
	## Main loop of EP
	##

	# We check for an initial solution

	damping <- 0.5
	convergence <- F
	iteration <- 1
	while ((!convergence && iteration <= 1000) || iteration <= 3) {

		aOld <- a

		##
		## We refine the first approximate factor
		##

		mOld <- f2Hat$mu
		vOld <- f2Hat$vu
		Y <- problemInfo$ratingMatrix

		logZ <- pnorm(Y * mOld / sqrt(vOld + 1), log.p = T)
		ratio <- exp(-logZ + dnorm(mOld / sqrt(vOld + 1), log = T))
		alpha <- ratio * Y / sqrt(vOld + 1)
		beta <- -ratio * (Y * mOld / sqrt(vOld + 1) + ratio) / (1 + vOld)
                eta2HatNew <- -beta / (1 + beta * vOld)
                eta1HatNew <- (alpha - mOld * beta) / (1 + beta * vOld)

		vuHatNew <- eta2HatNew^-1
		muHatNew <- eta1HatNew / eta2HatNew

		muHatNew[ is.infinite(vuHatNew) ] <- 0
		vuHatNew[ is.infinite(vuHatNew) ] <- 1e300

		index <- which(vuHatNew < 0)
		muHatNew[ index ] <- f1Hat$mu[ index ]
		vuHatNew[ index ] <- f1Hat$vu[ index ]

			# We do damping

		f1Hat$mu <- damping * muHatNew / vuHatNew + (1 - damping) * f1Hat$mu / f1Hat$vu
		f1Hat$vu <- (damping * vuHatNew^-1 + (1 - damping) * f1Hat$vu^-1)^-1
		f1Hat$mu <- f1Hat$vu * f1Hat$mu

		##
	      	## We refine the second approximate factor
		##

			# We create the rating entries for the variational method

		ratingEntries <- cbind(problemInfo$userIds, problemInfo$ratingIds,
			f1Hat$mu[ cbind(problemInfo$userIds, problemInfo$ratingIds) ],
			f1Hat$vu[ cbind(problemInfo$userIds, problemInfo$ratingIds) ])
		
			# We call the optimization method using the previous solution if iteration > 1

		if (is.null(variationalSolution))
			ret <- bvpcadFast(ratingEntries, d, f3Hat$mw, f3Hat$vw, t(f4Hat$mh), t(f4Hat$vh), 0)
		else
			ret <- bvpcadFast(ratingEntries, d, f3Hat$mw, f3Hat$vw, t(f4Hat$mh), t(f4Hat$vh), 0, variationalSolution)

		variationalSolution <- ret
		vbLowerBound <- ret$bound

			# We update the second approximate factor

		muNew <- ret$mP %*% t(ret$mQ)
		vuNew <- ret$mP^2 %*% t(ret$vQ) + ret$vP %*% t(ret$mQ^2) + ret$vP %*% t(ret$vQ)

		vuHatNew <- 1 / (1 / vuNew - 1 / f1Hat$vu)
		muHatNew <- vuHatNew * (muNew / vuNew - f1Hat$mu / f1Hat$vu)
		vwHatNew <- 1 / (1 / ret$vP - 1 / f3Hat$vw)
		mwHatNew <- vwHatNew * (ret$mP / ret$vP - f3Hat$mw / f3Hat$vw)
		vhHatNew <- 1 / (1 / t(ret$vQ) - 1 / f4Hat$vh)
		mhHatNew <- vhHatNew * (t(ret$mQ / ret$vQ) - f4Hat$mh / f4Hat$vh)

			# We update only those terms with positive variances

		index <- which(vuHatNew < 0)
		vuHatNew[ index ] <- f2Hat$vu[ index ]
		muHatNew[ index ] <- f2Hat$mu[ index ]
		index <- which(vwHatNew < 0)
		vwHatNew[ index ] <- f2Hat$vw[ index ]
		mwHatNew[ index ] <- f2Hat$mw[ index ]
		index <- which(vhHatNew < 0)
		vhHatNew[ index ] <- f2Hat$vh[ index ]
		mhHatNew[ index ] <- f2Hat$mh[ index ]

			# We do damping

		f2Hat$mu <- damping * muHatNew / vuHatNew + (1 - damping) * f2Hat$mu / f2Hat$vu
		f2Hat$vu <- (damping * vuHatNew^-1 + (1 - damping) * f2Hat$vu^-1)^-1
		f2Hat$mu <- f2Hat$vu * f2Hat$mu

		f2Hat$mw <- damping * mwHatNew / vwHatNew + (1 - damping) * f2Hat$mw / f2Hat$vw
		f2Hat$vw <- (damping * vwHatNew^-1 + (1 - damping) * f2Hat$vw^-1)^-1
		f2Hat$mw <- f2Hat$vw * f2Hat$mw

		f2Hat$mh <- damping * mhHatNew / vhHatNew + (1 - damping) * f2Hat$mh / f2Hat$vh
		f2Hat$vh <- (damping * vhHatNew^-1 + (1 - damping) * f2Hat$vh^-1)^-1
		f2Hat$mh <- f2Hat$vh * f2Hat$mh

		##
		## We refine the second approximate factor
		##

		for (i in 1 : d) {

			# We refine the approximate factor for the Gaussian process using the FITC approximation

			ret <- computeTitledDistribution(problemInfo$fitcItems$D, problemInfo$fitcItems$P, problemInfo$fitcItems$R,
				problemInfo$fitcItems$PRt, f2Hat$mh[ i, ], f2Hat$vh[ i, ])
			mhNew <- ret$mNew
			vhNew <- ret$vNew

			# We update the fourth approximate factor
		
			vhHatNew <- 1 / (1 / vhNew - 1 / f2Hat$vh[ i, ])
			mhHatNew <- vhHatNew * (mhNew / vhNew - f2Hat$mh[ i, ] / f2Hat$vh[ i, ])
			
			# We only update those terms with positive variances

			index <- which(vhHatNew < 0)
			vhHatNew[ index ] <- f4Hat$vh[ i, index ]
			mhHatNew[ index ] <- f4Hat$mh[ i, index ]

			# We do damping

			f4Hat$mh[ i, ] <- damping * mhHatNew / vhHatNew + (1 - damping) * f4Hat$mh[ i, ] / f4Hat$vh[ i, ]
			f4Hat$vh[ i, ] <- (damping * vhHatNew^-1 + (1 - damping) * f4Hat$vh[ i, ]^-1)^-1
			f4Hat$mh[ i, ] <- f4Hat$vh[ i, ] * f4Hat$mh[ i, ]
		}

		##
		## We refine the third approximate factor
		##

		for (i in 1 : d) {

			# We refine the approximate factor for the Gaussian process using the FITC approximation

			ret <- computeTitledDistribution(problemInfo$fitcUsers$D, problemInfo$fitcUsers$P, problemInfo$fitcUsers$R,
				problemInfo$fitcUsers$PRt, f2Hat$mw[ , i ], f2Hat$vw[ , i ])
			mwNew <- ret$mNew
			vwNew <- ret$vNew

			# We update the fourth approximate factor
		
			vwHatNew <- 1 / (1 / vwNew - 1 / f2Hat$vw[ , i ])
			mwHatNew <- vwHatNew * (mwNew / vwNew - f2Hat$mw[ , i ] / f2Hat$vw[ , i ])
			
			# We only update those terms with positive variances

			index <- which(vwHatNew < 0)
			vwHatNew[ index ] <- f3Hat$vw[ index, i ]
			mwHatNew[ index ] <- f3Hat$mw[ index, i ]

			# We do damping

			f3Hat$mw[ , i ] <- damping * mwHatNew / vwHatNew + (1 - damping) * f3Hat$mw[ , i ] / f3Hat$vw[ , i ]
			f3Hat$vw[ , i ] <- (damping * vwHatNew^-1 + (1 - damping) * f3Hat$vw[ , i ]^-1)^-1
			f3Hat$mw[ , i ] <- f3Hat$vw[ , i ] * f3Hat$mw[ , i ]
		}

		# We update the posterior approximation

		a$vu <- 1 / (1 / f1Hat$vu + 1 / f2Hat$vu)
		a$mu <- a$vu * (f1Hat$mu / f1Hat$vu + f2Hat$mu / f2Hat$vu)

		a$vw <- 1 / (1 / f3Hat$vw + 1 / f2Hat$vw)
		a$mw <- a$vw * (f3Hat$mw / f3Hat$vw + f2Hat$mw / f2Hat$vw)

		a$vh <- 1 / (1 / f4Hat$vh + 1 / f2Hat$vh)
		a$mh <- a$vh * (f4Hat$mh / f4Hat$vh + f2Hat$mh / f2Hat$vh)

		# We update the damping parameter

		damping <- damping * 0.95

		# We check for convergence

		change <- max(abs(a$mu[ cbind(problemInfo$userIds, problemInfo$ratingIds) ] -
			aOld$mu[ cbind(problemInfo$userIds, problemInfo$ratingIds) ]))
		if (change < 1e-2)
			convergence <- T
		else 
			convergence <- F

		cat(iteration, change, "\n")

		iteration <- iteration + 1
	}

	# We estimate the evidence

	evidence <- computeEvidence(a, problemInfo, f1Hat, f2Hat, f3Hat, f4Hat, vbLowerBound, nRatings, nUsers)

	# We return the posterior approximation

	list(a = a, evidence = evidence, problemInfo = problemInfo,
		f1Hat = f1Hat, f2Hat = f2Hat, f3Hat = f3Hat, f4Hat = f4Hat,
		vbLowerBound = vbLowerBound, variationalSolution = variationalSolution)
}

##
# Function that evaluates the EP approximation of the evidence
#

computeEvidence <- function(a, problemInfo, f1Hat, f2Hat, f3Hat, f4Hat, vbLowerBound, nRatings, nUsers) {

	# We compute the log s for the first factor

	logZ <- pnorm(problemInfo$ratingMatrix * f2Hat$mu / sqrt(f2Hat$vu + 1), log.p = T)
	
	logShat1 <- logZ + 0.5 * log(2 * pi) + 0.5 * log(f1Hat$vu) + 0.5 * log(f2Hat$vu) - 0.5 * log(a$vu) -
		0.5 * a$mu^2 / a$vu + 0.5 * f1Hat$mu^2 / f1Hat$vu + 0.5 * f2Hat$mu^2 / f2Hat$vu
	logShat1 <- sum(logShat1[ cbind(problemInfo$userIds, problemInfo$ratingIds) ])

	# We compute the log s for the second factor

	logZ <- vbLowerBound

	auxSum <- 0.5 * log(2 * pi) + 0.5 * log(f1Hat$vu) + 0.5 * log(f2Hat$vu) - 0.5 * log(a$vu) -
		0.5 * a$mu^2 / a$vu + 0.5 * f1Hat$mu^2 / f1Hat$vu + 0.5 * f2Hat$mu^2 / f2Hat$vu
	auxSum <- sum(auxSum[ cbind(problemInfo$userIds, problemInfo$ratingIds) ])

	logShat2 <- logZ + auxSum + sum(0.5 * log(2 * pi) + 0.5 * log(f4Hat$vh) + 0.5 * log(f2Hat$vh) - 0.5 * log(a$vh) -
		0.5 * a$mh^2 / a$vh + 0.5 * f4Hat$mh^2 / f4Hat$vh + 0.5 * f2Hat$mh^2 / f2Hat$vh) +
		sum(0.5 * log(2 * pi) + 0.5 * log(f3Hat$vw) + 0.5 * log(f2Hat$vw) - 0.5 * log(a$vw) -
		0.5 * a$mw^2 / a$vw + 0.5 * f3Hat$mw^2 / f3Hat$vw + 0.5 * f2Hat$mw^2 / f2Hat$vw)

	# We compute the log shat for the third factor

	logZ <- 0
	for (i in 1 : problemInfo$d)
		logZ <- logZ + getFITCevidence(problemInfo$fitcUsers$D, problemInfo$fitcUsers$P, problemInfo$fitcUsers$R,
				problemInfo$fitcUsers$PRt, f2Hat$mw[ , i ], f2Hat$vw[ , i ])

	logShat3 <- logZ + sum(0.5 * log(2 * pi) + 0.5 * log(f3Hat$vw) + 0.5 * log(f2Hat$vw) - 0.5 * log(a$vw) -
		0.5 * a$mw^2 / a$vw + 0.5 * f3Hat$mw^2 / f3Hat$vw + 0.5 * f2Hat$mw^2 / f2Hat$vw)

	# We compute the log shat for the fourth factor

		# First we evaluate the normalization constant Z. See Rasmussen's book

	logZ <- 0
	for (i in 1 : problemInfo$d)
		logZ <- logZ + getFITCevidence(problemInfo$fitcItems$D, problemInfo$fitcItems$P, problemInfo$fitcItems$R,
				problemInfo$fitcItems$PRt, f2Hat$mh[ i, ], f2Hat$vh[ i, ])

	logShat4 <- logZ + sum(0.5 * log(2 * pi) + 0.5 * log(f4Hat$vh) + 0.5 * log(f2Hat$vh) - 0.5 * log(a$vh) -
		0.5 * a$mh^2 / a$vh + 0.5 * f4Hat$mh^2 / f4Hat$vh + 0.5 * f2Hat$mh^2 / f2Hat$vh)

	# We evaluate the EP approximation of the evidence

	auxSum <- -0.5 * log(2 * pi) - 0.5 * log(f1Hat$vu) - 0.5 * log(f2Hat$vu) + 0.5 * log(a$vu) +
		0.5 * a$mu^2 / a$vu - 0.5 * f1Hat$mu^2 / f1Hat$vu - 0.5 * f2Hat$mu^2 / f2Hat$vu
	auxSum <- sum(auxSum[ cbind(problemInfo$userIds, problemInfo$ratingIds) ])

	logEvidence <- logShat1 + logShat2 + logShat3 + logShat4 + 
		auxSum + sum(-0.5 * log(2 * pi) - 0.5 * log(f4Hat$vh) - 0.5 * log(f2Hat$vh) + 0.5 * log(a$vh) +
		0.5 * a$mh^2 / a$vh - 0.5 * f4Hat$mh^2 / f4Hat$vh - 0.5 * f2Hat$mh^2 / f2Hat$vh) +
		sum(-0.5 * log(2 * pi) - 0.5 * log(f3Hat$vw) - 0.5 * log(f2Hat$vw) + 0.5 * log(a$vw) +
		0.5 * a$mw^2 / a$vw - 0.5 * f3Hat$mw^2 / f3Hat$vw - 0.5 * f2Hat$mw^2 / f2Hat$vw)

	logEvidence
}

##
# Function that generates a prediction for the preference of a user on several pairs of items
# in a vectorized form.
#

predictMPL <- function(ret, idUser, itemFeaturesA, itemFeaturesB) {

	# Matrices where to store the mean and variances of the conditional
	# Gaussian process distribution.

	m <- matrix(0, ret$problemInfo$d, nrow(itemFeaturesA))
	v <- matrix(0, ret$problemInfo$d, nrow(itemFeaturesA))

	# We compute the FITC prediction

	Xtest <- cbind(itemFeaturesA, itemFeaturesB)
	Xtest <- Xtest / ret$problemInfo$lengthScaleItems
	pStar <- t(spgpComputeKmnPreference(Xtest, ret$problemInfo$fitcItems$Xbar))
	dStar <- spgpComputeDiagKnPreference(Xtest) - rep(1, ret$problemInfo$fitcItems$m) %*% (ret$problemInfo$fitcItems$R %*% t(pStar))^2

	for (i in 1 : ret$problemInfo$d) {
		aux <- predictFITC(ret$problemInfo$fitcItems$D, ret$problemInfo$fitcItems$P, ret$problemInfo$fitcItems$R,
			ret$problemInfo$fitcItems$PRt, ret$f2Hat$mh[ i, ], ret$f2Hat$vh[ i, ], dStar, pStar)
		m[ i, ] <- aux$m
		v[ i , ] <- aux$v
	}

	mFinal <- apply(ret$a$mw[ idUser, ] * t(m), 1, sum)
	vFinal <- apply(ret$a$mw[ idUser, ]^2 * t(v), 1, sum) + apply(ret$a$vw[ idUser, ] * t(m)^2, 1, sum) + apply(ret$a$vw[ idUser, ] * t(v), 1, sum)

	list(m = mFinal, v = vFinal)
}

preprocessProblemInfo <- function(problemInfo) {
	
	# We extract each field from the list, selecting the item features of the rated items

	itemFeaturesA <- problemInfo$itemFeaturesA
	itemFeaturesB <- problemInfo$itemFeaturesB
	userIds <- problemInfo$userIds
	itemIdsA <- problemInfo$itemIdsA
	itemIdsB <- problemInfo$itemIdsB
	ratings <- problemInfo$ratings
	nTotalItems <- problemInfo$nTotalItems

	# The id of item A should always be smaller than the id of item B

	index <- which(itemIdsA > itemIdsB)
	if (length(index) != 0) {
		swap <- itemIdsA[ index, ]
		itemIdsA[ index, ] <- itemIdsB[ index, ]
		itemIdsB[ index, ] <- swap

		swap <- itemFeaturesA[ index, ]
		itemFeaturesA[ index, ] <- itemFeaturesB[ index, ]
		itemFeaturesB[ index, ] <- swap

		ratings[ index ] <- -ratings[ index ]
	}
	
	# We obtain the id of each combination of two items

	allItemIds <- matrix(seq(1, nTotalItems), nTotalItems, 1)
	n <- nrow(allItemIds)
	X <- cbind(apply(allItemIds, 2, function(x) sapply(x, function(x) kronecker(matrix(1, n, 1), x))), kronecker(matrix(1, n, 1), allItemIds))
	X <- X[ apply(X, 1, function(x) x[ 1 ] < x[ 2 ]), ]

	ratingIds <- apply(cbind(itemIdsA, itemIdsB), 1, function(x) which(X[ ,1 ] == x[ 1 ] & X[ , 2 ] == x[ 2 ]))

	# We compute the kernel matrix for the shared latent functions and its inverse.

	ratingIdsKernelMatrix <- sort(unique(ratingIds))
	index <- tapply(ratingIdsKernelMatrix, ratingIdsKernelMatrix, function(x) min(which(ratingIds == x)))
	itemFeaturesA <- itemFeaturesA[ index, ]
	itemFeaturesB <- itemFeaturesB[ index, ]

	# We map each rating id to its corresponding entry in the kernel matrix

	map <- rep(0, max(ratingIdsKernelMatrix))
	map[ ratingIdsKernelMatrix ] <- seq(1, length(ratingIdsKernelMatrix))
	ratingIds <- map[ ratingIds ]

	# We compute the rating matrix

	ratingMatrix <- matrix(0, length(unique(userIds)), length(unique(ratingIds)))
	ratingMatrix[ cbind(userIds, ratingIds) ] <- ratings

	# We return all the data

	list(itemFeaturesA = itemFeaturesA, itemFeaturesB = itemFeaturesB, itemIdsA = itemIdsA, itemIdsB = itemIdsB, ratingIds = ratingIds,
	     userIds = userIds, ratingMatrix = ratingMatrix,
	     d = problemInfo$d, nPseudoInputUsers = problemInfo$nPseudoInputUsers,
	     nPseudoInputItems = problemInfo$nPseudoInputItems, userFeatures = problemInfo$userFeatures)
}

addFITCauxiliarVariables <- function(problemInfo) {

	cholInverse <- function(M) rot180(forwardsolve(t(chol(rot180(M))), diag(nrow(M))))

	# First for the users

	fitcUsers <- list()
	fitcUsers$X <- problemInfo$userFeatures
	fitcUsers$Xbar <- problemInfo$userFeatures[ problemInfo$indexPseudoInputUsers, ]
	fitcUsers$m <- nrow(fitcUsers$Xbar)
	fitcUsers$n <- nrow(fitcUsers$X)
	fitcUsers$X <- fitcUsers$X / problemInfo$lengthScaleUsers
	fitcUsers$Xbar <- fitcUsers$Xbar / problemInfo$lengthScaleUsers
	fitcUsers$Km <- spgpComputeKm(fitcUsers$Xbar)
	fitcUsers$Kmn <- spgpComputeKmn(fitcUsers$X, fitcUsers$Xbar)
	fitcUsers$P <- t(fitcUsers$Kmn)
	fitcUsers$R <- cholInverse(fitcUsers$Km)
	fitcUsers$PRt <- fitcUsers$P %*% t(fitcUsers$R)
	fitcUsers$diagKn <- spgpComputeDiagKn(fitcUsers$X)
	fitcUsers$D <- fitcUsers$diagKn - rep(1, fitcUsers$m) %*% (fitcUsers$R %*% t(fitcUsers$P))^2

	# Second for the items

	fitcItems <- list()
	fitcItems$X <- cbind(problemInfo$itemFeaturesA, problemInfo$itemFeaturesB)
	fitcItems$Xbar <- fitcItems$X[ problemInfo$indexPseudoInputItems, ]
	fitcItems$m <- nrow(fitcItems$Xbar)
	fitcItems$n <- nrow(fitcItems$X)
	fitcItems$X <- fitcItems$X / problemInfo$lengthScaleItems
	fitcItems$Xbar <- fitcItems$Xbar / problemInfo$lengthScaleItems
	fitcItems$Km <- spgpComputeKmPreference(fitcItems$Xbar)
	fitcItems$Kmn <- spgpComputeKmnPreference(fitcItems$X, fitcItems$Xbar)
	fitcItems$P <- t(fitcItems$Kmn)
	fitcItems$R <- cholInverse(fitcItems$Km)
	fitcItems$PRt <- fitcItems$P %*% t(fitcItems$R)
	fitcItems$diagKn <- spgpComputeDiagKnPreference(fitcItems$X)
	fitcItems$D <- fitcItems$diagKn - rep(1, fitcItems$m) %*% (fitcItems$R %*% t(fitcItems$P))^2

	# We add the information to problemInfo

	problemInfo$fitcUsers <- fitcUsers
	problemInfo$fitcItems <- fitcItems

	problemInfo
}

##
# Function that implements the EP-VB algorithm for the multitask model
#
# @param	problemInfo	A list with the following fields:
#				
#				itemFeaturesA		-> n x q matrix with the features of the first rated item.
#				itemFeaturesB		-> n x q matrix with the features of the second rated item.
#				itemIdsA		-> n-dimensional vector with the id of the first item compared by the users.
#				itemIdsB		-> n-dimensional vector with the id of the first item compared by the users.
#				userIds			-> n-dimensional vector with the ids of the users.
#				ratings			-> n-dimensional vector with the ratings given by the users.
#				d			-> Number of latent functions.
#				nTotalItems		-> Number of total items.
#				userFeatures		-> U x p matrix with the user features.
#				nPseudoInputUsers	-> Number of pseudo inputs for the user weights.
#				nPseudoInputItems	-> Number of pseudo inputs for the latent functions.
#

epMPL <- function(problemInfo, lengthScaleItems, lengthScaleUsers) {

	# We preprocess the problemInfo structure

	problemInfo <- preprocessProblemInfo(problemInfo)

	# We select the pseudo inputs for the users and for the items

	problemInfo$indexPseudoInputUsers <- sample(1 : length(unique(problemInfo$userIds)), problemInfo$nPseudoInputUsers)
	problemInfo$indexPseudoInputItems <- sample(1 : length(unique(problemInfo$ratingIds)), problemInfo$nPseudoInputItems)

	# Target function to optimize the evidence

	target <- function(x) {
		problemInfo$lengthScaleItems <- exp(x[ 1 ])
		problemInfo$lengthScaleUsers <- exp(x[ 2 ])
		-epMPLinternal(problemInfo)$evidence
	}

#	ret <- optim(c(0, 0), target, method = "Nelder-Mead", control = list(trace = T, maxit = 20))

#	problemInfo$lengthScaleItems <- exp(ret$par[ 1 ])
#	problemInfo$lengthScaleUsers <- exp(ret$par[ 2 ])

	problemInfo$lengthScaleItems <- lengthScaleItems
	problemInfo$lengthScaleUsers <- lengthScaleUsers

	ret <- epMPLinternal(problemInfo)

	ret
}

##
# Function that computes the kernel matrix for the pseudo inputs.
#
# @param Xbar        -> Pseudo inputs (m x d matrix).
#
# @return -> Kernel matrix (m x m matrix).
#

spgpComputeKmPreference <- function(Xbar) {

	# We obtain the features of the different items

	m <- nrow(Xbar)
	d <- ncol(Xbar) / 2
	itemsA <- matrix(Xbar[ , 1 : d ], m, d)
	itemsB <- matrix(Xbar[ , (d + 1) : (2 * d) ], m, d)

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

	K + diag(1e-5, m)
}

##
# Function that computes the kernel matrix for the pseudo inputs and the data points.
#
# @param X    -> Features (n x d matrix).
# @param Xbar -> Pseudo inputs (m x d matrix).
#
# @return -> Kernel matrix (m x n matrix).
#

spgpComputeKmnPreference <- function(X, Xbar) {

	# We obtain the features of the different items

	n <- nrow(X)
	d <- ncol(X) / 2
	itemsA <- matrix(X[ , 1 : d ], n, d)
	itemsB <- matrix(X[ , (d + 1) : (2 * d) ], n, d)

	m <- nrow(Xbar)
	itemsAbar <- matrix(Xbar[ , 1 : d ], m, d)
	itemsBbar <- matrix(Xbar[ , (d + 1) : (2 * d) ], m, d)

	# We compute the kernel matrix

	QA <- matrix(apply(itemsA^2, 1, sum), m, n, byrow = T)
	QB <- matrix(apply(itemsB^2, 1, sum), m, n, byrow = T)

	QAbar <- matrix(apply(itemsAbar^2, 1, sum), m, n)
	QBbar <- matrix(apply(itemsBbar^2, 1, sum), m, n)

	distanceA <- QAbar + QA - 2 * itemsAbar %*% t(itemsA)
	KA <- exp(-0.5 * distanceA)

	distanceB <- QBbar + QB - 2 * itemsBbar %*% t(itemsB)
	KB <- exp(-0.5 * distanceB)

	distanceAB <- QAbar + QB - 2 * itemsAbar %*% t(itemsB)
	KAB <- exp(-0.5 * distanceAB)

	distanceBA <- QBbar + QA - 2 * itemsBbar %*% t(itemsA)
	KBA <- exp(-0.5 * distanceBA)

	KA + KB - KAB - KBA
}

##
# Function that computes the diagonal of the kernel matrix for the data points.
#
# @param X -> Features (n x d matrix).
#
# @return -> Diagonal of the kernel matrix (n-dimensional vector).
#

spgpComputeDiagKnPreference <- function(X) {

	n <- nrow(X)
	d <- ncol(X) / 2
	itemsA <- matrix(X[ , 1 : d ], n, d)
	itemsB <- matrix(X[ , (d + 1) : (2 * d) ], n, d)

	# We compute the kernel matrix using the preference kernel

	QA <- apply(itemsA^2, 1, sum)
	QB <- apply(itemsB^2, 1, sum)

	KA <- 1
	KB <- 1

	distanceAB <- QA + QB - 2 * apply(itemsA * itemsB, 1, sum)
	KAB <- exp(-0.5 * distanceAB)

	distanceBA <- QB + QA - 2 * apply(itemsB * itemsA, 1, sum)
	KBA <- exp(-0.5 * distanceBA)

	K <- KA + KB - KAB - KBA

	K + 1e-5
}

##
# Function that computes the kernel matrix for the pseudo inputs.
#
# @param Xbar        -> Pseudo inputs (m x d matrix).
#
# @return -> Kernel matrix (m x m matrix).
#

spgpComputeKm <- function(Xbar) {

	# We obtain the features of the different items

	m <- nrow(Xbar)
	Q <- matrix(apply(Xbar^2, 1, sum), m, m)

	distance <- Q + t(Q) - 2 * Xbar %*% t(Xbar)
	K <- exp(-0.5 * distance)

	K + diag(1e-5, m) + diag(0.1, nrow(K))
}

##
# Function that computes the kernel matrix for the pseudo inputs and the data points.
#
# @param X    -> Features (n x d matrix).
# @param Xbar -> Pseudo inputs (m x d matrix).
#
# @return -> Kernel matrix (m x n matrix).
#

spgpComputeKmn <- function(X, Xbar) {

	# We obtain the features of the different items

	n <- nrow(X)
	m <- nrow(Xbar)

	# We compute the kernel matrix

	Q <- matrix(apply(X^2, 1, sum), m, n, byrow = T)
	Qbar <- matrix(apply(Xbar^2, 1, sum), m, n)

	distance <- Qbar + Q - 2 * Xbar %*% t(X)
	K <- exp(-0.5 * distance)

	K
}

##
# Function that computes the diagonal of the kernel matrix for the data points.
#
# @param X -> Features (n x d matrix).
#
# @return -> Diagonal of the kernel matrix (n-dimensional vector).
#

spgpComputeDiagKn <- function(X) {

	n <- nrow(X)

	rep(1, n) + 1e-5 + 0.1
}
