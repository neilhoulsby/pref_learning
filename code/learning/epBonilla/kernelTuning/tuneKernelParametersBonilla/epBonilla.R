#
# R script that implements the EP method for multitask preference learning.
#
# Author: Jose Miguel Hernandez Lobato
#
# Date: 23 Dec 2011
#

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

epBonillaInternal <- function(problemInfo, initialization = NULL) {

	# We add the additional variables to the problemInfo structure

	problemInfo <- addAuxiliarVariables(problemInfo)

	# Total number of different pairs of items rated and the number of users

	nRatings <- length(problemInfo$ratings)

	# We initialize the posterior approximation and the factor approximations

	# The first approximate factor

	f1Hat <- list(mu = rep(0, nRatings), vu = rep(Inf, nRatings))

	# The second approximate factor

	f2Hat <- list(mu = rep(0, nRatings), vu = rep(Inf, nRatings))

	# The posterior approximation

	a <- list(mu = rep(0, nRatings), vu = rep(Inf, nRatings))

	##
	## We start refinint the factors for the first iteration
	## 

	# We refine the second approximate factor

	a$vu <- f2Hat$vu <- diag(problemInfo$K)

	# We check for an initial solution

	if (!is.null(initialization)) {
		f1Hat <- initialization$f1Hat
		f2Hat <- initialization$f2Hat
		a <- initialization$a
	}

	##
	## Main loop of EP
	##

	# We check for an initial solution

	damping <- 0.5
	convergence <- F
	iteration <- 1
	while ((((!convergence && iteration <= 100) || iteration <= 1))) {

		aOld <- a

		##
		## We refine the first approximate factor
		##

		mOld <- f2Hat$mu
		vOld <- f2Hat$vu
		Y <- problemInfo$ratings

		logZ <- pnorm(Y * mOld / sqrt(vOld + 1), log.p = T)
		ratio <- exp(-logZ + dnorm(mOld / sqrt(vOld + 1), log = T))
		alpha <- ratio * Y / sqrt(vOld + 1)
		beta <- -ratio * (Y * mOld / sqrt(vOld + 1) + ratio) / (1 + vOld)
                eta2HatNew <- -beta / (1 + beta * vOld)
                eta1HatNew <- (alpha - mOld * beta) / (1 + beta * vOld)

		vuHatNew <- eta2HatNew^-1
		muHatNew <- eta1HatNew / eta2HatNew

			# We do damping

		f1Hat$mu <- damping * muHatNew / vuHatNew + (1 - damping) * f1Hat$mu / f1Hat$vu
		f1Hat$vu <- (damping * vuHatNew^-1 + (1 - damping) * f1Hat$vu^-1)^-1
		f1Hat$mu <- f1Hat$vu * f1Hat$mu

		##
	      	## We refine the second approximate factor
		##

		# We refine the approximate factor for the Gaussian process using the FITC approximation

		Sigma <- solve(problemInfo$invK + diag(f1Hat$vu^-1, length(f1Hat$vu)))
		muNew <- as.double(Sigma %*% (f1Hat$mu / f1Hat$vu))
		vuNew <- diag(Sigma)

		# We update the fourth approximate factor
		
		vuHatNew <- 1 / (1 / vuNew - 1 / f1Hat$vu)
		muHatNew <- vuHatNew * (muNew / vuNew - f1Hat$mu / f1Hat$vu)
			
		# We do damping

		f2Hat$mu <- damping * muHatNew / vuHatNew + (1 - damping) * f2Hat$mu / f2Hat$vu
		f2Hat$vu <- (damping * vuHatNew^-1 + (1 - damping) * f2Hat$vu^-1)^-1
		f2Hat$mu <- f2Hat$vu * f2Hat$mu

		# We update the posterior approximation

		a$vu <- 1 / (1 / f1Hat$vu + 1 / f2Hat$vu)
		a$mu <- a$vu * (f1Hat$mu / f1Hat$vu + f2Hat$mu / f2Hat$vu)

		# We update the damping parameter

		damping <- damping * 0.95

		# We check for convergence

		change <- max(abs(a$mu - aOld$mu))

		if (change < 1e-3)
			convergence <- T
		else 
			convergence <- F

		cat(iteration, change, "\n")

		iteration <- iteration + 1
	}

	# We estimate the evidence

	evidence <- computeEvidence(a, problemInfo, f1Hat, f2Hat)

	# Some auxiliary variables

	problemInfo$Minv <- solve(diag(rep(1, nRatings)) + problemInfo$K * matrix(f1Hat$vu^-1, nRatings, nRatings))

	# We compute the gradient of the evidence with respect to the kernel parameters

	ret <- computeGradEvidence(a, problemInfo, f1Hat, f2Hat)
	gradLengthScaleUsers <- ret$gradLengthScaleUsers
	gradLengthScaleItems <- ret$gradLengthScaleItems
	gradNoiseUsers <- ret$gradNoiseUsers
	gradNoiseItems <- ret$gradNoiseItems

	# We return the posterior approximation

	list(a = a, evidence = evidence, problemInfo = problemInfo, f1Hat = f1Hat, f2Hat = f2Hat,
		gradLengthScaleUsers = gradLengthScaleUsers, gradLengthScaleItems = gradLengthScaleItems, gradNoiseUsers = gradNoiseUsers,
		gradNoiseItems = gradNoiseItems)
}

##
# Function that evaluates the EP approximation of the evidence
#

computeEvidence <- function(a, problemInfo, f1Hat, f2Hat) {

	# See Rasmussen's book

	logZ <- pnorm(problemInfo$ratings * f2Hat$mu / sqrt(f2Hat$vu + 1), log.p = T)

	# Some required variables

	K <- problemInfo$K
	invK <- problemInfo$invK
	f1Hateta2 <- f1Hat$vu^-1
	f1Hateta1 <- f1Hat$mu * f1Hat$vu^-1
	f2Hateta2 <- f2Hat$vu^-1
	f2Hateta1 <- f2Hat$mu * f2Hat$vu^-1
	aeta2 <- a$vu^-1
	aeta1 <- a$mu * a$vu^-1

	n <- nrow(K)

	Sigma <- solve(invK + diag(as.double(f1Hateta2), n))

	m <- as.double(Sigma %*% f1Hateta1)

	M <- matrix(f1Hateta2, n, n) * K + diag(rep(1, n))

	logZ <- 0.5 * as.double(t(m) %*% (invK + diag(as.double(f1Hateta2), n)) %*% m) -
		0.5 * determinant(M)$modulus[[ 1 ]] + sum(logZ) +
		sum(0.5 * log(f1Hateta2 / f2Hateta2 + 1) - 0.5 * aeta1^2 / aeta2 + 0.5 * f2Hateta1^2 / f2Hateta2)

	logZ
}

##
# Function that computes the gradient of the EP approximation of the evidence
# with respect to the kernel hyper-parameters.
#

computeGradEvidence <- function(a, problemInfo, f1Hat, f2Hat) {

	gradLengthScaleUsers <- rep(0, problemInfo$d1)
	gradLengthScaleItems <- rep(0, problemInfo$d2)
	gradNoiseUsers <- 0
	gradNoiseItems <- 0

	# We compute the gradient with respect to the user lenthScale parameters

	K <- problemInfo$K
	invK <- problemInfo$invK
	f1Hateta2 <- f1Hat$vu^-1
	f1Hateta1 <- f1Hat$mu * f1Hat$vu^-1
	f2Hateta2 <- f2Hat$vu^-1
	f2Hateta1 <- f2Hat$mu * f2Hat$vu^-1
	aeta2 <- a$vu^-1
	aeta1 <- a$mu * a$vu^-1

	n <- nrow(K)
	Minv <- solve(diag(rep(1, n)) + K * matrix(f1Hateta2, n, n))
	mAuxPost <- as.double(Minv %*% f1Hateta1)

	for (i in 1 : problemInfo$d1)
		gradLengthScaleUsers[ i ] <- as.double(-0.5 * sum(Minv * t(matrix(f1Hateta2, n, n) * problemInfo$dKdLengthScaleUsers[[ i ]])) +
			0.5 * t(mAuxPost) %*% problemInfo$dKdLengthScaleUsers[[ i ]] %*% mAuxPost)

	for (i in 1 : problemInfo$d2)
		gradLengthScaleItems[ i ] <- as.double(-0.5 * sum(Minv * t(matrix(f1Hateta2, n, n) * problemInfo$dKdLengthScaleItems[[ i ]])) +
			0.5 * t(mAuxPost) %*% problemInfo$dKdLengthScaleItems[[ i ]] %*% mAuxPost)
			
	gradNoiseUsers <- as.double(-0.5 * sum(Minv * t(matrix(f1Hateta2, n, n) * problemInfo$dKdNoiseUsers)) +
		0.5 * t(mAuxPost) %*% problemInfo$dKdNoiseUsers %*% mAuxPost)

	gradNoiseItems <- as.double(-0.5 * sum(Minv * t(matrix(f1Hateta2, n, n) * problemInfo$dKdNoiseItems)) +
		0.5 * t(mAuxPost) %*% problemInfo$dKdNoiseItems %*% mAuxPost)
	
	list(gradLengthScaleUsers = gradLengthScaleUsers,
		gradLengthScaleItems = gradLengthScaleItems, gradNoiseUsers = gradNoiseUsers, gradNoiseItems = gradNoiseItems)
}

##
# Function that generates a prediction for the preference of a user on several pairs of items
# in a vectorized form.
#

predictBonilla <- function(ret, userFeatures, itemFeaturesA, itemFeaturesB) {

	# We compute the FITC prediction

	d1 <- ret$problemInfo$d1
	d2 <- ret$problemInfo$d2
	Xtest <- cbind(userFeatures, itemFeaturesA, itemFeaturesB)
	pStar <- t(spgpComputeKmnPreference(Xtest, ret$problemInfo$X, d1, d2, ret$problemInfo$lengthScaleItems, ret$problemInfo$lengthScaleUsers))
	dStar <- spgpComputeDiagKnPreference(Xtest, d1, d2, ret$problemInfo$lengthScaleItems,
		ret$problemInfo$lengthScaleUsers, ret$problemInfo$noiseUsers, ret$problemInfo$noiseItems)

	invK <- ret$problemInfo$invK
	f1Hateta2 <- ret$f1Hat$vu^-1
	f1Hateta1 <- ret$f1Hat$mu * ret$f1Hat$vu^-1
	n <- nrow(ret$problemInfo$K)

	m <- as.double(pStar %*% ret$problemInfo$Minv %*% f1Hateta1)
	v <- as.double(dStar - ((pStar %*% (ret$problemInfo$Minv * matrix(f1Hateta2, n, n, byrow = T))) * pStar) %*% rep(1, ncol(pStar)))

	list(m = m, v = v)
}

addAuxiliarVariables <- function(problemInfo) {

	problemInfo$lengthScaleItems <- exp(problemInfo$lengthScaleItems)
	problemInfo$lengthScaleUsers <- exp(problemInfo$lengthScaleUsers)
	problemInfo$noiseUsers <- exp(problemInfo$noiseUsers)
	problemInfo$noiseItems <- exp(problemInfo$noiseItems)

	problemInfo$d1 <- ncol(problemInfo$userFeatures)
	problemInfo$d2 <- ncol(problemInfo$itemFeaturesA)

	problemInfo$X <- cbind(problemInfo$userFeatures, problemInfo$itemFeaturesA, problemInfo$itemFeaturesB)

	ret <- spgpComputeKmPreference(problemInfo$X, problemInfo$d1,
		problemInfo$d2, problemInfo$lengthScaleItems, problemInfo$lengthScaleUsers, problemInfo$noiseUsers, problemInfo$noiseItems)

	problemInfo$K <- ret$K
	problemInfo$dKdLengthScaleItems <- ret$dKdLengthScaleItems
	problemInfo$dKdLengthScaleUsers <- ret$dKdLengthScaleUsers
	problemInfo$dKdNoiseUsers <- ret$dKdNoiseUsers
	problemInfo$dKdNoiseItems <- ret$dKdNoiseItems
	problemInfo$invK <- chol2inv(chol(problemInfo$K))

	problemInfo
}

##
#

epBonilla <- function(problemInfo) {

	# We initialize the lengthScale parameters

	itemFeatures <- rbind(problemInfo$itemFeaturesA, problemInfo$itemFeaturesB)
	n <- nrow(itemFeatures)
	Q <- matrix(apply(itemFeatures^2, 1, sum), n, n)
	distance <- Q + t(Q) - 2 * itemFeatures %*% t(itemFeatures)
	problemInfo$lengthScaleItems <- rep(log(sqrt(0.5 * (median(distance[ upper.tri(distance) ])))), ncol(itemFeatures))
	problemInfo$noiseItems <- -2

	userFeatures <- problemInfo$userFeatures
	n <- nrow(userFeatures)
	Q <- matrix(apply(userFeatures^2, 1, sum), n, n)
	distance <- Q + t(Q) - 2 * userFeatures %*% t(userFeatures)
	problemInfo$lengthScaleUsers <- rep(log(sqrt(0.5 * (median(distance[ upper.tri(distance) ])))), ncol(userFeatures))
	problemInfo$noiseUsers <- -2

	ret <- epBonillaInternal(problemInfo)

	cat("\n", 0, "New evidence:", ret$evidence, "\n\n")

	eps <- 0.1
	convergence <- F
	iteration <- 1
	while (!convergence && iteration < 100) {

		problemInfo$lengthScaleUsers <- problemInfo$lengthScaleUsers + eps * ret$gradLengthScaleUsers
		problemInfo$lengthScaleItems <- problemInfo$lengthScaleItems + eps * ret$gradLengthScaleItems
		problemInfo$noiseUsers <- problemInfo$noiseUsers + eps * ret$gradNoiseUsers
		problemInfo$noiseItems <- problemInfo$noiseItems + eps * ret$gradNoiseItems

		retNew <- epBonillaInternal(problemInfo, ret)

		if (is.nan(retNew$evidence))
			return(ret)

		if (abs(retNew$evidence - ret$evidence) < 1e-2)
			convergence <- T

		if (retNew$evidence < ret$evidence)
			eps <- eps * 0.5
		else
			eps <- eps * 1.1

		cat("\n", iteration, "New evidence:", retNew$evidence, "eps:", eps, "\n\n")

		ret <- retNew

		iteration <- iteration + 1
	}

	ret
}

##
# Function that computes the kernel matrix for the pseudo inputs.
#
# @param Xbar        -> Pseudo inputs (m x d matrix).
#
# @return -> Kernel matrix (m x m matrix).
#

spgpComputeKmPreference <- function(Xbar, d1, d2, lengthScaleItems, lengthScaleUsers, noiseUsers, noiseItems) {

	# We obtain the features of the different items

	m <- nrow(Xbar)
	userFeatures <- matrix(Xbar[ , 1 : d1 ], m, d1)
	userFeatures <- userFeatures / matrix(lengthScaleUsers, m, d1, byrow = T)

	itemsA <- matrix(Xbar[ , (d1 + 1) : (d1 + d2) ], m, d2)
	itemsA <- itemsA / matrix(lengthScaleItems, m, d2, byrow = T)
	itemsB <- matrix(Xbar[ , (d1 + d2 + 1) : (d1 + 2 * d2) ], m, d2)
	itemsB <- itemsB / matrix(lengthScaleItems, m, d2, byrow = T)

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

	Kitems <- KA + KB - KAB - KBA + diag(noiseItems, m)

	Q <- matrix(apply(userFeatures^2, 1, sum), m, m)
	distance <- Q + t(Q) - 2 * userFeatures %*% t(userFeatures)
	Kusers <- exp(-0.5 * distance) + diag(noiseUsers, m)

	K <- Kusers * (Kitems + diag(1e-4, m))

	# We compute the derivatives of the kernel matrix with respect to the lengthscale parameters

	dKdLengthScaleItems <- list()
	for (i in 1 : d2) {
		QA <- matrix(itemsA[ , i ]^2, m, m)
		dKA <- KA * (QA + t(QA) - 2 * itemsA[ , i ] %*% t(itemsA[ , i ]))
		QB <- matrix(itemsB[ , i ]^2, m, m)
		dKB <- KB * (QB + t(QB) - 2 * itemsB[ , i ] %*% t(itemsB[ , i ]))
		dKAB <- KAB * (QA + t(QB) - 2 * itemsA[ , i ] %*% t(itemsB[ , i ]))
		dKBA <- KBA * (QB + t(QA) - 2 * itemsB[ , i ] %*% t(itemsA[ , i ]))
		dKdLengthScaleItems[[ i ]] <- Kusers * (dKA + dKB - dKAB - dKBA)
	}

	dKdLengthScaleUsers <- list()
	for (i in 1 : d1) {
		Q <- matrix(userFeatures[ , i ]^2, m, m)
		dKdLengthScaleUsers[[ i ]] <- (Kitems + diag(1e-4, m)) *
			Kusers * (Q + t(Q) - 2 * userFeatures[ , i ] %*% t(userFeatures[ , i ]))
	}

	dKdNoiseUsers <- (Kitems + diag(1e-4, m)) * diag(noiseUsers, m)

	dKdNoiseItems <- Kusers * diag(noiseItems, m)

	list(K = K, dKdLengthScaleItems = dKdLengthScaleItems,
		dKdLengthScaleUsers = dKdLengthScaleUsers, dKdNoiseUsers = dKdNoiseUsers, dKdNoiseItems = dKdNoiseItems)
}

##
# Function that computes the kernel matrix for the pseudo inputs and the data points.
#
# @param X    -> Features (n x d matrix).
# @param Xbar -> Pseudo inputs (m x d matrix).
#
# @return -> Kernel matrix (m x n matrix).
#

spgpComputeKmnPreference <- function(X, Xbar, d1, d2, lengthScaleItems, lengthScaleUsers) {

	# We obtain the features of the different items

	n <- nrow(X)
	userFeatures <- matrix(X[ , 1 : d1 ], n, d1)
	userFeatures <- userFeatures / matrix(lengthScaleUsers, n, d1, byrow = T)
	itemsA <- matrix(X[ , (d1 + 1) : (d1 + d2) ], n, d2)
	itemsA <- itemsA / matrix(lengthScaleItems, n, d2, byrow = T)
	itemsB <- matrix(X[ , (d1 + d2 + 1) : (d1 + 2 * d2) ], n, d2)
	itemsB <- itemsB / matrix(lengthScaleItems, n, d2, byrow = T)

	m <- nrow(Xbar)
	userFeaturesBar <- matrix(Xbar[ , 1 : d1 ], m, d1)
	userFeaturesBar <- userFeaturesBar / matrix(lengthScaleUsers, m, d1, byrow = T)
	itemsAbar <- matrix(Xbar[ , (d1 + 1) : (d1 + d2) ], m, d2)
	itemsAbar <- itemsAbar / matrix(lengthScaleItems, m, d2, byrow = T)
	itemsBbar <- matrix(Xbar[ , (d1 + d2 + 1) : (d1 + 2 * d2) ], m, d2)
	itemsBbar <- itemsBbar / matrix(lengthScaleItems, m, d2, byrow = T)

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

	Kitems <- KA + KB - KAB - KBA

	Q <- matrix(apply(userFeatures^2, 1, sum), m, n, byrow = T)
	Qbar <- matrix(apply(userFeaturesBar^2, 1, sum), m, n)
	distance <- Qbar + Q - 2 * userFeaturesBar %*% t(userFeatures)
	Kusers <- exp(-0.5 * distance)

	Kusers * Kitems
}

##
# Function that computes the diagonal of the kernel matrix for the data points.
#
# @param X -> Features (n x d matrix).
#
# @return -> Diagonal of the kernel matrix (n-dimensional vector).
#

spgpComputeDiagKnPreference <- function(X, d1, d2, lengthScaleItems, lengthScaleUsers, noiseUsers, noiseItems) {

	m <- nrow(X)
	userFeatures <- matrix(X[ , 1 : d1 ], m, d1)
	userFeatures <- userFeatures / matrix(lengthScaleUsers, m, d1, byrow = T)
	itemsA <- matrix(X[ , (d1 + 1) : (d1 + d2) ], m, d2)
	itemsA <- itemsA / matrix(lengthScaleItems, m, d2, byrow = T)
	itemsB <- matrix(X[ , (d1 + d2 + 1) : (d1 + 2 * d2) ], m, d2)
	itemsB <- itemsB / matrix(lengthScaleItems, m, d2, byrow = T)

	# We compute the kernel matrix using the preference kernel

	QA <- apply(itemsA^2, 1, sum)
	QB <- apply(itemsB^2, 1, sum)

	KA <- 1
	KB <- 1

	distanceAB <- QA + QB - 2 * apply(itemsA * itemsB, 1, sum)
	KAB <- exp(-0.5 * distanceAB)

	distanceBA <- QB + QA - 2 * apply(itemsB * itemsA, 1, sum)
	KBA <- exp(-0.5 * distanceBA)

	Kitems <- KA + KB - KAB - KBA + noiseItems

	Kitems * (1 + noiseUsers) + 1e-4
}
