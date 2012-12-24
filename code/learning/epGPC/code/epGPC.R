#
# R script that implements the EP method for multitask preference learning.
#
# Author: Jose Miguel Hernandez Lobato
#
# Date: 23 Dec 2011
#

epGPCinternal <- function(itemsA, itemsB, userIds, ratings, lengthScale, paStart = NULL) {

	# We use only the first two observations

	d <- 1

	nRatings <- length(ratings)
	nUsers <- length(unique(userIds))

	# We initialize the posterior approximation and the factor approximations

	# The first approximate factor

	if (is.null(paStart))
		f1Hat <- list(mu = rep(0, nRatings), vu = rep(Inf, nRatings))
	else
		f1Hat <- paStart$f1Hat

	# The second approximate factor

	if (is.null(paStart))
		f2Hat <- list(mu = matrix(0, d, nRatings), vu = matrix(Inf, d, nRatings))
	else
		f2Hat <- paStart$f2Hat

	# We compute the kernel matrix

	ret <- computeKernel(itemsA, itemsB, lengthScale)

	# The posterior approximation with all the interesting information

	pa <- list(mu = rep(0, nRatings), vu = rep(Inf, nRatings),
		f1Hat = f1Hat, f2Hat = f2Hat, 
		auxM1 = NULL, auxM2 = NULL, Sigma = NULL, B = NULL,
		itemsA = itemsA, itemsB = itemsB, ratings = ratings,
		userIds = userIds, lengthScale = lengthScale, d = d,
		K = ret$K, diffKdiffLengthScale = ret$diffKdiffLengthScale,
		invK = chol2inv(chol(ret$K)))

	##
	## We start refinint the factors for the first iteration
	## 

	# We refine the second approximate factor

	if (is.null(paStart)) {
		pa$f2Hat$vu <- matrix(diag(pa$K), d, nRatings, byrow = T)
		pa$vu <- pa$f2Hat$vu
	}

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

		# We refine the approximate term for the Gaussian process. See Rasmusen's book for numerical stable
		# implementations of the Gaussian process.
		
#		pa$B <- diag(rep(1, nRatings)) + matrix(sqrt(pa$f1Hat$vu^-1), nRatings, nRatings) * pa$K * matrix(sqrt(pa$f1Hat$vu^-1), nRatings, nRatings, byrow = T)
#		pa$auxM1 <- matrix(sqrt(pa$f1Hat$vu^-1), nRatings, nRatings) * chol2inv(chol(pa$B)) * matrix(sqrt(pa$f1Hat$vu^-1), nRatings, nRatings, byrow = T)
#		pa$auxM2 <- pa$auxM1 %*% pa$K
#		pa$Sigma <- pa$K - pa$K %*% pa$auxM2

		pa$Sigma <- chol2inv(chol(pa$invK + diag(as.double(pa$f1Hat$vu^-1), nRatings)))

		vuNew <- diag(pa$Sigma)
		muNew <- as.double(pa$Sigma %*% t(pa$f1Hat$mu * pa$f1Hat$vu^-1))

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

		damping <- damping * 0.95

		iteration <- iteration + 1
	}

	# We evaluate the evidence in a robust way

	pa$evidence <- computeEvidenceRobust(pa)
	
	# We compute the gradient of the evidence with respect to the kernel lengthscale
	# See Rasmusen's book.

	pa$B <- diag(rep(1, nRatings)) + matrix(sqrt(pa$f1Hat$vu^-1), nRatings, nRatings) * pa$K * matrix(sqrt(pa$f1Hat$vu^-1), nRatings, nRatings, byrow = T)
	pa$auxM1 <- matrix(sqrt(pa$f1Hat$vu^-1), nRatings, nRatings) * chol2inv(chol(pa$B)) * matrix(sqrt(pa$f1Hat$vu^-1), nRatings, nRatings, byrow = T)
	pa$auxM2 <- pa$auxM1 %*% pa$K
	b <- (diag(rep(1, nRatings)) - pa$auxM2) %*% t(pa$f1Hat$mu * pa$f1Hat$vu^-1)
	pa$gradient <- 0.5 * sum(((b %*% t(b) - pa$auxM1) * pa$diffKdiffLengthScale) %*% rep(1, nRatings))

	# We are done

	pa
}

##
# Function that computes the preference learning kernel matrix.
#

computeKernel <- function(itemsA, itemsB, lengthScale) {

	n <- nrow(itemsA)

	# We compute the kernel matrix

	QA <- matrix(apply(itemsA^2, 1, sum), n, n)
	distanceA <- QA + t(QA) - 2 * itemsA %*% t(itemsA)
	KA <- exp(-0.5 * distanceA / lengthScale^2)

	QB <- matrix(apply(itemsB^2, 1, sum), n, n)
	distanceB <- QB + t(QB) - 2 * itemsB %*% t(itemsB)
	KB <- exp(-0.5 * distanceB / lengthScale^2)

	distanceAB <- QA + t(QB) - 2 * itemsA %*% t(itemsB)
	KAB <- exp(-0.5 * distanceAB / lengthScale^2)

	distanceBA <- QB + t(QA) - 2 * itemsB %*% t(itemsA)
	KBA <- exp(-0.5 * distanceBA / lengthScale^2)

	K <- KA + KB - KAB - KBA

	K <- K + diag(rep(1e-5, n))

	# We compute the derivative of K with respect to the lengthScale parameter

	diffKAdiffLengthScale <- KA * distanceA * lengthScale^-3
	diffKBdiffLengthScale <- KB * distanceB * lengthScale^-3
	diffKABdiffLengthScale <- KAB * distanceAB * lengthScale^-3
	diffKBAdiffLengthScale <- KBA * distanceBA * lengthScale^-3

	diffKdiffLengthScale <- diffKAdiffLengthScale + diffKBdiffLengthScale -
		diffKABdiffLengthScale - diffKBAdiffLengthScale

	list(K = K, diffKdiffLengthScale = diffKdiffLengthScale)
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
# Function that generates a prediction for the preference of
# a user on a particular pair of items.
#

predictGPC <- function(pa, newItemA, newItemB) {

	# We compute the covariance between the new pair of items
	# and all the other pairs of items in the trainint set.

	k <- function(x, y) exp(-0.5 * sum((x - y) * (x - y)) / pa$lengthScale^2)

	kStarA <- apply(pa$itemsA, 1, function(x) k(newItemA, x))
	kStarB <- apply(pa$itemsB, 1, function(x) k(newItemB, x))
	kStarAB <- apply(pa$itemsB, 1, function(x) k(newItemA, x))
	kStarBA <- apply(pa$itemsA, 1, function(x) k(newItemB, x))

	kStar <- kStarA + kStarB - kStarAB - kStarBA

	kNewNew <- k(newItemA, newItemA) + k(newItemB, newItemB) - k(newItemA, newItemB) - k(newItemA, newItemB) + 1e-5

	# Formulas for prediction from Rasmussen's book

	m <- t(kStar) %*% (diag(length(kStar)) - pa$auxM2) %*% t(pa$f1Hat$mu * pa$f1Hat$vu^-1)
	v <- as.double(kNewNew - t(kStar) %*% pa$auxM1 %*% kStar)

	# Given the new mean and variance vectors we generate a prediction for the user

	c(m, v)
}

##
# Function that computes the EP approximation of the evidence
#

computeEvidenceRobust <- function(pa) {

	# See Rasmussen's book

	logZ <- pnorm(pa$ratings * pa$f2Hat$mu / sqrt(pa$f2Hat$vu + 1), log.p = T)

	n <- nrow(pa$K)

	L <- chol(diag(rep(1, n)) + matrix(sqrt(pa$f1Hat$vu^-1), n, n) * pa$K * matrix(sqrt(pa$f1Hat$vu^-1), n, n, byrow = T))

	fourthAndFirst <- 0.5 * sum(log(1 + pa$f1Hat$vu^-1 * pa$f2Hat$vu)) - sum(log(diag(L)))
	cuadraticFifthAndSecond <- 0.5 * (pa$f1Hat$mu * pa$f1Hat$vu^-1) %*%
		(pa$Sigma - diag(as.double((pa$f2Hat$vu^-1 + pa$f1Hat$vu^-1)^-1))) %*% t(pa$f1Hat$mu * pa$f1Hat$vu^-1)
	finalFifthTerm <- sum(0.5 * pa$f2Hat$mu * pa$f2Hat$vu^-1 * pa$vu *
		(pa$f1Hat$vu^-1 * pa$f2Hat$mu - 2 * pa$f1Hat$mu * pa$f1Hat$vu^-1))

	sum(logZ) + fourthAndFirst + cuadraticFifthAndSecond + finalFifthTerm
}

##
# Function that fits the model selecting the kernelwidth that maximizies the evidence using gradient descent.
#

epGPCgradient <- function(itemsA, itemsB, userIds, ratings) {

	convergence <- F
	stepSize <- 1e-3
	previousEvidence <- -Inf
	sigma <- 1
	ret <- NULL
	i <- 0
	while (!convergence && i < 25) {

		# We evaluate the gradient

		ret <- epGPCinternal(itemsA, itemsB, userIds, ratings, sigma, ret)

		cat("Evidence:", ret$evidence, "LengthScale:", sigma, "Tolerance:", ret$evidence - previousEvidence, "\n")

		if (abs(ret$evidence - previousEvidence) < 0.0001) {
			cat("Converged")
			convergence <- T
		} 

		if (ret$evidence > previousEvidence)
			stepSize <- stepSize * 1.1
		else
			stepSize <- stepSize * 0.5

		previousEvidence <- ret$evidence

		sigma <- sigma + stepSize * ret$gradient

		i <- i + 1
	}

	ret
}

##
# Function that fits the model selecting the kernelwidth that maximizies the evidence using the optim function
#

epGPCgrid <- function(itemsA, itemsB, userIds, ratings) {

	gridLengthScale <- exp(seq(-2, 2, length.out = 10))
	evidence <- rep(0, length(gridLengthScale))
	model <- list()
	for (i in 1 : length(gridLengthScale)) {
		model[[ i ]] <- epGPCinternal(itemsA, itemsB, userIds, ratings, gridLengthScale[ i ])
		evidence[ i ] <- model[[ i ]]$evidence
	}

	model[[ which.max(evidence) ]]
}

##
# Function that generates a prediction for the preference of a user on several pairs of items
# in a vectorized form.
#

predictGPCvectorized <- function(pa, newItemsA, newItemsB) {

	# We compute the covariance between the new pair of items
	# and all the other pairs of items in the trainint set.

	distanceA <- matrix(apply(pa$itemsA^2, 1, sum), nrow(pa$itemsA), nrow(newItemsA)) +
		matrix(apply(newItemsA^2, 1, sum), nrow(pa$itemsA), nrow(newItemsA), byrow = T) - 2 * pa$itemsA %*% t(newItemsA)
	KA <- exp(-0.5 * distanceA / pa$lengthScale^2)
	KA <- KA

	distanceB <- matrix(apply(pa$itemsB^2, 1, sum), nrow(pa$itemsB), nrow(newItemsB)) +
		matrix(apply(newItemsB^2, 1, sum), nrow(pa$itemsB), nrow(newItemsB), byrow = T) - 2 * pa$itemsB %*% t(newItemsB)
	KB <- exp(-0.5 * distanceB / pa$lengthScale^2)
	KB <- KB

	distanceAB <- matrix(apply(pa$itemsA^2, 1, sum), nrow(pa$itemsA), nrow(newItemsB)) +
		matrix(apply(newItemsB^2, 1, sum), nrow(pa$itemsA), nrow(newItemsB), byrow = T) - 2 * pa$itemsA %*% t(newItemsB)
	KAB <- exp(-0.5 * distanceAB / pa$lengthScale^2)
	KAB <- KAB

	distanceBA <- matrix(apply(pa$itemsB^2, 1, sum), nrow(pa$itemsB), nrow(newItemsA)) +
		matrix(apply(newItemsA^2, 1, sum), nrow(pa$itemsB), nrow(newItemsA), byrow = T) - 2 * pa$itemsB %*% t(newItemsA)
	KBA <- exp(-0.5 * distanceBA / pa$lengthScale^2)
	KBA <- KBA

	kStar <- KA + KB - KAB - KBA

	distanceNewNew <- apply(newItemsA^2, 1, sum) + apply(newItemsB^2, 1, sum) - 2 * apply(newItemsA * newItemsB, 1, sum)
	kNewNew <- exp(-0.5 * distanceNewNew / pa$lengthScale^2)

	kNewNew <- 1 + 1 - 2 * kNewNew + 1e-5

	# Formulas for prediction from Rasmussen's book

	m <- as.double(t(kStar) %*% (diag(nrow(kStar)) - pa$auxM2) %*% t(pa$f1Hat$mu * pa$f1Hat$vu^-1))
	v <- as.double(kNewNew - apply((t(kStar) %*% pa$auxM1) * t(kStar), 1, sum))

	list(m = m, v = v)
}
