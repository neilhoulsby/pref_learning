#
# Script that implements the operations required to run EP with the FTIC approximation
#

##
# Function that makes predictions using the FITC approximation.
#

predictFITC <- function(D0, P0, R0, P0R0t, mHat, vHat, dStar, pStar) {

	n <- nrow(P0)
	m <- ncol(P0)

	tauTilde <- vHat^-1
	muTilde <- mHat / vHat

	Dnew <- D0 / (1 + D0 * tauTilde)
	Pnew <- matrix(1 / (1 + D0 * tauTilde), n, m) * P0

	Rnew <- backsolve(rot180(t(chol(rot180(diag(m) + t(P0R0t) %*% (matrix(tauTilde / (1 + D0 * tauTilde), n, m) * P0R0t))))), R0)

	gammaNew <- t(Rnew) %*% (Rnew %*% (t(Pnew) %*% muTilde))

	# We obtain the new marginals

	mPrediction <- pStar %*% gammaNew
	vPrediction <- dStar + rep(1, m) %*% (Rnew %*% t(pStar))^2

#	Knew <- diag(as.double(Dnew)) + (Pnew %*% t(Rnew)) %*% t(Pnew %*% t(Rnew))
#	K0 <- diag(as.double(D0)) + (P0 %*% t(R0)) %*% t(P0 %*% t(R0))
#	kStar <- P0 %*% t(R0) %*% R0 %*% t(pStar)
#	kStarStar <- dStar + diag(pStar %*% t(R0) %*% R0 %*% t(pStar))
#	mPrediction2 <- t(kStar) %*% solve(K0 + diag(as.double(vHat))) %*% mHat
#	vPrediction2 <- kStarStar - diag(t(kStar) %*% solve(K0 + diag(as.double(vHat))) %*% kStar)

	list(m = mPrediction, v = vPrediction)
}


##
# Function that evaluates the normalization constant of the product
# of the FITC prior and a multivariate Gaussian density with diagonal correlation matrix.
#

getFITCevidence <- function(D0, P0, R0, P0R0t, mHat, vHat) {

	n <- nrow(P0)
	m <- ncol(P0)

	tauTilde <- vHat^-1
	muTilde <- mHat / vHat

	Dnew <- D0 / (1 + D0 * tauTilde)
	Pnew <- matrix(1 / (1 + D0 * tauTilde), n, m) * P0

	Rnew <- backsolve(rot180(t(chol(rot180(diag(m) + t(P0R0t) %*% (matrix(tauTilde / (1 + D0 * tauTilde), n, m) * P0R0t))))), R0)

	aNew <- Dnew * muTilde
	gammaNew <- t(Rnew) %*% (Rnew %*% (t(Pnew) %*% muTilde))

	# We obtain the new marginals

	vNew <- as.double(Dnew + rep(1, m) %*% (Rnew %*% t(Pnew))^2)
	mNew <- as.double(aNew) + as.double(Pnew %*% gammaNew)

	# We obtain the cavity marginals and its natural parameters

#	Knew <- diag(as.double(Dnew)) + (Pnew %*% t(Rnew)) %*% t(Pnew %*% t(Rnew))
#	K0 <- diag(as.double(D0)) + (P0 %*% t(R0)) %*% t(P0 %*% t(R0))
#	logZ3 <- -n / 2 * log(2 * pi) + 0.5 * 2 * sum(log(diag(chol(Knew)))) - 0.5 * 2 * sum(log(diag(chol(K0)))) -
#		0.5 * sum(log(vHat)) - 0.5 * sum(mHat^2 / vHat) + 0.5 * t(mNew) %*% solve(Knew) %*% mNew

	logZ <- -n / 2 * log(2 * pi) + sum(log(diag(Rnew))) - sum(log(diag(R0))) -
		0.5 * sum(log(1 + tauTilde * D0)) - 0.5 * sum(log(vHat)) + 0.5 * sum(muTilde * mNew) - 0.5 * sum(mHat^2 / vHat)

	logZ
}


##
# Function that computes the marginal means and variances of the product
# of the FITC prior and a multivariate Gaussian density with diagonal correlation matrix.
#

computeTitledDistribution <- function(D0, P0, R0, P0R0t, mHat, vHat) {

	n <- nrow(P0)
	m <- ncol(P0)

	Ttilde <- vHat^-1
	muTilde <- mHat / vHat

	Dnew <- D0 / (1 + D0 * Ttilde)
	Pnew <- matrix(1 / (1 + D0 * Ttilde), n, m) * P0

	Rnew <- backsolve(rot180(t(chol(rot180(diag(m) + t(P0R0t) %*% (matrix(Ttilde / (1 + D0 * Ttilde), n, m) * P0R0t))))), R0)

	aNew <- Dnew * muTilde
	gammaNew <- t(Rnew) %*% (Rnew %*% (t(Pnew) %*% muTilde))

	# We obtain the new marginals

	vNew <- as.double(Dnew + rep(1, m) %*% (Rnew %*% t(Pnew))^2)
	mNew <- as.double(aNew) + as.double(Pnew %*% gammaNew)

	list(mNew = mNew, vNew = vNew)
}

##
# Function which rotates a matrix 180 degreees.
#

rot180 <- function(M) {

	matrix(rev(as.double(M)), nrow(M), ncol(M))
}
