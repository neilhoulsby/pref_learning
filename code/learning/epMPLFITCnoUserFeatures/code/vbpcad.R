library(SparseM)

library(irlba)

library(Matrix)

dyn.load("vbpcad.so")

##
# Function that initializes the current posterior approximation calling the svd method
#
# @param k       -> Number of latent features to use.
# @param mR      -> Sparse matrix with the mean of the ratings in csr format.
# @param mPprior -> Matrix with the mean componentes of the prior for the entries of matrix P.
# @param vPprior -> Matrix with the variance componentes of the prior for the entries of matrix P.
# @param mQprior -> Matrix with the mean componentes of the prior for the entries of matrix Q.
# @param vQprior -> Matrix with the variance componentes of the prior for the entries of matrix Q.
# @param nP      -> Number of rows in matrix P.
# @param nQ      -> Number of rows in matrix Q.
#
# @return        -> A list with the parameteres of the posterior approximation.

initializePA <- function(k, mR, mPprior, vPprior, mQprior, vQprior, nP, nQ) {

	# We initialize the learning rate

	lrate <- 1e-2

	# We initialize the posterior means calling the svd method.

	mR <- as.matrix.coo(mR)
	mR <- sparseMatrix(i = mR@ia, j = mR@ja, x = mR@ra, dims = mR@dimension)

	if (k + 3 <= min(nrow(mR), ncol(mR))) {
		ret <- irlba(mR, nu = k, nv = k)
		mP <- ret$u %*% diag(ret$d, k)
		mQ <- ret$v
	} else {
		mP <- matrix(rnorm(nP * k), nP, k) * sqrt(vPprior)
		mQ <- matrix(rnorm(nQ * k), nQ, k) * sqrt(vQprior)
	}

	# We initialize the posterior variances to the prior.

	vP <- vPprior
	vQ <- vQprior

	# We initialize the prior.

	mPprior <- mPprior
	vPprior <- vPprior
	mQprior <- mQprior
	vQprior <- vQprior

	# We return the posterior approximation.

	list(mP = mP, vP = vP, mQ = mQ, vQ = vQ, mPprior = mPprior, vPprior = vPprior, mQprior = mQprior, vQprior = vQprior, lrate = lrate, bound = -Inf)
}

##
# The VBPCAd method, calling the C routine for the optimization process.
#
# @param r  -> Matrix with the mean and variance of the ratings, one row per rating. Rating format (user_id, movie_id, rating mean, rating variance).
# @param k  -> Number of latent features to use.
# @param mPprior -> Matrix with the mean componentes of the prior for the entries of matrix P.
# @param vPprior -> Matrix with the variance componentes of the prior for the entries of matrix P.
# @param mQprior -> Matrix with the mean componentes of the prior for the entries of matrix Q.
# @param vQprior -> Matrix with the variance componentes of the prior for the entries of matrix Q.
# @param pa -> Initial point for the optimization process.
#
# @return   -> A list with the parameters of the posterior approximation.
#

bvpcadFast <- function(r, k, mPprior, vPprior, mQprior, vQprior, optimizePrior = 0, pa = NULL) {

	# We store the mean and the variance of the ratings in a sparse matrix with dimension nUsers x nMovies.
	
	mR <- as.matrix.csr(new("matrix.coo", ia = as.integer(r[ , 1 ]),
		ja = as.integer(r[ , 2 ]), ra = r[ , 3 ], dimension = as.integer(c(max(r[ , 1 ]), max(r[ , 2 ])))))
	vR <- as.matrix.csr(new("matrix.coo", ia = as.integer(r[ , 1 ]),
		ja = as.integer(r[ , 2 ]), ra = r[ , 4 ], dimension = as.integer(c(max(r[ , 1 ]), max(r[ , 2 ])))))

	# We initialize the posterior approximation.

	if (is.null(pa))
		pa <- initializePA(k, mR, mPprior, vPprior, mQprior, vQprior, nrow(mR), ncol(mR))
	else {
		pa$mPprior <- mPprior
		pa$vPprior <- vPprior
		pa$mQprior <- mQprior
		pa$vQprior <- vQprior
	}

	# We initialize the auxiliary parameters to pass to the C optimizer.

	aux <- list(k = as.integer(k), nP = nrow(mR), nQ = ncol(mR), nR = length(mR@ra), optimizePrior = as.integer(optimizePrior))

	# We obtain the representation of the mean and the variance of the ratings in column and row compressed form.

	mRcsc <- as.matrix.csc(mR)
	mR <- list(csr_ra = mR@ra, csr_ja = mR@ja, csr_ia = mR@ia, csc_ra = mRcsc@ra, csc_ja = mRcsc@ja, csc_ia = mRcsc@ia)
	vRcsc <- as.matrix.csc(vR)
	vR <- list(csr_ra = vR@ra, csc_ra = vRcsc@ra)

	# We call the C optimizer.

	ret <- .Call("mainOptimizationLoop", mR, vR, pa, aux)

	ret
}
