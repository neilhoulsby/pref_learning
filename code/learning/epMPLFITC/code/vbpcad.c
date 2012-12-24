#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_PREVIOUS_BOUNDS 10

/* Auxiliary functions */

SEXP getListElement(SEXP list, const char *str);
double evaluateLowerBound(SEXP R_mR, SEXP R_vR, SEXP R_pa, SEXP R_aux, double *e, int *eP, int *eQ, double *avgTrE, double *avgPV);
void computeGradient(SEXP R_mR, SEXP R_vR, SEXP R_pa, SEXP R_aux, double *e, int *eP, int *eQ, double *grP, double *grQ);
void computeNewPosteriorVariances(SEXP R_mR, SEXP R_vR, SEXP R_pa, SEXP R_aux, double *e, int *eP, int *eQ, double *vPnew, double *vQnew);
void updatePosteriorVariances(SEXP R_ret, SEXP R_aux, double *vPnew, double *vQnew);
void restorePosteriorMeansAndVariances(SEXP R_ret, SEXP R_aux, double *vPold, double *vQold, double *grP, double *grQ);
void updatePosteriorMeans(SEXP R_ret, SEXP R_aux, double *vPnew, double *vQnew, double *grP, double *grQ, double lrate);
void updatePriorAndNoiseVariances(SEXP R_ret, SEXP R_aux);

double abs2(double x) {

	if (x < 0)
		x = -x;

	return x;
}

/**
 * Main function that implements the optimization process in C 
 *
 * @param R_mR  -> List with the representation of the means of the ratings as a sparse binary matrix.
 * @param R_vR  -> List with the representation of the variances of the ratings as a sparse binary matrix.
 * @param R_pa  -> The initial value of the posterior approximation.
 * @param R_aux -> List with auxiliary variables for the optimization process.
 *
 */

SEXP mainOptimizationLoop(SEXP R_mR, SEXP R_vR, SEXP R_pa, SEXP R_aux) {

	SEXP R_ret; 
	double bound, boundOld, lrate, avgTrE, avgPV, diffBounds;
	int convergence, iter, i;

	/* We reserve memory for auxiliary data structures */

	double *e = (double *) malloc(sizeof(double) * *INTEGER_POINTER(getListElement(R_aux, "nR")));
	int *eP = (int *) malloc(sizeof(int) * *INTEGER_POINTER(getListElement(R_aux, "nR")));
	int *eQ = (int *) malloc(sizeof(int) * *INTEGER_POINTER(getListElement(R_aux, "nR")));

	double *grP = (double *) malloc(sizeof(double) * *INTEGER_POINTER(getListElement(R_aux, "nP")) * *INTEGER_POINTER(getListElement(R_aux, "k")));
	double *grQ = (double *) malloc(sizeof(double) * *INTEGER_POINTER(getListElement(R_aux, "nQ")) * *INTEGER_POINTER(getListElement(R_aux, "k")));

	double *vPbackup = (double *) malloc(sizeof(double) * *INTEGER_POINTER(getListElement(R_aux, "nP")) * *INTEGER_POINTER(getListElement(R_aux, "k")));
	double *vQbackup = (double *) malloc(sizeof(double) * *INTEGER_POINTER(getListElement(R_aux, "nQ")) * *INTEGER_POINTER(getListElement(R_aux, "k")));

	double *previousBounds = (double *) malloc(sizeof(double) * MAX_PREVIOUS_BOUNDS);

	/* Protection macros for the parameters */

	PROTECT(R_mR = AS_LIST(R_mR));
	PROTECT(R_vR = AS_LIST(R_vR));
	PROTECT(R_pa = AS_LIST(R_pa));
	PROTECT(R_aux = AS_LIST(R_aux));

	/* We create a copy of the current approximation */

	PROTECT(R_ret = duplicate(R_pa));

	/* We evaluate the lower bound */

//	fprintf(stdout, "Starting the optimization process.\n");
//	fflush(stdout);
	boundOld = evaluateLowerBound(R_mR, R_vR, R_ret, R_aux, e, eP, eQ, &avgTrE, &avgPV);

	/* We update the posterior variances */

	computeNewPosteriorVariances(R_mR, R_vR, R_ret, R_aux, e, eP, eQ, vPbackup, vQbackup);
	updatePosteriorVariances(R_ret, R_aux, vPbackup, vQbackup);

	/* We evaluate the lower bound */

	boundOld = evaluateLowerBound(R_mR, R_vR, R_ret, R_aux, e, eP, eQ, &avgTrE, &avgPV);

	/* Main loop of the algorithm */

	lrate = *NUMERIC_POINTER(getListElement(R_pa, "lrate"));
	iter = 1;
	convergence = 0;
	while (!convergence && iter < 500) {

		/* We compute the new posterior variances and the gradient */

		computeNewPosteriorVariances(R_mR, R_vR, R_ret, R_aux, e, eP, eQ, vPbackup, vQbackup);
		computeGradient(R_mR, R_vR, R_ret, R_aux, e, eP, eQ, grP, grQ);

		/* We update the posterior means and then the variances */

		updatePosteriorMeans(R_ret, R_aux, vPbackup, vQbackup, grP, grQ, lrate);
		updatePosteriorVariances(R_ret, R_aux, vPbackup, vQbackup);

		/* We evaluate the lower bound */

		bound = evaluateLowerBound(R_mR, R_vR, R_ret, R_aux, e, eP, eQ, &avgTrE, &avgPV);

		if (bound > boundOld) {
			previousBounds[ (iter - 1) %  MAX_PREVIOUS_BOUNDS ] = bound;
			boundOld = bound;
			lrate = lrate * 1.1;

			/* We update the variances for the prior and for the level noise */

			updatePriorAndNoiseVariances(R_ret, R_aux);

		} else {
			previousBounds[ (iter - 1) %  MAX_PREVIOUS_BOUNDS ] = boundOld;
//			fprintf(stdout, "Slowing down the learning rate.\n");
//			fflush(stdout);
			lrate = lrate * 0.5;
			restorePosteriorMeansAndVariances(R_ret, R_aux, vPbackup, vQbackup, grP, grQ);
			bound = evaluateLowerBound(R_mR, R_vR, R_ret, R_aux, e, eP, eQ, &avgTrE, &avgPV);
		}
		
		/* We check for convergence */

		if (iter - 1 >= MAX_PREVIOUS_BOUNDS) {
			diffBounds = 0;
			for (i = 0 ; i < MAX_PREVIOUS_BOUNDS - 1 ; i++) {
				diffBounds += abs2(previousBounds[ i ] - previousBounds[ i + 1 ]);
			}
			diffBounds /= (MAX_PREVIOUS_BOUNDS - 1);
			diffBounds /= abs2(bound);

			if (diffBounds < 1e-6)
				convergence = 1;

//			fprintf(stdout, "%4d : %4f : %4f : %6f\n", iter, bound, sqrt(avgTrE), diffBounds); 
//			fflush(stdout);
		} else {
//			fprintf(stdout, "%4d : %4f : %4f\n", iter, bound, sqrt(avgTrE)); 
//			fflush(stdout);
		}

		iter++;
	}

	*NUMERIC_POINTER(getListElement(R_ret, "lrate")) = lrate;
	*NUMERIC_POINTER(getListElement(R_ret, "bound")) = bound;

	UNPROTECT(5);

	free(e); free(grP); free(grQ); free(eP); free(eQ); free(vPbackup); free(vQbackup);

	return R_ret;
}

/* Auxiliary function to get the list element named str, or return NULL */
     
SEXP getListElement(SEXP list, const char *str) {

	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

	int i;
     
	for (i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
	}

	return elmt;
}

/**
 * Function that evaluates the lower bound.
 */

double evaluateLowerBound(SEXP R_mR, SEXP R_vR, SEXP R_pa, SEXP R_aux, double *e, int *eP, int *eQ, double *avgTrE, double *avgPV) {

	double term1, term2, term3, *mPprior, *vPprior, *mP, *vP, *mQprior, *vQprior, *mQ, *vQ, *mR, *vR, p, v;
	int i, j, l, k, nP, nQ, nR, idxP, idxQ, *ia, *ja;
	
	k = *INTEGER_POINTER(getListElement(R_aux, "k"));
	nP = *INTEGER_POINTER(getListElement(R_aux, "nP"));
	nQ = *INTEGER_POINTER(getListElement(R_aux, "nQ"));

	/* We compute the value of the first term */

	mPprior = NUMERIC_POINTER(getListElement(R_pa, "mPprior"));
	vPprior = NUMERIC_POINTER(getListElement(R_pa, "vPprior"));
	mP = NUMERIC_POINTER(getListElement(R_pa, "mP"));
	vP = NUMERIC_POINTER(getListElement(R_pa, "vP"));
	
	term1 = 0;
	for (i = 0 ; i < nP * k; i++)
		term1 += 0.5 + 0.5 * log(vP[ i ] / vPprior[ i ]) - (mP[ i ] * mP[ i ] +
			vP[ i ] + mPprior[ i ] * mPprior[ i ] - 2 * mPprior[ i ] * mP[ i ]) / (2 * vPprior[ i ]);

	/* We compute the value of the second term */

	mQprior = NUMERIC_POINTER(getListElement(R_pa, "mQprior"));
	vQprior = NUMERIC_POINTER(getListElement(R_pa, "vQprior"));
	mQ = NUMERIC_POINTER(getListElement(R_pa, "mQ"));
	vQ = NUMERIC_POINTER(getListElement(R_pa, "vQ"));
	
	term2 = 0;
	for (i = 0 ; i < nQ * k; i++)
		term2 += 0.5 + 0.5 * log(vQ[ i ] / vQprior[ i ]) - (mQ[ i ] * mQ[ i ] +
			vQ[ i ] + mQprior[ i ] * mQprior[ i ] - 2 * mQprior[ i ] * mQ[ i ]) / (2 * vQprior[ i ]);

	/* We compute the value of the third term */

	ia = INTEGER_POINTER(getListElement(R_mR, "csc_ia"));
	ja = INTEGER_POINTER(getListElement(R_mR, "csc_ja"));
	mR = NUMERIC_POINTER(getListElement(R_mR, "csc_ra"));
	vR = NUMERIC_POINTER(getListElement(R_vR, "csc_ra"));

	/* We run through the movies */

	*avgTrE = 0;
	*avgPV = 0;
	term3 = 0;
	for (i = 0 ; i < nQ ; i++) {
		
		/* We run through the users that rated the movie */

		j = ia[ i ] - 1;
		while (j != ia[ i + 1 ] - 1) {

			/* We compute the prediction and the variance of the prediction and update the value of the third term */

			p = v = 0;
			for (l = 0 ; l < k ; l++) {
				idxP = ja[ j ] - 1 + l * nP;
				idxQ = i + l * nQ;
				p += mP[ idxP ] * mQ[ idxQ ];
				v += mP[ idxP ] * mP[ idxP ] * vQ[ idxQ ] + mQ[ idxQ ] * mQ[ idxQ ] * vP[ idxP ] + vQ[ idxQ ] * vP[ idxP ];
			}

			/* We store the error value and update the value of the third term */
	
			e[ j ] = mR[ j ] - p; eP[ j ] = ja[ j ] - 1; eQ[ j ] = i;
			term3 += (e[ j ] * e[ j ] + v) / (2 * vR[ j ]) + 0.5 * log(2 * PI * vR[ j ]);

			/* We update the average training error and the average predictive variance */
		
			*avgTrE += e[ j ] * e[ j ];
			*avgPV += v;

			j++;

		} /* while */

	} /* for */

	/* We normalize the sum of errors and predictive variances */

	nR = *INTEGER_POINTER(getListElement(R_aux, "nR"));

	*avgTrE /= nR;
	*avgPV /= nR;

	return term1 + term2 - term3;
}

/**
 * Function that computes the gradient of lower bound with respect to the posterior means.
 */

void computeGradient(SEXP R_mR, SEXP R_vR, SEXP R_pa, SEXP R_aux, double *e, int *eP, int *eQ, double *grP, double *grQ) {

	double *mPprior, *vPprior, *mQprior, *vQprior, *mP, *vP, *mQ, *vQ, *vR;
	int i, j, k, nP, nQ, nR;

	k  = *INTEGER_POINTER(getListElement(R_aux, "k"));
	nP = *INTEGER_POINTER(getListElement(R_aux, "nP"));
	nQ = *INTEGER_POINTER(getListElement(R_aux, "nQ"));
	nR = *INTEGER_POINTER(getListElement(R_aux, "nR"));

	mP = NUMERIC_POINTER(getListElement(R_pa, "mP"));
	vP = NUMERIC_POINTER(getListElement(R_pa, "vP"));
	mQ = NUMERIC_POINTER(getListElement(R_pa, "mQ"));
	vQ = NUMERIC_POINTER(getListElement(R_pa, "vQ"));

	/* We compute the gradient for the users */

	mPprior = NUMERIC_POINTER(getListElement(R_pa, "mPprior"));
	vPprior = NUMERIC_POINTER(getListElement(R_pa, "vPprior"));

	for (i = 0 ; i < nP * k ; i++)
		grP[ i ] = (mPprior[ i ] - mP[ i ]) / vPprior[ i ];

	vR = NUMERIC_POINTER(getListElement(R_vR, "csc_ra"));

	/* We run through the errors and the features */

	for (i = 0 ; i < nR ; i++) {
		for (j = 0 ; j < k ; j++) {
			grP[ eP[ i ] + j * nP ] += (e[ i ] * mQ[ eQ[ i ] + j * nQ ] - vQ[ eQ[ i ] + j * nQ ] * mP[ eP[ i ] + j * nP ]) / vR[ i ];
		}
	}

	/* We compute the gradient for the movies */

	mQprior = NUMERIC_POINTER(getListElement(R_pa, "mQprior"));
	vQprior = NUMERIC_POINTER(getListElement(R_pa, "vQprior"));

	for (i = 0 ; i < nQ * k ; i++)
		grQ[ i ] = (mQprior[ i ] - mQ[ i ]) / vQprior[ i ];

	/* We run through the errors and the features */

	for (i = 0 ; i < nR ; i++) {
		for (j = 0 ; j < k ; j++) {
			grQ[ eQ[ i ] + j * nQ ] += (e[ i ] * mP[ eP[ i ] + j * nP ] - vP[ eP[ i ] + j * nP ] * mQ[ eQ[ i ] + j * nQ ]) / vR[ i ];
		}
	}
}

/**
 * Function that computes the update for the posterior variances.
 */

void computeNewPosteriorVariances(SEXP R_mR, SEXP R_vR, SEXP R_pa, SEXP R_aux, double *e, int *eP, int *eQ, double *vPnew, double *vQnew) {

	double *vPprior, *vQprior, *mP, *vP, *mQ, *vQ, *vR;
	int i, j, k, nP, nQ, nR;

	k  = *INTEGER_POINTER(getListElement(R_aux, "k"));
	nP = *INTEGER_POINTER(getListElement(R_aux, "nP"));
	nQ = *INTEGER_POINTER(getListElement(R_aux, "nQ"));
	nR = *INTEGER_POINTER(getListElement(R_aux, "nR"));

	mP = NUMERIC_POINTER(getListElement(R_pa, "mP"));
	vP = NUMERIC_POINTER(getListElement(R_pa, "vP"));
	mQ = NUMERIC_POINTER(getListElement(R_pa, "mQ"));
	vQ = NUMERIC_POINTER(getListElement(R_pa, "vQ"));

	/* We compute the new variances for the users */

	vPprior = NUMERIC_POINTER(getListElement(R_pa, "vPprior"));

	for (i = 0 ; i < nP * k ; i++)
		vPnew[ i ] = 1 / vPprior[ i ];

	vR = NUMERIC_POINTER(getListElement(R_vR, "csc_ra"));

	/* We run through the errors and the features */

	for (i = 0 ; i < nR ; i++)
		for (j = 0 ; j < k ; j++)
			vPnew[ eP[ i ] + j * nP ] += (mQ[ eQ[ i ] + j * nQ ] * mQ[ eQ[ i ] + j * nQ ] + vQ[ eQ[ i ] + j * nQ ]) / vR[ i ];

	/* We compute the inverse of the variance */

	for (i = 0 ; i < nP * k ; i++)
		vPnew[ i ] = 1 / vPnew[ i ];

	/* We compute the new variances for the movies */

	vQprior = NUMERIC_POINTER(getListElement(R_pa, "vQprior"));

	for (i = 0 ; i < nQ * k ; i++)
		vQnew[ i ] = 1 / vQprior[ i ];

	/* We run through the errors and the features */

	for (i = 0 ; i < nR ; i++)
		for (j = 0 ; j < k ; j++)
			vQnew[ eQ[ i ] + j * nQ ] += (mP[ eP[ i ] + j * nP ] * mP[ eP[ i ] + j * nP ] + vP[ eP[ i ] + j * nP ]) / vR[ i ];

	/* We compute the inverse of the variance */

	for (i = 0 ; i < nQ * k ; i++)
		vQnew[ i ] = 1 / vQnew[ i ];

	/* The process is repeated because one update depends on the other */

	vPprior = NUMERIC_POINTER(getListElement(R_pa, "vPprior"));
	for (i = 0 ; i < nP * k ; i++)
		vPnew[ i ] = 1 / vPprior[ i ];
	for (i = 0 ; i < nR ; i++)
		for (j = 0 ; j < k ; j++)
			vPnew[ eP[ i ] + j * nP ] += (mQ[ eQ[ i ] + j * nQ ] * mQ[ eQ[ i ] + j * nQ ] + vQnew[ eQ[ i ] + j * nQ ]) / vR[ i ];
	for (i = 0 ; i < nP * k ; i++)
		vPnew[ i ] = 1 / vPnew[ i ];

	vQprior = NUMERIC_POINTER(getListElement(R_pa, "vQprior"));
	for (i = 0 ; i < nQ * k ; i++)
		vQnew[ i ] = 1 / vQprior[ i ];
	for (i = 0 ; i < nR ; i++)
		for (j = 0 ; j < k ; j++)
			vQnew[ eQ[ i ] + j * nQ ] += (mP[ eP[ i ] + j * nP ] * mP[ eP[ i ] + j * nP ] + vPnew[ eP[ i ] + j * nP ]) / vR[ i ];
	for (i = 0 ; i < nQ * k ; i++)
		vQnew[ i ] = 1 / vQnew[ i ];
}

/**
 * Function that updates the value of the posterior variances.
 */

void updatePosteriorVariances(SEXP R_ret, SEXP R_aux, double *vPnew, double *vQnew) {

	double *vP, *vQ, backup;
	int i, k, nP, nQ;

	k  = *INTEGER_POINTER(getListElement(R_aux, "k"));
	nP = *INTEGER_POINTER(getListElement(R_aux, "nP"));
	nQ = *INTEGER_POINTER(getListElement(R_aux, "nQ"));

	vP = NUMERIC_POINTER(getListElement(R_ret, "vP"));
	vQ = NUMERIC_POINTER(getListElement(R_ret, "vQ"));

	/* We update the posterior variances */

	for (i = 0 ; i < nP * k ; i++) {
		backup = vP[ i ];
		vP[ i ] = vPnew[ i ];
		vPnew[ i ] = backup;
	}
	for (i = 0 ; i < nQ * k ; i++) {
		backup = vQ[ i ];
		vQ[ i ] = vQnew[ i ];
		vQnew[ i ] = backup;
	}
}

/**
 * Function that updates the posterior means.
 */

void updatePosteriorMeans(SEXP R_ret, SEXP R_aux, double *vPnew, double *vQnew, double *grP, double *grQ, double lrate) {

	double *mP, *mQ, backup, alpha = 2.0 / 3;
	int i, k, nP, nQ;

	k  = *INTEGER_POINTER(getListElement(R_aux, "k"));
	nP = *INTEGER_POINTER(getListElement(R_aux, "nP"));
	nQ = *INTEGER_POINTER(getListElement(R_aux, "nQ"));

	mP = NUMERIC_POINTER(getListElement(R_ret, "mP"));
	mQ = NUMERIC_POINTER(getListElement(R_ret, "mQ"));

	/* We update the posterior means */

	for (i = 0 ; i < nP * k ; i++) {
		backup = mP[ i ];
		mP[ i ] += lrate * pow(vPnew[ i ], alpha) * grP[ i ];
		grP[ i ] = backup;
	}

	for (i = 0 ; i < nQ * k ; i++) {
		backup = mQ[ i ];
		mQ[ i ] += lrate * pow(vQnew[ i ], alpha) * grQ[ i ];
		grQ[ i ] = backup;
	}
}

/**
 * Function that restores the psterior means and variances.
 */

void restorePosteriorMeansAndVariances(SEXP R_ret, SEXP R_aux, double *vPold, double *vQold, double *grP, double *grQ) {

	double *mP, *mQ, *vP, *vQ;
	int i, k, nP, nQ;

	k  = *INTEGER_POINTER(getListElement(R_aux, "k"));
	nP = *INTEGER_POINTER(getListElement(R_aux, "nP"));
	nQ = *INTEGER_POINTER(getListElement(R_aux, "nQ"));

	mP = NUMERIC_POINTER(getListElement(R_ret, "mP"));
	mQ = NUMERIC_POINTER(getListElement(R_ret, "mQ"));
	vP = NUMERIC_POINTER(getListElement(R_ret, "vP"));
	vQ = NUMERIC_POINTER(getListElement(R_ret, "vQ"));

	/* We restore the posterior means and variances */

	for (i = 0 ; i < nP * k ; i++) {
		mP[ i ] = grP[ i ];
		vP[ i ] = vPold[ i ];
	}

	for (i = 0 ; i < nQ * k ; i++) {
		mQ[ i ] = grQ[ i ];
		vQ[ i ] = vQold[ i ];
	}
}

/**
 * Function updates the variances of the prior and the level of noise.
 */

void updatePriorAndNoiseVariances(SEXP R_ret, SEXP R_aux) {

	double *mP, *mQ, *vP, *vQ, *mPprior, *mQprior, *vPprior, *vQprior, newVariance;
	int i, j, k, nP, nQ, optimizePrior;

	mP = NUMERIC_POINTER(getListElement(R_ret, "mP"));
	mQ = NUMERIC_POINTER(getListElement(R_ret, "mQ"));
	vP = NUMERIC_POINTER(getListElement(R_ret, "vP"));
	vQ = NUMERIC_POINTER(getListElement(R_ret, "vQ"));
	mPprior = NUMERIC_POINTER(getListElement(R_ret, "mPprior"));
	mQprior = NUMERIC_POINTER(getListElement(R_ret, "mQprior"));
	vPprior = NUMERIC_POINTER(getListElement(R_ret, "vPprior"));
	vQprior = NUMERIC_POINTER(getListElement(R_ret, "vQprior"));

	nP = *INTEGER_POINTER(getListElement(R_aux, "nP"));
	nQ = *INTEGER_POINTER(getListElement(R_aux, "nQ"));
	k  = *INTEGER_POINTER(getListElement(R_aux, "k"));

	optimizePrior  = *INTEGER_POINTER(getListElement(R_aux, "optimizePrior"));
	
	if (optimizePrior) {

		/* We update the variance of the prior for the coefficients */
		
		for (i = 0 ; i < k ; i++) {
	
			newVariance = 0;
			for (j = 0 ; j < nQ ; j++)
				newVariance += vQ[ j + i * nQ ] + (mQ[ j + i * nQ ] - mQprior[ j + i * nQ ]) * (mQ[ j + i * nQ ] - mQprior[ j + i * nQ ]);
			newVariance /= nQ;

			for (j = 0 ; j < nQ ; j++)
				vQprior[ j + i * nQ ] = newVariance;
		}
	}
}
