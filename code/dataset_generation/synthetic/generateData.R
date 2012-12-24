#
# R script that generates synthetic data from a multitask preference learning model
#
# Author: Jose Miguel Hernandez Lobato
#
# Date: 23 Dec 2011
#

system("rm -f *.txt")

set.seed(1)

# Function that computes the kernel matrix for the user weights

computeKernel <- function(X, lengthScale) {

	n <- nrow(X)

	# We compute the kernel matrix

	Q <- matrix(apply(X^2, 1, sum), n, n)
	distance <- Q + t(Q) - 2 * X %*% t(X)
	K <- exp(-0.5 * distance / lengthScale^2)

	K + diag(rep(1e-5), n)
}

# Function that computes the preference kernel matrix

computePreferenceKernel <- function(itemsA, itemsB, lengthScale) {

	n <- nrow(itemsA)

	# We compute the kernel matrix

	QA <- matrix(apply(itemsA^2, 1, sum), n, n)
	distanceA <- QA + t(QA) - 2 * itemsA %*% t(itemsA)
	KA <- exp(-0.5 * distanceA / lengthScale^2)
	KA <- KA

	QB <- matrix(apply(itemsB^2, 1, sum), n, n)
	distanceB <- QB + t(QB) - 2 * itemsB %*% t(itemsB)
	KB <- exp(-0.5 * distanceB / lengthScale^2)
	KB <- KB

	distanceAB <- QA + t(QB) - 2 * itemsA %*% t(itemsB)
	KAB <- exp(-0.5 * distanceAB / lengthScale^2)
	KAB <- KAB

	distanceBA <- QB + t(QA) - 2 * itemsB %*% t(itemsA)
	KBA <- exp(-0.5 * distanceBA / lengthScale^2)
	KBA <- KBA

	K <- KA + KB - KAB - KBA

	K + diag(rep(1e-5), n)
}

# Number of items, dimension of item features and dimension of user features

nItems <- 10
itemDim <- 2
userDim <- 2

# We generate the item features from a uniform distribution with zero mean and unit standard deviation

itemFeatures <- matrix((runif(itemDim * nItems) - 0.5) / sqrt(1 / 12), nItems, itemDim)

# We store the item features

write.table(itemFeatures, "itemFeatures.txt", col.names = F, row.names = F)

# We store the number of items

write.table(nItems, "nItems.txt", col.names = F, row.names = F)

# We generate all the item pairs

allItemPairs <- matrix(0, nItems * (nItems -1) / 2, 2 * itemDim)
allItemIds <- matrix(0, nItems * (nItems -1) / 2, 2)
counter <- 1
for (i in 1 : (nItems - 1)) {
	for (j in (i + 1) : nItems) {
		allItemPairs[ counter, ] <- c(itemFeatures[ i, ], itemFeatures[ j, ])
		allItemIds[ counter, ] <- c(i, j)
		counter <- counter + 1
	}
}

# Number of users for which preferences will be obtained.

nUsers <- 1000

# We generate the data for each user.

dataset <- c()
for (i in 1 : nUsers) {
	dataset <- rbind(dataset, cbind(allItemPairs, rep(i, nrow(allItemPairs)), rep(0, nrow(allItemPairs))))

	print(i)
}

# Number of latent Gaussian process factors.

d <- 5

# We generate the user features from a uniform distribution with zero mean and unit standard deviation

userFeatures <- matrix((runif(userDim * nUsers) - 0.5) / sqrt(1 / 12), nUsers, userDim)

# We store the user features

write.table(userFeatures, "userFeatures.txt", col.names = F, row.names = F)

# The kernel lengthscale for the covariance functions

lengthScale <- 1

# We generate the kernel for the users

K <- computeKernel(userFeatures, lengthScale)

# We sample the mixing coefficients for the different users

cholK <- chol(K)
w <- matrix(0, nUsers, d)
for (i in 1 : d) {
	w[ , i ] <- t(rnorm(nrow(cholK))) %*% cholK
}

# We construct the preference kernels

K <- computePreferenceKernel(allItemPairs[ , 1 : itemDim ], allItemPairs[ , (itemDim + 1) : (2 * itemDim) ], lengthScale)

# We sample the values of the Gaussian processes h

cholK <- chol(K)
h <- matrix(0, d, nrow(cholK))
for (i in 1 : d) {
	h[ i, ] <- t(rnorm(nrow(cholK))) %*% cholK
}

# We generate the observations for each user

for (i in 1 : nUsers) {
	currentRatings <- which(dataset[ , 2 * itemDim + 1 ] == i)
	dataset[ currentRatings, 2 * itemDim + 2 ] <- sign(matrix(w[ i, ], 1, d) %*% h + rnorm(ncol(h)))
}

# We save the data

write.table(dataset, "dataset.txt", col.names = F, row.names = F)

# We save the matrix W used to generate the data

write.table(w, "w.txt", col.names = F, row.names = F)

# We save the matrix h used to generate the data

write.table(h, "h.txt", col.names = F, row.names = F)

# We map the item features to item ids

dataset <- cbind(rep(allItemIds[ , 1 ], nUsers), rep(allItemIds[ ,2 ], nUsers), dataset[ , 2 * itemDim + 1 ], dataset[ , 2 * itemDim + 2 ])
itemDim <- 1

# We split the dataset in training and test pairs.
# Each training set contains 2 / 3 of the ratings by each user.

for (i in 1 : 25) {

	trainSet <- matrix(0, 0, 2 * itemDim + 2)
	testSet <- matrix(0, 0, 2 * itemDim + 2)
	poolSet <- matrix(0, 0, 2 * itemDim + 2)
	for (j in 1 : nUsers) {

		index <- which(dataset[ , 2 * itemDim + 1 ] == j)

		# The training set will have 10 elements

		indexTrain <- sample(index, 5)

		# The pool set will have 20 elements

		candidatesPool <- setdiff(index, indexTrain)
		indexPool <- sample(candidatesPool, 35)

		# The test set will have 15 elements

		candidatesTest <- setdiff(index, c(indexTrain, indexPool))
		indexTest <- sample(candidatesTest, 5)

		trainSet <- rbind(trainSet, dataset[ indexTrain, ])
		poolSet <- rbind(poolSet, dataset[ indexPool, ])
		testSet <- rbind(testSet, dataset[ indexTest, ])
	}

	write.table(trainSet, paste("train", i, ".txt", sep = ""), col.names = F, row.names = F)
	write.table(poolSet, paste("pool", i, ".txt", sep = ""), col.names = F, row.names = F)
	write.table(testSet, paste("test", i, ".txt", sep = ""), col.names = F, row.names = F)

	print(i)
}
