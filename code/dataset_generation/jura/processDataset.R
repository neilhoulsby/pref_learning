#
# R script that extract the features and creates the train test sets
#

# Indicator function that maps a categorical value to its corresponding binary vector

indicator <- function(x, nValues) { ret <- rep(0, nValues); ret[ x ] <- 1 ; ret }

# We delete the existing results

system("rm -r *.txt")

write.table(7, "nItems.txt", col.names = F, row.names = F)

# We load the dataset

juraData <- as.matrix(read.table("jura.dat"))

# We normalize the heavy metal measurements

measurements <- juraData[ , 5 : ncol(juraData) ]
meanAux <- apply(measurements, 2, mean)
sdAux <- apply(measurements, 2, sd)
measurements <- (measurements - matrix(meanAux, nrow(measurements), ncol(measurements), byrow = T)) /
	matrix(sdAux, nrow(measurements), ncol(measurements), byrow = T)

# We create the feature matrix for each user

userFeatures <- juraData[ , c(1, 2) ]

# We binarized the categorical features

juraData[ , 3 ] <- as.integer(as.factor(juraData[ , 3 ]))
nValues <- length(unique(juraData[ , 3 ]))
if (nValues > 2) {
	binarizedFeature <- t(apply(matrix(juraData[ , 3 ], nrow(juraData), 1), 1, function(x) indicator(x, nValues)))
} else {
	binarizedFeature <- as.integer(juraData[ , 3 ] == juraData[ 1, 3 ])
}

userFeatures <- cbind(userFeatures, as.matrix(binarizedFeature))

juraData[ , 4 ] <- as.integer(as.factor(juraData[ , 4 ]))
nValues <- length(unique(juraData[ , 4 ]))
if (nValues > 2) {
	binarizedFeature <- t(apply(matrix(juraData[ , 4 ], nrow(juraData), 1), 1, function(x) indicator(x, nValues)))
} else {
	binarizedFeature <- as.integer(juraData[ , 4 ] == juraData[ 1, 4 ])
}

userFeatures <- cbind(userFeatures, as.matrix(binarizedFeature))

# We create the item features

indexItemFeatures <- sample(1 : nrow(juraData), 20)
featureMatrix <- t(measurements[ indexItemFeatures, ])

# We store the faeture matrix

write.table(featureMatrix, "itemFeatures.txt", col.names = F, row.names = F)

# We load the preferences of the different users

userPreferences <- measurements

# We randomly pick 349 users which have not been used for tuning the kernel parameters

candidates <- seq(1, nrow(userPreferences))[ -indexItemFeatures ]

nUsers <- 339
index <- sample(candidates, nUsers)
userPreferencesSplit <- t(apply(userPreferences[ index, ], 1, rank))
userFeaturesFinal <- userFeatures[ index, ]

# We store the user features

write.table(userFeaturesFinal, "userFeatures.txt", col.names = F, row.names = F)

# We create all the data for the users

dataset <- c()
for (i in 1 : nrow(userPreferencesSplit)) {

	auxMatrix <- matrix(0, 7 * 6 / 2, 4)
	auxMatrix[ , 3 ] <- i

	counter <- 1
	for (j in 1 : 6) {
		for (k in (j + 1) : 7) {
			auxMatrix[ counter, 1 ] <- j
			auxMatrix[ counter, 2 ] <- k

			if (which(userPreferencesSplit[ i, ] == j) < which(userPreferencesSplit[ i, ] == k))
				auxMatrix[ counter, 4 ] <- 1
			else
				auxMatrix[ counter, 4 ] <- -1

			counter <- counter + 1
		}
	}

	dataset <- rbind(dataset, auxMatrix)
}

itemDim <- 1

# We split the data in training, test and pool set

for (i in 1 : 25) {

	trainSet <- matrix(0, 0, 2 * itemDim + 2)
	testSet <- matrix(0, 0, 2 * itemDim + 2)
	poolSet <- matrix(0, 0, 2 * itemDim + 2)
	for (j in 1 : nUsers) {

		index <- which(dataset[ , 2 * itemDim + 1 ] == j)

		# The training set will have 3 elements

		indexTrain <- sample(index, 3)

		# The pool set will have 10 elements

		candidatesPool <- setdiff(index, indexTrain)
		indexPool <- sample(candidatesPool, 15)

		# The test set will have 2 elements

		candidatesTest <- setdiff(index, c(indexTrain, indexPool))
		indexTest <- sample(candidatesTest, 3)

		trainSet <- rbind(trainSet, dataset[ indexTrain, ])
		poolSet <- rbind(poolSet, dataset[ indexPool, ])
		testSet <- rbind(testSet, dataset[ indexTest, ])
	}

	write.table(trainSet, paste("train", i, ".txt", sep = ""), col.names = F, row.names = F)
	write.table(poolSet, paste("pool", i, ".txt", sep = ""), col.names = F, row.names = F)
	write.table(testSet, paste("test", i, ".txt", sep = ""), col.names = F, row.names = F)

	print(i)
}
