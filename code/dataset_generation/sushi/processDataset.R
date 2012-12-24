#
# R script that extract the sushi features and creates the train test sets
#

# We delete the existing results

system("rm -r *.txt")

# We load the sushi item descriptions

sushiData <- read.table("sushi3.idata")

# The ids of the items we are interested in

itemIds <- c(0, 1, 2, 3, 4, 15, 7, 8, 26, 29)

# We store the number of items

write.table(length(itemIds), "nItems.txt", col.names = F, row.names = F)

# We select those items we are interested in

sushiData <- sushiData[ itemIds + 1, ]

# We replace the categorical variables for the third feature

categories <- unique(sushiData[ , 5 ])
for (i in 1 : nrow(sushiData)) {
	sushiData[ i , 5 ] <- which(categories == sushiData[ i, 5 ])
}

# We create the feature matrix for each item

featureMatrix <- matrix(0, 10, 15)
for (i in 1 : nrow(sushiData)) {
	if (sushiData[ i, 3 ] == 1)
		featureMatrix[ i, 1 ] <- 1
	else
		featureMatrix[ i, 2 ] <- 1

	if (sushiData[ i, 4 ] == 1)
		featureMatrix[ i, 3 ] <- 1
	else
		featureMatrix[ i, 4 ] <- 1

	featureMatrix[ i, 4 + sushiData[ i, 5 ] ] <- 1

	for (j in 1 : 4) 
		featureMatrix[ i, 4 + length(categories) + j ] <- sushiData[ i, 5 + j ]
}

# We store the faeture matrix

write.table(featureMatrix, "itemFeatures.txt", col.names = F, row.names = F)

# We load the features of the different users

userFeatures <- as.matrix(read.table("sushi3.udata"))[ , -c(1, 4, 5, 6, 7, 11) ]

# Indicator function that maps a categorical value to its corresponding binary vector

indicator <- function(x, nValues) { ret <- rep(0, nValues); ret[ x ] <- 1 ; ret }

# We map each categorical variable to its binary representation

userFeaturesFinal <- matrix(0, nrow(userFeatures), 0)
for (i in 1 : ncol(userFeatures)) {

	# We map the month to its category

	userFeatures[ , i ] <- as.integer(as.factor(userFeatures[ , i ]))
	nValues <- length(unique(userFeatures[ , i ]))
	if (nValues > 2)
		binarizedFeature <- t(apply(matrix(userFeatures[ , i ], nrow(userFeatures), 1), 1, function(x) indicator(x, nValues)))
	else
		binarizedFeature <- as.integer(userFeatures[ , i ] == userFeatures[ 1, i ])

	userFeaturesFinal <- cbind(userFeaturesFinal, as.matrix(binarizedFeature))

	print(i)
}

# We load the preferences of the different users

userPreferences <- read.table("sushi3a.5000.10.order")[ , 3 : 12 ] + 1

# We randomly pick 5000 users which have not been used for tuning the kernel parameters

candidates <- seq(1, nrow(userPreferences))

nUsers <- 100
index <- sample(candidates, nUsers)
userPreferencesSplit <- userPreferences[ index, ]
userFeaturesFinal <- userFeaturesFinal[ index, ]

# We store the user features

write.table(userFeaturesFinal, "userFeatures.txt", col.names = F, row.names = F)

# We create all the data for the users

dataset <- c()
for (i in 1 : nrow(userPreferencesSplit)) {

	auxMatrix <- matrix(0, 10 * 9 / 2, 4)
	auxMatrix[ , 3 ] <- i

	counter <- 1
	for (j in 1 : 9) {
		for (k in (j + 1) : 10) {
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

		# The training set will have 5 elements

		indexTrain <- sample(index, 5)

		# The pool set will have 35 elements

		candidatesPool <- setdiff(index, indexTrain)
		indexPool <- sample(candidatesPool, 35)

		# The test set will have 5 elements

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
