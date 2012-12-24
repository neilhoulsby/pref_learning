#
# R script that extract the election features and creates the train test sets
#

library(impute)

set.seed(1)

# We delete the existing results

system("rm -r *.txt")

# We load the data

dataElection <- as.matrix(read.table("all_polling_data.dat", sep = "\t"))[ , -1 ]
userFeatures <- dataElection[ , c(1 : 2) ]
auxMatrix <- dataElection[ , 3 : ncol(dataElection) ]

# We only keep those constituencies with at least 6 non-missing entries

index <- which(apply(!is.na(auxMatrix), 1, sum) >= 6)
userFeatures <- userFeatures[ index, ]
auxMatrix <- auxMatrix[ index, ]

# We impute missing ratings

auxMatrix <- impute.knn(auxMatrix, colmax = 0.99)$data
ratingMatrixFinal <- auxMatrix
ratingMatrixFinal <- ratingMatrixFinal +
	matrix(rnorm(nrow(ratingMatrixFinal) * ncol(ratingMatrixFinal)), nrow(ratingMatrixFinal), ncol(ratingMatrixFinal)) * 1e-5

# We obtain the item features

indexItemFeatures <- sample(1 : nrow(ratingMatrixFinal), 20)
itemFeatures <- t(ratingMatrixFinal[ indexItemFeatures, ])
ratingMatrixFinal <- ratingMatrixFinal[ -indexItemFeatures, ]
userFeatures <- userFeatures[ -indexItemFeatures, ]

write.table(itemFeatures, "itemFeatures.txt", col.names = F, row.names = F)
write.table(nrow(itemFeatures), "nItems.txt", col.names = F, row.names = F)

# We load the features of the different users

userPreferences <- t(apply(ratingMatrixFinal, 1, rank))

# We randomly pick 5000 users which have not been used for tuning the kernel parameters

candidates <- seq(1, nrow(userPreferences))

nUsers <- nrow(userPreferences)
index <- sample(candidates, nUsers)
userPreferencesSplit <- userPreferences[ index, ]
userFeaturesFinal <- userFeatures[ index, ]

# We store the user features

write.table(userFeaturesFinal, "userFeatures.txt", col.names = F, row.names = F)

# We create all the data for the users

dataset <- c()
for (i in 1 : nrow(userPreferencesSplit)) {

	auxMatrix <- matrix(0, 8 * 7 / 2, 4)
	auxMatrix[ , 3 ] <- i

	counter <- 1
	for (j in 1 : 7) {
		for (k in (j + 1) : 8) {
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

		indexTrain <- sample(index, 3)

		# The pool set will have 35 elements

		candidatesPool <- setdiff(index, indexTrain)
		indexPool <- sample(candidatesPool, 22)

		# The test set will have 5 elements

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
