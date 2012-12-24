#
# R script that extract the features and creates the train test sets
#

set.seed(1)

library(impute)

library(Matrix)

# We delete the existing results

system("rm -r *.txt")

# We load the rating matrix

ratings <- as.matrix(read.table("ratingsTranslated.dat"))
ratingMatrix <- sparseMatrix(i = ratings[ , 1 ], j = ratings[ , 2 ])

# We randomly sample 10 movies among the top 50 most rated movies

mostRatedMovies <- sort(apply(ratingMatrix, 2, sum), decreasing = T, index.return = T)$ix[ 1 : 50 ]
mostRatedMovies <- sample(mostRatedMovies, 10)

write.table(10, "nItems.txt", col.names = F, row.names = F)

index <- which(sapply(ratings[ , 2 ], function(x) is.element(x, mostRatedMovies)))
ratings <- ratings[ index, ]
aux <- seq(1, ncol(ratingMatrix))
newIdsMostRatedMovies <- sort(mostRatedMovies, index.return = T)$ix
aux[ mostRatedMovies ] <- newIdsMostRatedMovies
ratings[ , 2 ] <- aux[ ratings[ , 2 ] ]

# We select those users with 7 or more ratings

auxMatrix <- sparseMatrix(i = ratings[ , 1 ], j = ratings[ , 2 ], x = 1)
indexValidUsers <- which(apply(auxMatrix, 1, sum) > 6)
auxMatrix <- matrix(NA, max(ratings[ , 1 ]), max(ratings[ , 2 ]))
auxMatrix[ ratings[ , c(1, 2) ] ] <- ratings[ , 3 ]
auxMatrix <- auxMatrix[ indexValidUsers, ]

# We impute missing ratings

auxMatrix <- impute.knn(auxMatrix)$data
ratingMatrixFinal <- auxMatrix
ratingMatrixFinal <- ratingMatrixFinal +
	matrix(rnorm(nrow(ratingMatrixFinal) * ncol(ratingMatrixFinal)), nrow(ratingMatrixFinal), ncol(ratingMatrixFinal)) * 1e-5

# We obtain the item features

movieFeatures <- as.matrix(read.table("finalMoviesTableTranslated.dat"))[ , -1 ]
movieFeatures <- movieFeatures[ mostRatedMovies, ]

write.table(movieFeatures, "itemFeatures.txt", col.names = F, row.names = F)

# We load the features of the different users

userFeatures <- as.matrix(read.table("finalUsersTable.dat"))[ , -1 ]

# We load the preferences of the different users

userPreferences <- t(apply(ratingMatrixFinal, 1, rank))

# We randomly pick 5000 users which have not been used for tuning the kernel parameters

candidates <- seq(1, nrow(userPreferences))

nUsers <- 1000
index <- sample(candidates, nUsers)
userPreferencesSplit <- userPreferences[ index, ]
userFeaturesFinal <- userFeatures[ indexValidUsers[ index ], ]

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
