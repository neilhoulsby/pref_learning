#
# R script that evaluates the performance of the multitask preference learning model on synthetic data.
#
# Author: Jose Miguel Hernandez Lobato
#
# Date: 24 Dec 2011
#

set.seed(1)

source("epGPC.R")

# Choose active learning method from {bald, entropy, random}

active_method <- "bald"

# Choose datataset

datapath <- "../../../dataset_generation/sushi/"

# We eliminate the previous results

system("rm -f results/*/*.txt")

# We load the simulation number

simNumber <- read.table("simulationNumber.txt")$V1

# We repeat for the current train/test partition

for (i in simNumber) {

	# We load the item descriptions

	itemFeatures <- as.matrix(read.table(paste(datapath, "/itemFeatures.txt", sep="")))

	# We load the training, pool and test sets

	train <- as.matrix(read.table(paste(datapath, "/train", i, ".txt", sep = "")))
	test <- as.matrix(read.table(paste(datapath, "/test", i, ".txt", sep = "")))
	pool <- as.matrix(read.table(paste(datapath, "/pool", i, ".txt", sep = "")))

	# We do nQueries iterations of active learning

	nQueries <- 30
	for (j in 1 : (nQueries + 1)) {

		# We standardize the item features so that they have zero mean and unit standard deviation in the training set

		meanItemFeatures <- apply(itemFeatures[ unique(c(train[ , 1 ], train[ , 2 ])), ], 2, mean)
		sdItemFeatures <- apply(itemFeatures[ unique(c(train[ , 1 ], train[ , 2 ])), ], 2, sd)

		itemFeaturesStandardized <- (itemFeatures - matrix(meanItemFeatures, nrow(itemFeatures), ncol(itemFeatures), byrow = T)) /
			matrix(sdItemFeatures, nrow(itemFeatures), ncol(itemFeatures), byrow = T)

		# We create the training data to pass to the EP method

		itemsAtrain <- itemFeaturesStandardized[ train[ , 1 ], ]
		itemsBtrain <-itemFeaturesStandardized[ train[ , 2 ], ]
		userIdsTrain <- train[ , 3 ]
		ratingsTrain <- train[ , 4 ]

		# We obtain the length scale parameter

		itemFeaturesAux <- itemFeaturesStandardized[ unique(c(train[ , 1 ], train[ , 2 ])), ]
		n <- nrow(itemFeaturesAux)
		Q <- matrix(apply(itemFeaturesAux^2, 1, sum), n, n)
		distance <- Q + t(Q) - 2 * itemFeaturesAux %*% t(itemFeaturesAux)
		lengthScale <- sqrt(0.5 * (median(distance[ upper.tri(distance) ])))

		# We fit a GPC for each userId

		gpc <- list()
		totalTime <- 0
		for (k in unique(userIdsTrain)) {

			index <- which(userIdsTrain == k)

			itemsAlocal <- matrix(itemsAtrain[ index, ], length(index), ncol(itemFeatures))
			itemsBlocal <- matrix(itemsBtrain[ index, ], length(index), ncol(itemFeatures))
			userIdsLocal <- userIdsTrain[ index ]
			ratingsLocal <- ratingsTrain[ index ]

			time <- system.time(gpc[[ k ]] <- epGPCinternal(itemsAlocal, itemsBlocal, userIdsLocal, ratingsLocal, lengthScale))
			totalTime <- totalTime + time[[ 1 ]] + time[[ 2 ]]
		}
		write.table(totalTime, paste("results/", active_method, "/time", j - 1, ".txt", sep = ""), col.names = F, row.names = F, append = T)

		# We eavaluate the performance of the method on the test dataset

		itemsAtest <- itemFeaturesStandardized[ test[ , 1 ], ]
		itemsBtest <- itemFeaturesStandardized[ test[ , 2 ], ]
		userIdsTest <- test[ , 3 ]
		ratingsTest <- test[ , 4 ]

		pred <- rep(0, nrow(test))
		for (k in unique(userIdsTest)) {
			index <- which(userIdsTest == k)
			pred[ index ] <- sign(predictGPCvectorized(gpc[[ k ]], itemsAtest[ index, ], itemsBtest[ index, ])$m)
		}

		error <- mean(pred != ratingsTest)

		write.table(error, paste("results/", active_method, "/error", j - 1, ".txt", sep = ""), col.names = F, row.names = F, append = T)

		# We make predictions in the pool set

		itemsApool <- itemFeaturesStandardized[ pool[ , 1 ], ]
		itemsBpool <- itemFeaturesStandardized[ pool[ , 2 ], ]
		userIdsPool <- pool[ , 3 ]

		toAdd <- c()
		toRemove <- c()
		for (k in unique(userIdsPool)) {
			index <- which(userIdsPool == k)
			pred <- predictGPCvectorized(gpc[[ k ]], itemsApool[ index, ], itemsBpool[ index, ])

			# We evaluate the bald score for each point

			prob <- pnorm(pred$m / sqrt(pred$v + 1))
			firstTerm <- -(1 - prob) * log(1 - prob) - prob * log(prob)
			C <- sqrt(pi * log(2)  / 2)
			secondTerm <- C / (sqrt(pred$v) + C^2) * exp(-pred$m^2 / (2 * (pred$v + C^2)))

      if (active_method == "bald") {
        al <- firstTerm - secondTerm
      } else if (active_method == "entropy") {
        al <- firstTerm
      } else if (active_methods == "random") {
        al <- rnorm(length(firstTerm))
      } else {
        print("Invalid active learning scheme.")
        stop()
      }

			# We select the most informative point

			toAdd <- rbind(toAdd, pool[ index[ which.max(al) ], ])
			toRemove <- c(toRemove, index[ which.max(al) ])
		}

		# We include these items in the train dataset and eliminate them from the pool set

		train <- rbind(train, toAdd)
		pool <- pool[ -toRemove, ]

		print(j)
	}
}
