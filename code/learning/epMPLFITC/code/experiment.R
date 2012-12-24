#
# R script that evaluates the performance of the multitask preference learning model on synthetic data.
#
# Author: Jose Miguel Hernandez Lobato
#
# Date: 24 Dec 2011
#

set.seed(1)

source("epMPL.R")

# Choose active learning method from {bald, entropy, random}

active_method <- "bald"

# Choose datataset

datapath <- "../../../dataset_generation/synthetic/"

# We eliminate the previous results

system("rm -f results/*/*.txt")

# We load the item descriptions

itemFeatures <- as.matrix(read.table(paste(datapath, "/itemFeatures.txt", sep="")))

# We load the user descriptions

userFeatures <- as.matrix(read.table(paste(datapath, "/userFeatures.txt", sep="")))

# We standardize the user features so that they have zero mean and unit standard deviation

meanUserFeatures <- apply(userFeatures, 2, mean)
sdUserFeatures <- apply(userFeatures, 2, sd)
sdUserFeatures[ sdUserFeatures == 0 ] <- 1

userFeatures <- (userFeatures - matrix(meanUserFeatures, nrow(userFeatures), ncol(userFeatures), byrow = T)) /
	matrix(sdUserFeatures, nrow(userFeatures), ncol(userFeatures), byrow = T)

# We load the number of different items

nItems <- read.table(paste(datapath, "/nItems.txt", sep=""))$V1

# We load the simulation number

simNumber <- read.table("simulationNumber.txt")$V1

# We repeat for the current train/test partition

for (i in simNumber) {

	# We load the training, pool and test sets

	train <- as.matrix(read.table(paste(datapath, "/train", i, ".txt", sep = "")))
	test <- as.matrix(read.table(paste(datapath, "/test", i, ".txt", sep = "")))
	pool <- as.matrix(read.table(paste(datapath, "/pool", i, ".txt", sep = "")))

	# We do nQueries iterations of active learning

	nQueries <- 1
	for (j in 1 : (nQueries + 1)) {

		# We standardize the item features so that they have zero mean and unit standard deviation in the training set

		meanItemFeatures <- apply(itemFeatures[ unique(c(train[ , 1 ], train[ , 2 ])), ], 2, mean)
		sdItemFeatures <- apply(itemFeatures[ unique(c(train[ , 1 ], train[ , 2 ])), ], 2, sd)

		itemFeaturesStandardized <- (itemFeatures - matrix(meanItemFeatures, nrow(itemFeatures), ncol(itemFeatures), byrow = T)) /
			matrix(sdItemFeatures, nrow(itemFeatures), ncol(itemFeatures), byrow = T)

		# We create the training data to pass to the EP method

		itemFeaturesAtrain <- itemFeaturesStandardized[ train[ , 1 ], ]
		itemFeaturesBtrain <-itemFeaturesStandardized[ train[ , 2 ], ]
		userIdsTrain <- train[ , 3 ]
		ratingsTrain <- train[ , 4 ]

		# We select the kernel width

		itemFeaturesAux <- itemFeaturesStandardized[ unique(c(train[ , 1 ], train[ , 2 ])), ]
		n <- nrow(itemFeaturesAux)
		Q <- matrix(apply(itemFeaturesAux^2, 1, sum), n, n)
		distance <- Q + t(Q) - 2 * itemFeaturesAux %*% t(itemFeaturesAux)
		lengthScaleItems <- sqrt(0.5 * (median(distance[ upper.tri(distance) ])))

		itemFeaturesAux <- userFeatures
		n <- nrow(itemFeaturesAux)
		Q <- matrix(apply(itemFeaturesAux^2, 1, sum), n, n)
		distance <- Q + t(Q) - 2 * itemFeaturesAux %*% t(itemFeaturesAux)
		lengthScaleUsers <- sqrt(0.5 * (median(distance[ upper.tri(distance) ])))

		# We construct the problemInfo list

		problemInfo <- list(itemFeaturesA = itemFeaturesAtrain, itemFeaturesB = itemFeaturesBtrain,
			itemIdsA = train[ , 1 ], itemIdsB = train[ , 2 ],
			userIds = userIdsTrain, ratings = ratingsTrain, d = 20,
			nTotalItems = nItems, userFeatures = userFeatures, nPseudoInputUsers = 25, nPseudoInputItems = 25)

		# We call the EP method

		time <- system.time(ret <- epMPL(problemInfo, lengthScaleItems, lengthScaleUsers))
		write.table(time[[ 1 ]] + time[[ 2 ]], paste("results/", active_method, "/time", j - 1, ".txt", sep = ""), col.names = F, row.names = F, append = T)

		# We eavaluate the performance of the method on the test set

		userIdsTest <- test[ , 3 ]
		ratingsTest <- test[ , 4 ]
		itemFeaturesAtest <- itemFeaturesStandardized[ test[ , 1 ], ]
		itemFeaturesBtest <- itemFeaturesStandardized[ test[ , 2 ], ]

		pred <- sign(predictMPL(ret, userIdsTest, itemFeaturesAtest, itemFeaturesBtest)$m)
		error <- mean(pred != ratingsTest)
		write.table(error, paste("results/", active_method, "/error", j - 1, ".txt", sep = ""), col.names = F, row.names = F, append = T)

		# We make predictions in the pool set

		userIdsPool <- pool[ , 3 ]
		itemFeaturesApool <- itemFeaturesStandardized[ pool[ , 1 ], ]
		itemFeaturesBpool <- itemFeaturesStandardized[ pool[ , 2 ], ]

		pred <- predictMPL(ret, userIdsPool, itemFeaturesApool, itemFeaturesBpool)

		# We compute the bald score for each point in the pool set

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
    }

		# For each user we select the point with highest bald score

		pointsToInclude <- tapply(al, userIdsPool, which.max) +
			c(0, cumsum(tapply(al, userIdsPool, length))[ 1 : (length(unique(userIdsPool)) - 1) ])

		# We include these items in the train dataset and eliminate them from the pool set

		train <- rbind(train, pool[ pointsToInclude, ])
		pool <- pool[ -pointsToInclude, ]
	}
}
