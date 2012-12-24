#
# R script that evaluates the performance of the multitask preference learning model on synthetic data.
#
# Author: Jose Miguel Hernandez Lobato
#
# Date: 24 Dec 2011
#

set.seed(1)

source("epBonilla.R")

# We eliminate the previous results

system("rm -f results/*.txt")

# We load the item descriptions

itemFeatures <- as.matrix(read.table("../trainTestPartitions/itemFeatures.txt"))

# We load the user descriptions

userFeatures <- as.matrix(read.table("../trainTestPartitions/userFeatures.txt"))

# We standardize the user features so that they have zero mean and unit standard deviation

meanUserFeatures <- apply(userFeatures, 2, mean)
sdUserFeatures <- apply(userFeatures, 2, sd)
sdUserFeatures[ sdUserFeatures == 0 ] <- 1

userFeatures <- (userFeatures - matrix(meanUserFeatures, nrow(userFeatures), ncol(userFeatures), byrow = T)) /
	matrix(sdUserFeatures, nrow(userFeatures), ncol(userFeatures), byrow = T)

# We load the number of different items

nItems <- read.table("../trainTestPartitions/nItems.txt")$V1

# We load the simulation number

simNumber <- 1

# We repeat for the current train/test partition

for (i in simNumber) {

	# We load the training, pool and test sets

	train <- as.matrix(read.table(paste("../trainTestPartitions/train", i, ".txt", sep = "")))
	test <- as.matrix(read.table(paste("../trainTestPartitions/test", i, ".txt", sep = "")))
	pool <- as.matrix(read.table(paste("../trainTestPartitions/pool", i, ".txt", sep = "")))

	train <- train[ which(train[ , 3 ] <= 25), ]

	# We do 10 iterations of active learning

	nQueries <- 0
	for (j in 1 : (nQueries + 1)) {

		# We standardize the item features so that they have zero mean and unit standard deviation in the training set

		meanItemFeatures <- apply(itemFeatures[ unique(c(train[ , 1 ], train[ , 2 ])), ], 2, mean)
		sdItemFeatures <- apply(itemFeatures[ unique(c(train[ , 1 ], train[ , 2 ])), ], 2, sd)
		sdItemFeatures[ sdItemFeatures == 0 ] <- 1

		itemFeaturesStandardized <- (itemFeatures - matrix(meanItemFeatures, nrow(itemFeatures), ncol(itemFeatures), byrow = T)) / 
			matrix(sdItemFeatures, nrow(itemFeatures), ncol(itemFeatures), byrow = T)

		itemFeaturesStandardized <- itemFeatures

		# We create the training data to pass to the EP method

		itemFeaturesAtrain <- itemFeaturesStandardized[ train[ , 1 ], ]
		itemFeaturesBtrain <-itemFeaturesStandardized[ train[ , 2 ], ]
		userIdsTrain <- train[ , 3 ]
		ratingsTrain <- train[ , 4 ]

		# We construct the problemInfo list

		problemInfo <- list(itemFeaturesA = itemFeaturesAtrain, itemFeaturesB = itemFeaturesBtrain,
			ratings = ratingsTrain, userFeatures = userFeatures[ userIdsTrain, ])

		# We call the EP method

		ret <- epBonilla(problemInfo)

		# We store the kernel hyper-parameters

		write.table(ret$problemInfo$lengthScaleUsers, "lengthScaleUsers.txt", col.names = F, row.names = F)
		write.table(ret$problemInfo$lengthScaleItems, "lengthScaleItems.txt", col.names = F, row.names = F)
		write.table(ret$problemInfo$noiseUsers, "noiseUsers.txt", col.names = F, row.names = F)
		write.table(ret$problemInfo$noiseItems, "noiseItems.txt", col.names = F, row.names = F)

		# We eavaluate the performance of the method on the test set

		userIdsTest <- test[ , 3 ]
		ratingsTest <- test[ , 4 ]
		itemFeaturesAtest <- itemFeaturesStandardized[ test[ , 1 ], ]
		itemFeaturesBtest <- itemFeaturesStandardized[ test[ , 2 ], ]

		pred <- sign(predictBonilla(ret, userFeatures[ userIdsTest, ], itemFeaturesAtest, itemFeaturesBtest)$m)
		error <- mean(pred != ratingsTest)
		write.table(error, paste("results/error", j - 1, ".txt", sep = ""), col.names = F, row.names = F, append = T)

		# We make predictions in the pool set

		userIdsPool <- pool[ , 3 ]
		itemFeaturesApool <- itemFeaturesStandardized[ pool[ , 1 ], ]
		itemFeaturesBpool <- itemFeaturesStandardized[ pool[ , 2 ], ]

		pred <- predictBonilla(ret, userFeatures[ userIdsPool, ], itemFeaturesApool, itemFeaturesBpool)

		# We compute the bald score for each point in the pool set

		prob <- pnorm(pred$m / sqrt(pred$v + 1))
		firstTerm <- -(1 - prob) * log(1 - prob) - prob * log(prob)
		C <- sqrt(pi * log(2)  / 2)
		secondTerm <- C / (sqrt(pred$v) + C^2) * exp(-pred$m^2 / (2 * (pred$v + C^2)))

		bald <- firstTerm - secondTerm

		# For each user we select the point with highest bald score

		pointsToInclude <- tapply(bald, userIdsPool, which.max) +
			c(0, cumsum(tapply(bald, userIdsPool, length))[ 1 : (length(unique(userIdsPool)) - 1) ])

		# We include these items in the train dataset and eliminate them from the pool set

		train <- rbind(train, pool[ pointsToInclude, ])
		pool <- pool[ -pointsToInclude, ]
	}
}
