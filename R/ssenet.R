#' Semi-supervised Elastic net classification panel
#'
#' build an elastic net classification panel
#' @param xtrain nxp matrix - training dataset
#' @param ytrain factor - response variable
#' @param alpha = 1 (lasso), alpha = 0 (ridge), 0 < alpha < 1 (elastic net penalty)
#' @param lambda = strength of elastic net penalty
#' @param family "binomial" or "multinomial"
#' @param xtest nxp matrx  - test dataset
#' @param ytest factor - response variable
#' @param filter = "none" or "p.value"
#' @param topranked = 50 (top number of features to select and build a classifier)
#' @param keepVar - names of specific variable to keep in model
#' @param weights - observational weights; default to 1
#' @export
ssenet = function (xtrain, ytrain, alpha, lambda = NULL, lambda_nfolds=3, family, xtest = NULL,
  ytest = NULL, filter = "none", topranked = 50, keepVar = NULL, useObsWeights = FALSE, max.iter = 100, perc.full = 1, thr.conf = 0.5) {

  # Check that ytrain has missing data
  assertthat::assert_that(sum(is.na(ytrain)) > 0, msg = "ytrain must have some missing labels")

  # Check format of data input
  assertthat::assert_that(class(xtrain) == "matrix", msg = "xtrain must be a matrix!")
  assertthat::assert_that(ifelse(is.null(xtest), TRUE, class(xtest) == "matrix"), msg = "xtest must be a matrix!")
  assertthat::assert_that(class(ytrain) == "factor", msg = "ytrain must be a factor!")
  assertthat::assert_that(ifelse(is.null(ytest), TRUE, class(ytest) == "factor"), msg = "ytest must be a factor!")

  # Check that the number of samples in the minority classes has more than lambda_nfolds
  assertthat::assert_that(min(table(ytrain)) > lambda_nfolds, msg = "Minority class has less samples than lambda_nfolds!")

  # Check that all classes have at least 3 observation
  assertthat::assert_that(any(table(ytrain) > 3), msg = "3 obs must be present for a given cell-type")

  # if observations weights are zero to NULL, default to 1 for each obseration
  if(!useObsWeights){
    weights <- rep(1, length(ytrain))
  }

  # # lambda is preset, since cv.glmnet doesnt work if sample size in the individual classes is too low.
  # if (is.null(lambda)) {
  #   lambda = 0.01
  #   print("No lambda value provided; using a lambda value of 0.01")
  # } else {
  #   lambda = lambda
  # }

  # Pre-filtering xtrain (none, pvalue, or top features)
  if (filter == "none") {
    x1train <- xtrain
  }
  if (filter == "p.value") {
    design <- model.matrix(~ytrain)
    fit <- limma::eBayes(limma::lmFit(t(xtrain), design))
    top <- limma::topTable(fit, coef = 2, adjust.method = "BH", n = nrow(fit))
    x1train <- xtrain[, rownames(top)[1:topranked]]
  }
  if (is.null(keepVar)) {
    penalty.factor <- rep(1, ncol(x1train))
    x2train <- x1train
  } else {
    x1train <- x1train[, setdiff(colnames(x1train), keepVar)]
    x2train <- as.matrix(cbind(x1train, xtrain[, keepVar]))
    colnames(x2train) <- c(colnames(x1train), keepVar)
    penalty.factor <- c(rep(1, ncol(x1train)), rep(0, length(keepVar)))
  }

  if (is.null(lambda)) {
    cv.fit <- glmnet::cv.glmnet(x2train[!is.na(ytrain), ],
                                ytrain[!is.na(ytrain)],
                                family="multinomial",
                                type.multinomial = "grouped")
    lambda = cv.fit$lambda.min
  } else {
    lambda = lambda
  }

  ytrainImputed <- ssenet::imputeLabels(x2train, ytrain, alpha, lambda, lambda_nfolds, family, max.iter, perc.full, thr.conf, useObsWeights)
  assertthat::assert_that(any(table(ytrainImputed$pred) > 3), msg = "3 obs must be present for a given cell-type")

  # remove sample with 1 or 0 observations
  keepObs <- as.character(ytrainImputed$pred) %in% names(table(ytrainImputed$pred))[table(ytrainImputed$pred) > 1]

  # set model parameters
  glmnet_default_params <- list(x = xtrain[keepObs, ],
                                y = factor(ytrainImputed$pred[keepObs]),
                                alpha = alpha,
                                penalty.factor = penalty.factor,
                                weights = ytrainImputed$weights[keepObs])
  cvglmnet_default_params <- list(x = xtrain[keepObs, ],
                                  y = factor(ytrainImputed$pred[keepObs]),
                                  weights = ytrainImputed$weights[keepObs],
                                  nfolds = lambda_nfolds)
  if(family == "binomial"){
    args_glmnet <- list(family = "binomial")
    args_cvglmnet <- list(family = "binomial")
  }
  if(family == "multinomial"){
    args_glmnet <- list(family = "multinomial", type.multinomial = "grouped")
    args_cvglmnet <- list(family = "multinomial")
  }

  # Fit model
  fit <- do.call(glmnet::glmnet, c(glmnet_default_params, args_glmnet))
  # cv.fit <- do.call(glmnet::cv.glmnet, c(cvglmnet_default_params, args_cvglmnet))
  lambda <- ifelse(is.null(lambda), do.call(glmnet::cv.glmnet, c(cvglmnet_default_params, args_cvglmnet))$lambda.min, lambda)
  features <- extractFeatures(fit, lambda, family)

  # Apply model to test data if provided
  if (!is.null(xtest)) {
    probs <- glmnet::predict.multnet(fit, newx = xtest[, colnames(x2train)], s = lambda, type = "response")
    predictResponse <- unlist(glmnet::predict.multnet(fit, newx = xtest, s = lambda, type = "class"))
    if (family == "binomial") {
      perfTest <- ssenet::tperformance(weights = as.numeric(as.matrix(probs)), trueLabels = ytest)
    } else {
      mat <- table(factor(as.character(predictResponse),
        levels = levels(ytest)), ytest)
      mat2 <- mat
      diag(mat2) <- 0
      classError <- colSums(mat2)/colSums(mat)
      er <- sum(mat2)/sum(mat)
      ber <- mean(classError)
      perfTest <- c(classError, er, ber)
      names(perfTest) <- c(names(classError), "ER", "BER")
    }
  } else {
    perfTest <- predictResponse <- probs <- NA
  }
  result <- list(xtrain = xtrain, ytrain = ytrain, fit = fit, enet.panel = features$enet.panel,
    lambda = lambda, lambda_nfolds = lambda_nfolds, alpha = alpha, family = family, probs = probs,
    Active.Coefficients = features$Active.Coefficients, perfTest = perfTest, ytrainImputed=ytrainImputed,
    predictResponse = predictResponse, filter = filter, topranked = topranked, keepVar=keepVar,
    useObsWeights = useObsWeights, max.iter = max.iter, perc.full = perc.full, thr.conf = thr.conf)
  class(result) <- "ssenet"
  return(result)
}


# implemented using https://github.com/mabelc/SSC/blob/master/R/SelfTraining.R
#' @export
imputeLabels = function(x, y, alpha, lambda, lambda_nfolds, family, max.iter, perc.full, thr.conf, useObsWeights){
  ### Init variables ###

  # Identify the classes
  classes <- levels(y)
  nclasses <- length(classes)

  # Init variable to store the labels
  ynew <- y

  # Obtain the indexes of labeled and unlabeled instances
  labeled <- which(!is.na(y))
  unlabeled <- which(is.na(y))

  ## Check the labeled and unlabeled sets
  if(length(labeled) == 0){   # labeled is empty
    stop("The labeled set is empty. All the values in y parameter are NA.")
  }
  if(length(unlabeled) == 0){ # unlabeled is empty
    stop("The unlabeled set is empty. None value in y parameter is NA.")
  }

  ### Self Training algorithm ###

  # Count the examples per class
  cls.summary <- table(y[labeled])
  # Determine the total of instances to include per iteration
  cantClass <- round(cls.summary / min(cls.summary))
  totalPerIter <- sum(cantClass)
  # Compute count full
  count.full <- length(labeled) + round(length(unlabeled) * perc.full)

  iter <- 1
  while ((length(labeled) < count.full) && (length(unlabeled) >= totalPerIter) && (iter <= max.iter)) {

    if(!useObsWeights){
      weights <- rep(1, length(ynew[labeled]))
    } else {
      weights <- as.numeric((1/table(ynew[labeled]))[as.character(ynew[labeled])])
    }

    assertthat::assert_that(any(table(ynew[labeled]) > 3), msg = "3 obs must be present for a given class")

    # Train classifier
    if(family == "binomial"){
      model <- glmnet::glmnet(x[labeled, ], ynew[labeled], family = family, alpha = alpha, weights = weights[labeled])
      if (is.null(lambda)) {
        cvob1 = glmnet::cv.glmnet(x[labeled, ], ynew[labeled], family = family, weights = weights[labeled])
        lambda = cvob1$lambda.min
      } else {
        lambda = lambda
      }
      # Predict probabilities per classes of unlabeled examples
      prob <- checkProb(prob = glmnet::predict.lognet(model, x[unlabeled, ], s = lambda, type = "response")[,,1], ninstances = length(unlabeled), classes)
    }
    if(family == "multinomial"){
      model <- glmnet::glmnet(x[labeled, ],
                              ynew[labeled],
                              family=family,
                              alpha = alpha,
                              type.multinomial = "grouped",
                              weights = weights)
      if (is.null(lambda)) {
        cvob1 = glmnet::cv.glmnet(x[labeled, ],
                                  ynew[labeled],
                                  family=family,
                                  type.multinomial = "grouped",
                                  weights = weights)
        lambda = cvob1$lambda.min
      } else {
        lambda = lambda
      }
      # Predict probabilities per classes of unlabeled examples
      prob <- checkProb(prob = glmnet::predict.multnet(model, x[unlabeled, ],
                                                       s = lambda,
                                                       type = "response",
                                                       type.multinomial = "grouped")[,,1],
                        ninstances = length(unlabeled),
                        classes)
    }

    # Select the instances with better class probability
    pre.selection <- selectInstances(cantClass, prob)

    # Select the instances with probability grather than the theshold confidence
    indexes <- which(pre.selection$prob.cls > thr.conf)
    if(length(indexes) == 0){
      iter <- iter + 1
      next
    }
    selection <- pre.selection[indexes,]

    # Add selected instances to L
    labeled.prime <- unlabeled[selection$unlabeled.idx]
    sel.classes <- classes[selection$class.idx]
    ynew[labeled.prime] <- sel.classes
    labeled <- c(labeled, labeled.prime)

    # Delete selected instances from U
    unlabeled <- unlabeled[-selection$unlabeled.idx]

    iter <- iter + 1
  }

  if(!useObsWeights){
    weights <- rep(1, length(ynew[labeled]))
  } else {
    weights <- as.numeric((1/table(ynew[labeled]))[as.character(ynew[labeled])])
  }

  if(family == "binomial"){
    model <- glmnet::glmnet(x[labeled, ], ynew[labeled], family=family, alpha = alpha, weights = weights)
    if (is.null(lambda)) {
      cvob1 = glmnet::cv.glmnet(x[labeled, ], ynew[labeled], family = family, weights = weights)
      lambda = cvob1$lambda.min
    } else {
      lambda = lambda
    }
    pred <- predict(model, x, s = lambda, type = "class")
  }
  if(family == "multinomial"){
    model <- glmnet::glmnet(x[labeled, ],
                            ynew[labeled],
                            family=family,
                            alpha = alpha,
                            type.multinomial = "grouped",
                            weights = weights)
    if (is.null(lambda)) {
      cvob1 = glmnet::cv.glmnet(x[labeled, ],
                                ynew[labeled],
                                family=family,
                                type.multinomial = "grouped",
                                weights = weights)
      lambda = cvob1$lambda.min
    } else {
      lambda = lambda
    }
    pred <- glmnet::predict.multnet(model, x, s = lambda, type = "class", type.multinomial = "grouped")
  }
  #final weights
  weights <- as.numeric((1/table(pred))[pred])
  return(list(pred=pred, weights=weights))
}
