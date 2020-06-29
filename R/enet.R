#' Elastic net classification panel
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
enet = function (xtrain, ytrain, alpha, lambda = NULL, lambda_nfolds=3, family, xtest = NULL,
  ytest = NULL, filter = "none", topranked = 50, keepVar = NULL, weights = NULL) {

  # Check format of data input
  assertthat::assert_that(class(xtrain)[1] == "matrix", msg = "xtrain must be a matrix!")
  assertthat::assert_that(ifelse(is.null(xtest), TRUE, class(xtest) == "matrix"), msg = "xtest must be a matrix!")
  assertthat::assert_that(class(ytrain)[1] == "factor", msg = "ytrain must be a factor!")
  assertthat::assert_that(ifelse(is.null(ytest), TRUE, class(ytest) == "factor"), msg = "ytest must be a factor!")

  # Check that the number of samples in the minority classes has more than lambda_nfolds
  assertthat::assert_that(min(table(ytrain)) > lambda_nfolds, msg = "Minority class has less samples than lambda_nfolds!")

  # if observations weights are zero to NULL, default to 1 for each obseration
  if(is.null(weights)){
    weights <- rep(1, length(ytrain))
  }

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

  # set model parameters
  glmnet_default_params <- list(x = x2train, y = ytrain, alpha = alpha,
    penalty.factor = penalty.factor, weights = weights)
  cvglmnet_default_params <- list(x = x2train, y = ytrain, weights = weights, nfolds = lambda_nfolds)
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
  cv.fit <- do.call(glmnet::cv.glmnet, c(cvglmnet_default_params, args_cvglmnet))
  lambda <- ifelse(is.null(lambda), cv.fit$lambda.1se, lambda)
  features <- ssenet::extractFeatures(fit, lambda, family)


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
    lambda = lambda, lambda_nfolds=lambda_nfolds, alpha = alpha, family = family, probs = probs,
    Active.Coefficients = features$Active.Coefficients, perfTest = perfTest,
    predictResponse = predictResponse, filter = filter, topranked = topranked, keepVar=keepVar, weights = weights)
  class(result) <- "enet"
  return(result)
}
