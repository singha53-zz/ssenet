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
enet = function (xtrain, ytrain, alpha, lambda = NULL, family, xtest = NULL,
  ytest = NULL, filter = "p.value", topranked = 50, keepVar = NULL, weights = NULL) {

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

  # Fit glmnet model for binary response
  if (family == "binomial") {
    fit <- glmnet::glmnet(x2train, ytrain, family = "binomial", alpha = alpha,
      penalty.factor = penalty.factor, weights = weights)
    if (is.null(lambda)) {
      cv.fit <- glmnet::cv.glmnet(x2train, ytrain, family = "binomial", weights = weights)
      lambda = cv.fit$lambda.min
    } else {
      lambda = lambda
    }
    Coefficients <- coef(fit, s = lambda)
    Active.Index <- which(Coefficients[, 1] != 0)
    Active.Coefficients <- Coefficients[Active.Index, ]
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
  }

  # Fit glmnet model for multi-class response
  if (family == "multinomial") {
    fit <- glmnet::glmnet(x2train, ytrain, family = "multinomial", alpha = alpha,
      type.multinomial = "grouped", penalty.factor = penalty.factor, weights = weights)
    if (is.null(lambda)) {
      cv.fit <- glmnet::cv.glmnet(x2train, ytrain, family = "multinomial", weights = weights)
      lambda = cv.fit$lambda.min
    } else {
      lambda = lambda
    }
    Coefficients <- coef(fit, s = lambda)
    Active.Index <- which(Coefficients[[1]][, 1] != 0)
    Active.Coefficients <- Coefficients[[1]][Active.Index,]
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
  }

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
  result <- list(xtrain = xtrain, ytrain = ytrain, fit = fit, enet.panel = enet.panel,
    lambda = lambda, alpha = alpha, family = family, probs = probs,
    Active.Coefficients = Active.Coefficients, perfTest = perfTest,
    predictResponse = predictResponse, filter = filter, topranked = topranked, keepVar=keepVar, weights = weights)
  class(result) <- "enet"
  return(result)
}
