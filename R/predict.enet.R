#' iterxM-fold cross-validation function for Elastic net
#'
#' Estimate test error of elastic net panel
#' @param object an object of a specific machine learning method (ML) class
#' @param validation a string type of validation; either "Mfold" or "loo"
#' @param M an integer specifying the number of folds
#' @param iter an integer specifying the number of times to repeat the cross-validation
#' @param ncores an object of a specific machine learning method (ML) class
#' @param progressBar a boolean that specifies whether a progress bar should be displayed or not
#' @export
predict.enet = function (object, M = 5, iter = 10,
  ncores = 4, progressBar = TRUE) {
  assertthat::assert_that(M > 1, msg = "At least 2 folds (M) are required to run cross-validation!")
  X <- object$xtrain
  y <- object$ytrain
  n <- nrow(X)
  alpha <- object$alpha
  family <- object$family
  lambda <- object$lambda
  lambda_nfolds <- object$lambda_nfolds
  filter <- object$filter
  topranked  <- object$topranked
  keepVar  <- object$keepVar
  weights <- object$weights
  if (M > 1) {
    folds <- lapply(1:iter, function(i) caret::createFolds(y, k = M))
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", ncores))
    parallel::clusterCall(cl, function() library("ssenet"))
    parallel::clusterExport(cl, varlist = c("X", "y", "alpha", "lambda", 'lambda_nfolds', "folds", "progressBar",
      "family", "filter", "topranked", "keepVar", "weights"), envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi,
      X, y, alpha, lambda, lambda_nfolds, progressBar, family, filter,
      topranked, keepVar, weights) {
      ssenet::enetCV(X = X, y = y, alpha = alpha, lambda = lambda,
        lambda_nfolds = lambda_nfolds,
        folds = foldsi, progressBar = progressBar,
        family = family, filter = filter, topranked = topranked,
        keepVar=keepVar, weights = weights)
    }, X, y, alpha, lambda, lambda_nfolds, progressBar, family, filter,
      topranked, keepVar, weights) %>% ssenet::zip_nPure()
    parallel::stopCluster(cl)
    perf <- do.call(rbind, cv$perf) %>% as.data.frame %>%
      tidyr::gather(ErrName, Err) %>% dplyr::group_by(ErrName) %>%
      dplyr::summarise(Mean = mean(Err), SD = sd(Err))
  } else {
    assertthat::assert_that((min(table(y))/M) > 3, msg = "M value too large!")
  }
  result = list()
  result$folds = folds
  result$probs = cv$probs
  result$trueLabels = cv$trueLabels
  result$panels = cv$enet.panel
  result$perf = perf
  return(invisible(result))
}

#' M-fold cross-validation function for Elastic net
#'
#' Estimate test error of elastic net panel
#' @param X nxp matrix - training dataset
#' @param y factor - response variable
#' @param alpha = 1 (lasso), alpha = 0 (ridge), 0 < alpha < 1 (elastic net penalty)
#' @param lambda = strength of elastic net penalty
#' @param folds caret list specifying the indicies in the different folds
#' @param family "binomial" or "multinomial"
#' @param progressBar a boolean that specifies whether a progress bar should be displayed or not
#' @param filter = "none" or "p.value"
#' @param topranked = 50 (top number of features to select and build a classifier)
#' @param keepVar - names of specific variable to keep in model
#' @param weights - observational weights; default to 1
#' @export
enetCV = function (X, y, alpha, lambda, lambda_nfolds, folds, progressBar, family,
  filter, topranked, keepVar, weights) {
  M <- length(folds)
  probs <- predictResponseList <- enet.panel <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    xtrain = X[-omit, , drop = FALSE]
    ytrain = y[-omit]
    xtest <- X[omit, , drop = FALSE]
    ytest = y[omit]

    fit <- ssenet::enet(xtrain=xtrain, ytrain=ytrain, alpha = alpha, lambda = lambda,
      lambda_nfolds = lambda_nfolds,
      family = family, xtest = xtest, ytest = ytest, filter = filter,
      topranked = topranked, keepVar = keepVar, weights = weights[-omit])
    probs[[i]] <- fit$probs
    predictResponseList[[i]] <- fit$predictResponse
    enet.panel[[i]] <- fit$enet.panel
  }
  predictResponse <- unlist(predictResponseList)
  if (family == "binomial") {
    probs <- unlist(probs)
    trueLabels = y[unlist(folds)]
    perf <- ssenet::tperformance(weights = probs, trueLabels = trueLabels)
  }else {
    trueLabels = y[unlist(folds)]
    mat <- table(factor(predictResponse,
      levels(y)), factor(trueLabels, levels(y)))
    mat2 <- mat
    diag(mat2) <- 0
    classError <- colSums(mat2)/colSums(mat)
    er <- sum(mat2)/sum(mat)
    ber <- mean(classError)
    perf <- c(classError, er, ber)
    names(perf) <- c(names(classError), "ER", "BER")
  }
  return(list(probs = probs, trueLabels = trueLabels, perf = perf,
    enet.panel = enet.panel, predictResponse = predictResponse))
}
