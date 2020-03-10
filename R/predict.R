#' cross-validation function for all methods
#'
#' Estimate test error of elastic net panel
#' @param object an object of a specific machine learning method (ML) class
#' @param validation a string type of validation; either "Mfold" or "loo"
#' @param M an integer specifying the number of folds
#' @param iter an integer specifying the number of times to repeat the cross-validation
#' @param threads an object of a specific machine learning method (ML) class
#' @param progressBar a boolean that specifies whether a progress bar should be displayed or not
#' @export
predict <- function(object, ...){
  UseMethod("predict")
}
