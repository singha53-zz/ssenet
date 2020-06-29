checkProb <- function(prob, ninstances, classes){
  # Check probabilities matrix
  if(!is.matrix(prob)){
    stop(
      sprintf(
        paste0(
          "Predict function incorrect output.\n",
          "'prob' is an object of class %s.\n",
          "Expected an object of class matrix."
        ),
        class(prob)
      )
    )
  }
  if(ninstances != nrow(prob)){
    stop(
      sprintf(
        paste0(
          "Predict function incorrect output.\n",
          "The row number of 'prob' is %s.\n",
          "Expected a number equal to %i (value of 'ninstances')."
        ),
        nrow(prob),
        ninstances)
    )
  }
  if(length(classes) != ncol(prob)){
    stop(
      sprintf(
        paste0(
          "Predict function incorrect output.\n",
          "The column number of 'prob' is %s.\n",
          "Expected a number equal to %i (length of 'classes')."
        ),
        ncol(prob),
        length(classes))
    )
  }
  if(length(classes) != length(intersect(classes, colnames(prob)))){
    stop(
      paste0(
        "Predict function incorrect output.\n",
        "The columns names of 'prob' is a set not equal to 'classes' set."
      )
    )
  } else {
    # order columns by classes
    prob <- prob[, classes]
    if(!is.matrix(prob)){
      # when nrow of prob is 1
      prob <- matrix(prob, nrow = 1)
      colnames(prob) <- classes
    }
  }

  return(prob)
}


selectInstances <- function(cantClass, probabilities){
  len <- 0
  class.idx <- numeric()
  unlabeled.idx <- numeric()
  prob.cls <- numeric()

  for (k in 1:sum(cantClass)) { # buscar el mejor por clase y etiquetarlo
    best <- arrayInd(which.max(probabilities), dim(probabilities))
    i <- best[1] # fila (instancia)
    c <- best[2] # columna (clase)
    if (probabilities[i,c] == -1){
      break;
    }

    if (cantClass[c] > 0) {
      len <- len + 1
      class.idx[len] <- c
      unlabeled.idx[len] <- i
      prob.cls[len] <- probabilities[i, c]

      cantClass[c] <- cantClass[c] - 1
      probabilities[i,] <- -1 # para que no se repita la instancia
      if (cantClass[c] == 0)
        probabilities[,c] <- -1 # para que no se repita la clase
    }

  }

  r <- data.frame(class.idx = class.idx, unlabeled.idx = unlabeled.idx, prob.cls = prob.cls)
  return(r)
}

#' @export
customTheme = function (sizeStripFont, xAngle, hjust, vjust, xSize, ySize,
  xAxisSize, yAxisSize) {
  theme(strip.background = element_rect(colour = "black", fill = "white",
    size = 1), strip.text.x = element_text(size = sizeStripFont),
    strip.text.y = element_text(size = sizeStripFont), axis.text.x = element_text(angle = xAngle,
      hjust = hjust, vjust = vjust, size = xSize, color = "black"),
    axis.text.y = element_text(size = ySize, color = "black"),
    axis.title.x = element_text(size = xAxisSize, color = "black"),
    axis.title.y = element_text(size = yAxisSize, color = "black"),
    panel.background = element_rect(fill = "white", color = "black"))
}

#' @export
zip_nPure = function (.x, .fields = NULL, .simplify = FALSE) {
  if (length(.x) == 0)
    return(list())
  if (is.null(.fields)) {
    if (is.null(names(.x[[1]]))) {
      .fields <- seq_along(.x[[1]])
    }
    else {
      .fields <- stats::setNames(names(.x[[1]]), names(.x[[1]]))
    }
  }
  else {
    if (is.character(.fields) && is.null(names(.fields))) {
      names(.fields) <- .fields
    }
  }
  out <- lapply(.fields, function(i) lapply(.x, .subset2, i))
  if (.simplify)
    out <- lapply(out, simplify_if_possible)
  out
}

#' @export
tperformance = function (weights, trueLabels)
{
  df = data.frame(prob = as.numeric(weights), status = model.matrix(~factor(as.character(trueLabels),
    levels = levels(trueLabels)))[, 2])
  roc.score = roc(response = df$status, predictor = weights,
    plot = FALSE, percent = TRUE, na.rm = TRUE, direction = "<")
  optimal.cutpoint.Youden <- optimal.cutpoints(X = "prob",
    status = "status", tag.healthy = 0, methods = "Youden",
    data = df, control = control.cutpoints(), ci.fit = FALSE,
    conf.level = 0.95, trace = FALSE, pop.prev = 0.5)
  optimalValues <- round(c(summary(optimal.cutpoint.Youden)$p.table$Global$Youden[[1]][1:5,
    ], roc.score$auc/100), 3)
  names(optimalValues) <- c(names(optimalValues)[-length(names(optimalValues))],
    "AUC")
  optimalValues
}

#' @export
extractFeatures = function(fit, lambda, family){
  Coefficients <- coef(fit, s = lambda)
  if(family == "binomial"){
    Active.Index <- which(Coefficients[, 1] != 0)
    Active.Coefficients <- Coefficients[Active.Index, ]
  }
  if(family == "multinomial"){
    Active.Index <- which(Coefficients[[1]][, 1] != 0)
    Active.Coefficients <- Coefficients[[1]][Active.Index,]
  }
  enet.panel <- names(Active.Coefficients)[-1]
  enet.panel.length <- length(enet.panel)
  return(list(Active.Coefficients=Active.Coefficients, enet.panel=enet.panel))
}

#' @export
estimate_error <- function(truth, pred) {
  mat <- table(factor(as.character(pred), levels = levels(truth)), truth)
  mat2 <- mat
  diag(mat2) <- 0
  classError <- colSums(mat2)/colSums(mat)
  er <- sum(mat2)/sum(mat)
  ber <- mean(classError)
  perfTest <- c(classError, er, ber)
  names(perfTest) <- c(names(classError), "ER", "BER")
  perfTest
}

#' @export
set_na <- function(y, perc_na){
  unlist(lapply(levels(y), function(i){
    na_index <- sample(x = which(y == i), size = ceiling(length(which(y == i)) * perc_na))
  }))
}
