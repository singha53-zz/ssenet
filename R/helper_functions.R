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
