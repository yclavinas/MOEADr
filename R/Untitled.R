#' Restricted Neighborhood Replacement Update for MOEA/D
#'
#' Population update using the restricted neighborhood replacement
#' method for the MOEADr package.
#'
#' The restricted neighborhood replacement method behaves like the "standard"
#' replacement method, except that each individual can only be selected up to
#' `nr` times. After this limit has been reached, the next best individual in
#' the same neighborhood is selected.
#'
#' This update routine is intended to be used internally by the main [moead()]
#' function, and should not be called directly by the user.
#'
#' @param update List containing the population update parameters. See
#' Section `Update Strategies` of the [moead()] documentation for
#' details. `update` must contain a field `update$nr`, a positive integer that
#' determines the maximum number of copies of each candidate solution.
#' @param X Matrix of candidate solutions
#' @param Xt Matrix of incumbent solutions
#' @param Y Matrix of objective function values of `X`
#' @param Yt Matrix of objective function values of `Xt`
#' @param B Neighborhood matrix, generated by [define_neighborhood()].
#' @param V List object containing information about the constraint violations
#' of the candidate solutions, generated by [evaluate_population()]
#' @param Vt List object containing information about the constraint violations
#' of the incumbent solutions, generated by [evaluate_population()]
#' @param sel.indx matrix of selection indices, generated by
#' [order_neighborhood()]
#' @param ... other parameters (included for compatibility with generic call)
#'
#' @return List object containing the update population matrix (`X`),
#' and its corresponding matrix of objective function values (`Y`) and
#' constraint value list (`V`).
#'
#' @export
#'
#' @section References:
#' F. Campelo, L.S. Batista, C. Aranha (2020): The {MOEADr} Package: A
#' Component-Based Framework for Multiobjective Evolutionary Algorithms Based on
#' Decomposition. Journal of Statistical Software \doi{10.18637/jss.v092.i06}\cr
#'

updt_restricted <- function(update, X, Xt, Y, Yt, V, Vt, sel.indx, B, ...){
  
  # ========== Error catching and default value definitions
  assertthat::assert_that(
    assertthat::has_name(update,"nr"),
    assertthat::is.count(update$nr))
  
  nr            <- update$nr
  rest.sel.indx <- sel.indx
  
  # Function for returning the selected solution (variable or objectives space)
  # for a subproblem:
  # - i: subproblem index
  # - sel.indx: matrix of selection indices
  # - XY: matrix of candidate solutions (in variable or objective space)
  # - XYt: matrix of incumbent solutions (in variable or objective space)
  # - B: matrix of neighborhoods
  do.update <- function(i, sel.indx, XY, XYt, B) {
    for (j in sel.indx[i,]) { #each element in b_i, in fitness order
      if (j > ncol(B))
        return(XYt[i, , drop = FALSE])     # last row = incumbent solution
      else if (used[B[i, j]] < nr) { # tests if the current element is still available
        used[B[i, j]] <<- used[B[i, j]] + 1 # modifies count matrix in parent env
        return(XY[B[i, j], , drop = FALSE])
      }
    }
  }
  
  # Vector of indices (random permutation), and deshuffling vector
  I  <- sample.int(nrow(X))
  I2 <- order(I)
  # Counter of how many time each solution has been used
  
  ranks <- eaf::pareto_rank(Y)
  temp <- rest.sel.indx
  new.Xnext <- data.frame()
  new.Ynext <- data.frame()
  idxs <- max(ranks):1
  
  used.Y <- rep(0, nrow(X))
  used.X <- rep(0, nrow(X))
  
  for (r in idxs){
    idx <- which(ranks == r)
    rest.sel.indx <- matrix(rest.sel.indx[idx,],  nrow = length(idx))
    I  <- sample.int(length(idx))
    
    used <- used.X
    Xnext <- t(vapply(X         = I,
                      FUN       = do.update,
                      FUN.VALUE = numeric(ncol(X)),
                      sel.indx  = rest.sel.indx,
                      XY        = X,
                      XYt       = Xt,
                      B         = B,
                      USE.NAMES = FALSE))
    new.Xnext <- rbind(new.Xnext, Xnext)
    used.X <- used.X + used
    
    used <- used.Y
    Ynext <- t(vapply(X         = I,
                      FUN       = do.update,
                      FUN.VALUE = numeric(ncol(Y)),
                      sel.indx  = rest.sel.indx,
                      XY        = Y,
                      XYt       = Yt,
                      B         = B,
                      USE.NAMES = FALSE))
    new.Ynext <- rbind(new.Ynext, Ynext)
    used.Y <- used.Y + used
    rest.sel.indx <- temp
  }
  Xnext <- as.matrix(new.Xnext)
  Xnext <- Xnext[I2, ]
  Ynext <- as.matrix(new.Ynext)
  Ynext <- Ynext[I2, ]
  
  
  if(is.null(V)){
    Vnext <- NULL
  } else{
    Vnext <- list(Cmatrix = NULL, Vmatrix = NULL, v = NULL)

    Cmatrix <- V$Cmatrix
    Vmatrix <- V$Vmatrix
    v <- rep(0, length(V$v))
    used.1 <- rep(0, nrow(Y))
    used.2 <- rep(0, nrow(Y))
    for (r in idxs){
      idx <- which(ranks == r)
      rest.sel.indx <- matrix(rest.sel.indx[idx,],  nrow = length(idx))
      # I <- I[idx]
      I  <- sample.int(length(idx))
            
            ## 1: Cmatrix
            used <- used.1
            Vnext$Cmatrix <- t(vapply(X         = I,
                                      FUN       = do.update,
                                      FUN.VALUE = numeric(ncol(V$Cmatrix)),
                                      sel.indx  = rest.sel.indx,
                                      XY        = V$Cmatrix,
                                      XYt       = Vt$Cmatrix,
                                      B         = B,
                                      USE.NAMES = FALSE))
            Cmatrix[idx,] <- Vnext$Cmatrix
            used.1 <- used + used.1
            ## 2: Vmatrix
            used <- used.2
            Vnext$Vmatrix <- t(vapply(X         = I,
                                      FUN       = do.update,
                                      FUN.VALUE = numeric(ncol(V$Vmatrix)),
                                      sel.indx  = rest.sel.indx,
                                      XY        = V$Vmatrix,
                                      XYt       = Vt$Vmatrix,
                                      B         = B,
                                      USE.NAMES = FALSE))
            Vmatrix[idx,] <- Vnext$Vmatrix
            used.2 <- used + used.2
            ## 3: v
            Vnext$v <- rowSums(Vnext$Vmatrix)
            v[idx] <- Vnext$v
            rest.sel.indx <- temp
            used.2 <- used
    }
    # print(v)
    # exit()
    # print("terminou")
    Vnext$v <- v
    Vnext$Vmatrix <- Vmatrix
    Vnext$Cmatrix <- Cmatrix
    
    
  }
  # Output
  return(list(X = Xnext,
              Y = Ynext,
              V = Vnext))
}
