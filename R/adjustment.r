#' A function to calculate satistical adjustment weight
#'
#' This function allows to calculate satistical adjustment weight given some margins (equalities or inequalities).
#' @param data a data frame.
#' @param weight a nonnegative numeric vector (with positive sum) of length \code{nrow(data)}.
#' @param margins a list specifing the margins constrains.
#' @param weight_max a nonnegative numeric constant.
#' @param weight_min a numeric constant greater or equal to 1.
#' @keywords adjustment calibration weighting
#' @export
#' @examples
#' data <- iris
#' summary(data)
#' margins <- list(list(var_name = "Species",value=c(setosa=0.315),
#' min=c(versicolor=0.45),
#' max=c(versicolor=0.01)))

#' w <- adjustment(data = data,margins = margins)$w
#' data$w<-w
#' library(dplyr)
#' data%>%
#'   group_by(Species)%>%
#'   summarise(eff=sum(w))%>%
#'   merge(data%>%summarise(eff0=sum(w)))%>%
#'   mutate(freq = eff/eff0)
adjustment <- function(data, weight = NULL, margins = list(), 
  weight_min = 0.1, weight_max = 10, ...) {
  
  
  ## sanity check
  assert_that(is.data.frame(data))
  
  if (is.null(weight)) 
    weight <- rep(1, nrow(data))
  
  assert_that(is.numeric(weight))
  assert_that(sum(weight) > 0)
  assert_that(nrow(data) == length(weight))
  
  assert_that(is.numeric(weight_min))
  assert_that(length(weight_min) == 1)
  assert_that(weight_min >= 0)
  
  assert_that(is.numeric(weight_max))
  assert_that(length(weight_max) == 1)
  assert_that(weight_max >= 1)
  
  assert_that(is.list(margins))
  
  for (k in seq_along(margins)) {
    assert_that(is.list(margins[[k]]))
    
    assert_that(has_name(margins[[k]], "var_name"))
    assert_that(has_name(margins[[k]], "value") || has_name(margins[[k]], 
      "min") || has_name(margins[[k]], "max"))
    
    assert_that(has_name(data, margins[[k]]$var_name))
    
    if (has_name(margins[[k]], "value")) {
      assert_that(has_attr(margins[[k]]$value, "names"))
      assert_that(all(names(margins[[k]]$value) %in% data[, 
        margins[[k]]$var_name]))
      assert_that(all(is.numeric(margins[[k]]$value)))
      assert_that(all(margins[[k]]$value >= 0))
      assert_that(sum(margins[[k]]$value) <= 1)
    }
    
    if (has_name(margins[[k]], "min")) {
      assert_that(has_attr(margins[[k]]$min, "names"))
      assert_that(all(names(margins[[k]]$min) %in% data[, 
        margins[[k]]$var_name]))
      assert_that(all(is.numeric(margins[[k]]$min)))
      assert_that(all(margins[[k]]$min >= 0))
      assert_that(sum(margins[[k]]$min) <= 1)
    }
    
    
    if (has_name(margins[[k]], "max")) {
      assert_that(has_attr(margins[[k]]$max, "names"))
      assert_that(all(names(margins[[k]]$max) %in% data[, 
        margins[[k]]$var_name]))
      assert_that(all(is.numeric(margins[[k]]$max)))
      assert_that(all(margins[[k]]$max >= 0))
      # assert_that(sum(margins[[k]]$max)<=1)# this condition is
      # not valid if the number of specified margins is lower than
      # the number of modalities in the data
    }
    
  }
  
  
  ## Building matrix representations
  E <- Matrix(nrow = 0, ncol = nrow(data))
  F <- NULL
  G <- Matrix(nrow = 0, ncol = nrow(data))
  H <- NULL
  
  
  for (k in seq_along(margins)) {
    
    M <- sparse.model.matrix(~x - 1, data.frame(x = as.character(data[, 
      margins[[k]]$var_name])))
    colnames(M) <- gsub("^x", "", colnames(M))
    if (length(margins[[k]]$value) > 0) {
      E <- rbind(E, t(M[, names(margins[[k]]$value), drop = FALSE]))
      F <- c(F, margins[[k]]$value)
    }
    
    if (length(margins[[k]]$min) > 0) {
      G <- rbind(G, t(M[, names(margins[[k]]$min), drop = FALSE]))
      H <- c(H, margins[[k]]$min)
    }
    
    if (length(margins[[k]]$max) > 0) {
      G <- rbind(G, -t(M[, names(margins[[k]]$max), drop = FALSE]))
      H <- c(H, -margins[[k]]$max)
    }
  }
  
  E <- E/sum(weight)
  G <- G/sum(weight)
  
  
  G <- rbind(G, Diagonal(x = rep(1, nrow(data))), -Diagonal(x = rep(1, 
    nrow(data))))
  H <- c(H, rep(weight_min, nrow(data)), -rep(weight_max, nrow(data)))
  
  E <- rbind(E, matrix(1, nrow = 1, ncol = nrow(data)))
  F <- c(F, sum(weight))
  
  A <- Diagonal(x = rep(1, nrow(data)))
  B <- weight
  
  # E<-as(E,'dgTMatrix') E<-simple_triplet_matrix(i = E@i + 1,
  # j = E@j + 1, v = E@x, nrow=E@Dim[1], ncol=E@Dim[2])
  # G<-as(G,'dgTMatrix') G<-simple_triplet_matrix(i = G@i + 1,
  # j = G@j + 1, v = G@x, nrow=G@Dim[1], ncol=G@Dim[2])
  
  model <- lsei(A = A %>% as.matrix, B = B %>% as.matrix, E = E %>% 
    as.matrix, F = F %>% as.matrix, G = G %>% as.matrix, 
    H = H %>% as.matrix, ...)
  list(w = model$X, IsError = model$IsError)
}