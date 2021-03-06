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
#' @importFrom limSolve lsei
#' @import assertthat
#' @import Matrix
#' @import slam
#' @import magrittr
#' @examples
#' data <- iris
#' summary(data)
#' margins <- list(list(var_name = 'Species',value=c(setosa=0.315),
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
adjustment <- function(data, weight = NULL, margins = list(), weight_min = 0.1, weight_max = 10,
  ...) {


  ## sanity check
  assert_that(is.data.frame(data), msg = "data must be a data.frame")

  if (is.null(weight))
    weight <- rep(1, nrow(data))

  assert_that(is.numeric(weight), msg = "weight must be a numeric")
  assert_that(sum(weight) > 0, msg = "sum(weight) must be positive")
  assert_that(nrow(data) == length(weight), msg = "length(weight) must be equal to nrow(data)")

  assert_that(is.numeric(weight_min), msg = "weight_min must be numeric")
  assert_that(length(weight_min) == 1, msg = "length of weight_min must be equal to 1")
  assert_that(weight_min >= 0, msg = "weight_min must be nonnegative")

  assert_that(is.numeric(weight_max), msg = "weight_max must be numeric")
  assert_that(length(weight_max) == 1, msg = "length of weight_max must be equal to 1")
  assert_that(weight_max >= weight_min, msg = "weight_max must be greater or equal to weight_min")

  assert_that(is.list(margins), msg = "margins must be a list")

  for (k in seq_along(margins)) {
    assert_that(is.list(margins[[k]]), msg = "margins must be a list of list (or an empty list)")

    assert_that(has_name(margins[[k]], "var_name"), msg = paste0("margins[[", k, "]] must have a var_name entry"))
    assert_that(has_name(margins[[k]], "value") || has_name(margins[[k]], "min") || has_name(margins[[k]],
      "max"), msg = paste0("margins[[", k, "]] must have a value entry or a min entry or a max entry"))

    assert_that(has_name(data, margins[[k]]$var_name), msg = paste0("data must have a column of name margins[[",
      k, "]]$var_name"))

    if (has_name(margins[[k]], "value")) {
      assert_that(has_attr(margins[[k]]$value, "names"), msg = paste0("margins[[", k,
        "]]$value must have a name attribute"))
      assert_that(all(names(margins[[k]]$value) %in% data[, margins[[k]]$var_name]),
        msg = paste0("each name of margins[[", k, "]]$value  must appears in data[,margins[[",
          k, "]]$var_name]"))
      assert_that(all(is.numeric(margins[[k]]$value)), msg = paste0("margins[[", k, "]]$value must be numeric"))
      assert_that(all(margins[[k]]$value >= 0), msg = paste0("each element of margins[[",
        k, "]]$value must be nonnegative"))
      assert_that(sum(margins[[k]]$value) <= 1, msg = paste0("each element of margins[[",
        k, "]]$value must lower or equal to 1"))
    }

    if (has_name(margins[[k]], "min")) {
      assert_that(has_attr(margins[[k]]$min, "names"), msg = paste0("margins[[", k, "]]$min must have a name attribute"))
      assert_that(all(names(margins[[k]]$min) %in% data[, margins[[k]]$var_name]), msg = paste0("each name of margins[[",
        k, "]]$min  must appears in data[,margins[[", k, "]]$var_name]"))
      assert_that(all(is.numeric(margins[[k]]$min)), msg = paste0("margins[[", k, "]]$min must be numeric"))
      assert_that(all(margins[[k]]$min >= 0), msg = paste0("each element of margins[[",
        k, "]]$min must be nonnegative"))
      assert_that(sum(margins[[k]]$min) <= 1, msg = paste0("each element of margins[[",
        k, "]]$min must lower or equal to 1"))
    }

    if (has_name(margins[[k]], "max")) {
      assert_that(has_attr(margins[[k]]$max, "names"), msg = paste0("margins[[", k, "]]$max must have a name attribute"))
      assert_that(all(names(margins[[k]]$max) %in% data[, margins[[k]]$var_name]), msg = paste0("each name of margins[[",
        k, "]]$max  must appears in data[,margins[[", k, "]]$var_name]"))
      assert_that(all(is.numeric(margins[[k]]$max)), msg = paste0("margins[[", k, "]]$max must be numeric"))
      assert_that(all(margins[[k]]$max >= 0), msg = paste0("each element of margins[[",
        k, "]]$max must be nonnegative"))
      # assert_that(sum(margins[[k]]$max) <= 1,msg=paste0('each element of margins[[',k,']]$max
      # must lower or equal to 1'))# not relevant if the number of modality is not complete
    }





  }


  ## Building matrix representations
  E <- Matrix(nrow = 0, ncol = nrow(data))
  F <- NULL
  G <- Matrix(nrow = 0, ncol = nrow(data))
  H <- NULL

  for (k in seq_along(margins)) {

    M <- sparse.model.matrix(~x - 1, data.frame(x = as.character(data[, margins[[k]]$var_name])))
    colnames(M) <- gsub("^x", "", colnames(M))
    if (length(margins[[k]]$value) > 0) {
      E <- rBind(E, t(M[, names(margins[[k]]$value), drop = FALSE]))
      F <- c(F, margins[[k]]$value)
    }

    if (length(margins[[k]]$min) > 0) {
      G <- rBind(G, t(M[, names(margins[[k]]$min), drop = FALSE]))
      H <- c(H, margins[[k]]$min)
    }

    if (length(margins[[k]]$max) > 0) {
      G <- rBind(G, -t(M[, names(margins[[k]]$max), drop = FALSE]))
      H <- c(H, -margins[[k]]$max)
    }
  }

  E <- E/sum(weight)
  G <- G/sum(weight)


  G <- rBind(G, Diagonal(x = rep(1, nrow(data))), -Diagonal(x = rep(1, nrow(data))))
  H <- c(H, rep(weight_min, nrow(data)), -rep(weight_max, nrow(data)))

  E <- rBind(E, matrix(1, nrow = 1, ncol = nrow(data)))
  F <- c(F, sum(weight))

  A <- Diagonal(x = rep(1, nrow(data)))
  B <- weight

  # E<-as(E,'dgTMatrix') E<-simple_triplet_matrix(i = E@i + 1, j = E@j + 1, v = E@x,
  # nrow=E@Dim[1], ncol=E@Dim[2]) G<-as(G,'dgTMatrix') G<-simple_triplet_matrix(i = G@i + 1,
  # j = G@j + 1, v = G@x, nrow=G@Dim[1], ncol=G@Dim[2])

  model <- lsei(A = A %>% as.matrix, B = B %>% as.matrix, E = E %>% as.matrix, F = F %>%
    as.matrix, G = G %>% as.matrix, H = H %>% as.matrix, ...)

  list(w = model$X, IsError = model$IsError)
}
