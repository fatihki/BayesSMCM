#' @title Cumulative hazard and hazard functions under the piecewise exponential baseline hazard model
#'
#' @description
#' Calculates the cumulative hazard and hazard functions defined by the piecewise exponential baseline hazard model.

#' @name cumulative_hazard
#' 
#' @param t A vector of times to event.
#' @param interval_bounds Interval boundaries for the hazard distribution of piecewise exponential.
#' @param lambda Baseline hazard parameters of the piecewise exponential distribution.
#' 
#' @return A list with \code{h0} (hazard function values) and \code{H0} (cumulative hazard function values).
#' 
#' @export
cumulative_hazard <- function(t, interval_bounds, lambda) {
  
  n                 <- length(t)              # Number of subjects
  J                 <- length(lambda)         # Number of intervals
  
  cumulative_hazard <- c()
  hazard            <- c()
  
  for (i in 1:n) {
    cumulative_hazard[i]  <- 0
    for (j in 1:J) {
      if (t[i] > interval_bounds[j] && t[i] <= interval_bounds[j+1]) {
        # Partial contribution for the current interval
        cumulative_hazard[i] <- cumulative_hazard[i] + lambda[j] * (t[i] - interval_bounds[j])
        current_interval     <- j
        break
      } else if (t[i] > interval_bounds[j+1]) {
        # Full contribution for earlier intervals
        cumulative_hazard[i] <- cumulative_hazard[i] + lambda[j] * (interval_bounds[j+1] - interval_bounds[j])
      }
    }
    hazard[i] <- lambda[current_interval]
  }
  return(list( h0 = hazard, H0 = cumulative_hazard ) )
}