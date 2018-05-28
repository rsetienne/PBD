#' Checked and safer version of PBD::pbd_sim
#' @param erg extinction rate of a good species
#' @param eri extinction rate of an incipient species
#' @param scr speciation completion rate
#' @param sirg speciation initiation rate of a good species
#' @param siri speciation initiation rate of an incipient species
#' @param stem_age stem age. Set either the stem age or the crown age.
#' @param crown_age crown age. Set either the stem age or the crown age.
#' @param max_n_taxa maximum number of taxa. If this value is exceeded,
#'   the simulation is aborted and removed.
#' @param add_shortest_and_longest Gives the output of the new samplemethods
#'   'shortest' and 'longest'.
#' @export
pbd_sim_checked <- function(
  erg,
  eri,
  scr,
  sirg,
  siri,
  stem_age = NULL,
  crown_age = NULL,
  max_n_taxa = Inf,
  add_shortest_and_longest = FALSE
) {
  if (erg < 0.0) {
    stop("'erg' must be positive")
  }
  if (eri < 0.0) {
    stop("'eri' must be positive")
  }
  if (sirg < 0.0) {
    stop("'sirg' must be positive")
  }
  if (siri < 0.0) {
    stop("'siri' must be positive")
  }
  if (scr < 0.0) {
    stop("'scr' must be positive")
  }
  if (is.null(crown_age)) {
    if (is.null(stem_age)) {
      stop(
        "At least one of 'crown_age' or 'stem_age' ",
        "must be non-zero and positive"
      )
    } else if (stem_age <= 0.0) {
      stop("'stem_age' must be non-zero and positive")
    }
  } else if (crown_age <= 0.0) {
    stop("'crown_age' must be non-zero and positive")
  } else if (crown_age > 0.0) {
    if (!is.null(stem_age)) {
      stop("Must set either 'crown_age' or 'stem_age'")
    }
  }
  if (max_n_taxa < 0.0) {
    stop("'max_n_taxa' must be positive")
  }

  # age
  age <- stem_age
  if (is.null(age)) age <- crown_age
  testit::assert(!is.null(age))
  testit::assert(age > 0.0)

  # soc
  soc <- 1
  if (!is.null(crown_age)) soc <- 2
  testit::assert(!is.null(soc))
  testit::assert(soc == 1 || soc == 2)

  PBD::pbd_sim(
    pars = c(sirg, scr, siri, erg, eri),
    age = age,
    soc = soc,
    plotit = FALSE,
    limitsize = max_n_taxa,
    add_shortest_and_longest = add_shortest_and_longest
  )
}
