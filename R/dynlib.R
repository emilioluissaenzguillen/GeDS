#' @useDynLib GeDS
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix rankMatrix
#' @importFrom splines splineDesign splineOrder splineKnots polySpline
#' @importFrom stats knots .lm.fit influence as.formula coef glm predict.glm qnorm delete.response gaussian glm.fit lm lm.fit lm.wfit model.frame model.matrix model.response pchisq pnorm pt qchisq qt residuals summary.lm terms coefficients
#' @importFrom methods setMethod setClass as new
#' @importFrom grDevices dev.new devAskNewPage trans3d
#' @importFrom graphics legend lines.default plot.default persp points rug
#' @importFrom Rmpfr Bernoulli
#' @importFrom utils packageDescription
NULL


.onAttach <- function(libname, pkgname) {
  #nice basic info
  vers <- packageDescription("GeDS")[["Version"]]
  packageStartupMessage("##################################################################\n",
                        "\n",
                        "This is GeDS version ", vers, ". \n",
                        "See ",
                        sQuote("package?GeDS"), " for a brief introduction.\n",
                        "Type ", sQuote("citation('GeDS')"), " to learn how to cite this package.\n",
                        "\n",
                        "Please report any issue or bug to the authors (See the description\n",
                        "file)\n",
                        "\n",
                        "##################################################################\n",
                        appendLF = TRUE)
  return(TRUE)
}


################################################################################
################################################################################
################################### BaFe2As2 ###################################
################################################################################
################################################################################
#' @title Barium-Ferrum-Arsenide Powder Diffraction Data
#' @name BaFe2As2
#' @description
#' This dataset contains the results of a neutron diffraction experiment on
#' Barium-Ferrum-Arsenide (\eqn{\mathrm{Ba Fe_2 As_2}}) powder carried out by
#' Kimber et al. (2009) and used in Kaishev et al. (2016). The neutron
#' diffraction intensity was measured at 1,151 different dispersion angles in
#' order to model the diffraction profile.
#' 
#' @format A \code{data.frame} with 1151 cases and 2 variables:
#' \describe{
#'   \item{angle}{ the dispersion angle, viewed as the independent variable.}
#'   \item{intensity}{ the neutron diffraction intensity, viewed as the response variable.}
#' }
#' @usage data(BaFe2As2)
#' @source \href{https://openaccess.city.ac.uk/12418/}{openaccess.city.ac.uk}
#' @references
#' Kimber, S.A.J., Kreyssig, A., Zhang, Y.Z., Jeschke, H.O., Valenti, R.,
#' Yokaichiya, F., Colombier, E., Yan, J., Hansen, T.C., Chatterji, T.,
#' McQueeney, R.J., Canfield, P.C., Goldman, A.I. and Argyriou, D.N. (2009).
#' Similarities between structural distortions under pressure and chemical
#' doping in superconducting \eqn{\mathrm{Ba Fe_2 As_2}}. \emph{Nat Mater},
#' \strong{8}, 471--475.
#'
#' Kaishev, V.K., Dimitrova, D.S., Haberman, S. and Verrall, R.J. (2016).
#' Geometrically designed, variable knot regression splines.
#' \emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
#' DOI: \doi{10.1007/s00180-015-0621-7}
#' @examples \dontrun{
#' # to load the data
#' data('BaFe2As2')
#'
#' # fit a GeDS regression and produce a simple plot of the result. See ?NGeDS
#' # c.f. Kaishev et al. (2016), section 4.2
#' (Gmod <- NGeDS(intensity ~ f(angle), data = BaFe2As2, beta = 0.6, phi = 0.99,
#'                q = 3, show.iters = T))
#' plot(Gmod)
#' }
#' 
#' @docType data
NULL


################################################################################
################################################################################
################################## coalMining ##################################
################################################################################
################################################################################
#' @title Coal Mining Disasters data
#' @name coalMining
#' @description
#' A dataset with 112 entries containing annual numbers of accidents due to
#' disasters in British coal mines for years from 1850 to 1962, considered in
#' Carlin et al. (1992) and also in Eilers and Marx (1996).
#' 
#' @format A \code{data.frame} with 112 entries, corresponding to the
#' years from 1850 to 1962. Each entry has:
#' \describe{
#'   \item{accidents}{ number of severe accidents that have occurred each year.}
#'   \item{years}{ year during which accidents occurred.}
#' }
#' @usage data(coalMining)
#' @source \url{https://people.reed.edu/~jones/141/Coal.html}
#' @references
#' Carlin, B.P., Gelfand, A.E. and Smith, A.F.M. (1992).
#' Hierarchical Bayesian analysis of changepoint problems.
#' \emph{Applied Statistics}, \strong{41}(2), 389--405.
#'
#' Eilers, P.H.C. and Marx, B.D. (1996). Flexible Smoothing with B-splines
#' and Penalties. \emph{Statistical Science}, \strong{11}(2), 89--121.
#' 
#' @docType data
NULL

################################################################################
################################################################################
################################# EWmortality ##################################
################################################################################
################################################################################
#' @title Death counts in England and Wales
#' @name EWmortality
#' @description
#' The dataset consists of information about the mortality of the English and
#' Welsh male population aggregated over the years 2000, 2001 and 2002.
#'
#' @format A \code{data.frame} with 109 entries and 3 variables: \code{Age},
#' \code{Deaths} and \code{Exposure}. \code{Exposure} is a mid-year estimate of
#' the population exposed to risk.
#' @usage data(EWmortality)
#'
#' @docType data
NULL

