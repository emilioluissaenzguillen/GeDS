#' @useDynLib GeDS
#' @importFrom Rcpp  sourceCpp
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
################################### CrystalData ################################
################################################################################
################################################################################
#' @title Crystallographic Scattering Data
#' @name CrystalData
#' @aliases CrystalData10k CrystalData300k
#' @description
#' This dataset contains crystallographic measurements obtained from a particle 
#' accelerator experiment. The measurements correspond to the function \eqn{F(Q)} 
#' at various \eqn{Q} values, which are used to analyze the scattering properties 
#' of an unknown crystalline material. The dataset is available in two versions 
#' based on the precision of the measurements:
#'
#' - **`CrystalData10k`** (lower precision);
#' - **`CrystalData300k`** (higher precision).
#'
#' The goal of the experiment is to estimate \eqn{F(Q)} from noisy data using 
#' a GeDS model and compute its Fourier transform, which provides valuable insights 
#' into the structure of the material.
#'
#' @format A \code{data.frame} with 1721 observations and 2 variables:
#' \itemize{
#'   \item{\code{Q} (\eqn{\text{Å}^{-1}})}: The scattering vector, measured in inverse angstroms (\eqn{\text{Å}^{-1}}).
#'   \item{\code{FQ} (a.u.)}: The measured function \eqn{F(Q)}, given in arbitrary units (a.u.).
#' }
#' @usage data(CrystalData10k) 
#' @usage data(CrystalData300k)
#' @source Data collected from a particle accelerator experiment.
#' @examples \dontrun{
#' # Load the dataset (choose 10k or 300k version)
#' data('CrystalData10k')
#' 
#' # Fit a GeDS/GeDSboost model and compare how well the intensity peaks are captured
#' Gmod <- NGeDS(F_Q ~ f(Q), data = CrystalData10k, phi = 0.999, q = 3)
#' # for CrystalData300k set int.knots_init = 1, phi = 0.999, q = 4, instead
#' Gmodboost <- NGeDSboost(F_Q ~ f(Q), data = CrystalData10k, phi = 0.9975, q = 4) 
#' 
#' par(mfrow = c(1,2))
#' plot(Gmod, n = 2)
#' plot(Gmodboost, n = 2) 
#' }
#' 
#' @docType data
NULL


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
#' \itemize{
#'   \item{angle}: the dispersion angle, viewed as the independent variable.
#'   \item{intensity}: the neutron diffraction intensity, viewed as the response
#'   variable.
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
#' @title Coal Mining Disasters Data
#' @name coalMining
#' @description
#' A dataset with 112 entries containing annual numbers of accidents due to
#' disasters in British coal mines for years from 1850 to 1962, considered in
#' Carlin et al. (1992) and also in Eilers and Marx (1996).
#' 
#' @format A \code{data.frame} with 112 entries, corresponding to the
#' years from 1850 to 1962. Each entry has:
#' \itemize{
#'   \item{accidents}: number of severe accidents that have occurred each year.
#'   \item{years}: year during which the accidents occurred.
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
#' @title Death Counts in England and Wales
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

