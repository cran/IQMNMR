\name{IQMNMR}
\alias{IQMNMR}
\docType{package}
\title{
Identification And Quantification Of targeted Metaboloties By Using 1D 1H NMR FID
}
\description{
This package present a new automated technique based on RELAX algorithm, didital filter,and free induction decay (FID), which is detected in one-dimensional (1D) proton (1H) NMR experiment, for the unambiguous identification and quantification of targeted metabolites.
}
\details{
\tabular{ll}{
Package: \tab IQMNMR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-08-20\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This function is used to identify and quantify targetd metabolites from 1D 1H NMR FID.  Before using this function, users must give two input data which describe the character of expected metabolites.  The paradigm of these two input data can be obtained by using the code data(space_x) and data(lists_metabolites).
In order to validate the algorithm, the simulated NMR data are used firstly.  The function select_metabolites is designed to select metabolites which generate the simulated NMR signals.  The function NMR_experiment is designed to simulate NMR experiment.  The function NMR_spectra is designed to show the simulated NMR spectra.
The function idntification_quantification is used to identify and quantify targeted metabolites from 1D 1H NMR FID.

}
\references{
Cui, Q., Lewis, I. A., Hegeman, A. D., Anderson, M. E., Li, J., Schulte, C. F., Westler, W. M., Eghbalnia, H. R., Sussman, M. R. and Markley, J. L. (2008) \emph{Nature biotechnology}. \bold{26:2}, 162--164.

Bi, Z., Bruner, A. P., Li, J., Scott, K. N., Liu, Z. S., Stopka, C. B., Kim, H. W., Wilson, D. C. (1999) \emph{Journal of Magnetic Resonance}. \bold{140:1}, 108--119.
}
\author{
Maintainer: Xu Song <songxu@zzu.edu.cn>
}

\keyword{ package }
