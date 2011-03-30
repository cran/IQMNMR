\name{space_x}
\alias{space_x}
\docType{data}
\title{
The Chemical Shift And Amplitude Ratio Of Targeted Metabolites
}
\description{
Before identifying and quantifying metabolites,  the prior knowledge data set of target metabolites must be constructed.  The prior knowledge data set consists of two tables: lists_metabolites and  space_x.  The table "space_x" provides the information about the target metabolites' chemical shifts and the intramolecular resonance signals' intensity ratios.  These prior knowledge data were collected by Xu Song from Madison Metabolomics Consortium Database (MMCD).  User can also collect their data, and define their table according to this paradigm.
}
\usage{data(space_x)}
\format{
  A data frame with 248 observations on the following 3 variables.
  \describe{
    \item{\code{ID}}{The identification number of targeted metabolites.  Every metabolite has a unique ID.}
    \item{\code{Chemical_Shifts}}{The reference chemical shifts of targeted metabolites (relative to DSS).}
    \item{\code{amplitude_ratio}}{The reference intensity ratios of intramolecular resonance signals.  The intensity of the strongest signal is set at 1.}
  }
}
\details{
Before identifying and quantifying metabolites, users must fill in two input tables named as space_x and lists_metabolites.  This table give a paradigm of the table "space_x".  The table "space_x" contain three variables. 
These values can be derived from MMCD \url{http://mmcd.nmrfam.wisc.edu/} \url{http://www.bmrb.wisc.edu/metabolomics/}. 
}
\source{
\url{http://mmcd.nmrfam.wisc.edu/}

\url{http://www.bmrb.wisc.edu/metabolomics/}
}
\references{
Cui, Q., Lewis, I. A., Hegeman, A. D., Anderson, M. E., Li, J., Schulte, C. F., Westler, W. M., Eghbalnia, H. R., Sussman, M. R. and Markley, J. L. (2008) \emph{Nature biotechnology}. \bold{26:2}, 162--164.
}
\examples{
data(space_x)
}
\keyword{datasets}
