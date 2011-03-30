\name{lists_metabolites}
\alias{lists_metabolites}
\docType{data}
\title{
The NMR-related Information Of Targeted Metabolites
}
\description{
Before identifying and quantifying metabolites,  the prior knowledge data set of targeted metabolites must be constructed.  The prior knowledge data set consists of two tables: lists_metabolites and space_x.  The table "lists_metabolites" provides the information about the targeted metabolites' molecular constitution and NMR experimental conditions.  These prior knowledge data were collected by Xu Song from Madison Metabolomics Consortium Database (MMCD).  User can also collect their data, and define their table according to this paradigm.
}
\usage{data(lists_metabolites)}
\format{
  A data frame with 9 observations on the following 10 variables.
  \describe{
    \item{\code{ID}}{The identification number of metabolites. Every metabolite has a unique ID.}
    \item{\code{name}}{The name of metabolites.}
    \item{\code{Reference}}{The name of internal standard for the quantification of metabolites.}
    \item{\code{pH}}{The pH values of standard sample.}
    \item{\code{Solvent}}{The solvent of standard sample.}
    \item{\code{Field_strength}}{The magnetic field strength of NMR spectrometer. The units are MHz.}
    \item{\code{Pulse_sequence}}{The pulse sequences used to generate reference NMR spectra.}
    \item{\code{effective_hydrogen_atom}}{The total number of hydrogen nuclei generating the metabolite's resonance lines. It is directly proportional to the sum of integrated signal areas of a metabolite. This number is named as "effective_hydrogen_atom".}
    \item{\code{molecular_formula}}{The molecular formula of metabolites.}
  }
}
\details{
Before identifying and quantifying metabolites, users must fill in tables named as space_x and lists_metabolites.  This table give a paradigm of the lists_metabolites table. The lists_metabolites table contain nine variables. 
effective_hydrogen_atom: The total number of hydrogen nuclei generating the metabolite's resonance lines. It is directly proportional to the sum of integrated signal areas of a metabolite. This number is named as "effective_hydrogen_atom". This is the most important variable which decide the quantification of metabolites.  Many metabolites' effective hydrogen atom number can obtained from web database \url{http://mmcd.nmrfam.wisc.edu/} \url{http://www.bmrb.wisc.edu/metabolomics/}.  Some metabolites should, However, be given by users.

}
\source{
\url{http://mmcd.nmrfam.wisc.edu/}

\url{http://www.bmrb.wisc.edu/metabolomics/}
}
\references{
Cui, Q., Lewis, I. A., Hegeman, A. D., Anderson, M. E., Li, J., Schulte, C. F., Westler, W. M., Eghbalnia, H. R., Sussman, M. R. and Markley, J. L. (2008) \emph{Nature biotechnology}. \bold{26:2}, 162--164.
}
\examples{
data(lists_metabolites)
}
\keyword{datasets}
