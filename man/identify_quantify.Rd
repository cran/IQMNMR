\name{identify_quantify}
\alias{identify_quantify}
\title{
Identification And Quantification Of Metabolites
}
\description{
This function is used to identify and quantify metabolites from 1D 1H NMR FID.  Before using this function, users must give two data which are names as space_x and lists_metabolites.  These two data's paradigm can be obtained by type data(space_x) and data(lists_metabolites).  These two data describe the character of targeted metabolites.  These information can be obtained from http://www.bmrb.wisc.edu/ and http://mmcd.nmrfam.wisc.edu/.  
}
\usage{
identify_quantify(TSP_or_DSS, TSP_or_DSS_concentration, variance_frequency, variance_frequency_interval, TD, AQ, SWH, SFO1, O1, SNR, convergence_index, maximal_vector_angle, DE, data_x, data_y, FIDL, FID.file, threads)
}

\arguments{
  \item{TSP_or_DSS}{
The internal standard used in sample. 1 represents TSP, 2 represents DSS.
}
  \item{TSP_or_DSS_concentration}{
The concentration of internal standard used in sample. The unit is micromoles/L or moles/L.
}
  \item{variance_frequency}{
The resonance frequencies are not all the same as the reference frequencies, even for the same condition. This argument defines the variance of the resonance frequency.
}
  \item{variance_frequency_interval}{
An important characteristic parameter in NMR is the frequency interval of the resonance frequencies. This argument defines the variance of the frequency intervals.
}
  \item{TD}{
Time domain data size of sample. It is an NMR experiment parameter. Please refer to BRUKER NMR Guide.
}
  \item{AQ}{
The acquisition time of sample. It is an NMR experiment parameter. The unit is second. Please refer to BRUKER NMR Guide.
}
  \item{SWH}{
The spectral width of sample. The unit is Hz. It is an NMR experiment parameter. Please refer to BRUKER NMR Guide.
}
  \item{SFO1}{
The transmitter frequency. The unit is MHz. It is an NMR experiment parameter. Please refer to BRUKER NMR Guide.
}
  \item{O1}{
The transmitter frequency offset in Hz. It is an NMR experiment parameter. this value can be obtained from NMR spectra.
}
  \item{SNR}{
The signal to noise Ratio.  Note that this is peak SNR.  The iterations will stop at its value.
}
  \item{convergence_index}{
This argument defines the minimal relative change of the cost function in between two consecutive iterations.  The smaller the convergence_index, the more accurate and precise the result, but the more the required time.  The default value is 0.001. 
}
  \item{maximal_vector_angle}{
Vector angle is the angle between the vector of possible signals combination and the vector of standard references' signals.  A smaller vector angle indicates a high probability that there exist the targed metabolites.  'maximal_vector_angle' is the reasonable maximal vector angle that the concentration of targetd metabolite would be thought to be zero if the vector angles between all possible sinals conbinations and this metabolites' signals are more than the 'maximal_vector_angle'.  This parameter's value should be determined by Preliminary experiment.  For example,  setting it at pi/2; using two or more metabolites to generate FIDs. after processing FIDs, the value of function "identify_quantify" will give the vector angles. User can, then, determine the reasonable maximal vector angle by comparing the measured concentration to the real concentrations.
}
  \item{DE}{
The prescan delay. It is an NMR experiment parameter. The unit is second. Please refer to BRUKER NMR Guide.
}
  \item{data_x}{
This argument provides the information about the target metabolites' chemical shifts and the intramolecular resonance signals' intensity ratios.  It must be specified by users or the code "data(space_x)".
}
  \item{data_y}{
This argument provides the information about the target metabolites' molecular constitution and the experimental conditions of reference NMR spectra.  It must be specified by users or the code "data(lists_metabolites)".
}
  \item{FIDL}{
This argument defines the number of FID points used for IQMNMR. Commonly, this value is 6000 in our experiments, because the signal intensity of FIDs had attenuated after this point.
}
  \item{FID.file}{
This argument is the file path of FID.
}
  \item{threads}{
This parameter is the unit of processing that can be scheduled by an operating system.  Its value should be assigned according to the clusters of multicore multithreaded computer or multiple computers.
}
}
\details{
This function is used to identify and quantify targeted metabolites from 1D 1H NMR FID.  Before using this function, users must give two data which describe the character of targeted metabolites.  Two paradigm can be obtained by type data(space_x) and data(lists_metabolites).

The parameters (SWH, SFO1, TD, AQ, DE, O1) are the parameters of BRUKER NMR.

The FIDL, and SNR determine the performance of this function.  Commonly, the FIDL value is 6000 in our experiments, because the NMR FID had attenuated after this length.  the smaller the value of SNR, the more sensitive the algorithm.  the small the value of SNR, the more the time required for computing.
}
\value{
The result is a dataframe which list the names, concentrations and vector angles of targeted metabolites. the vector angle is the angle between the vector of possible signals and the vector of standard references' signals.  A smaller vector angle indicates a high probability that there exist the targed metabolites.  In practise, the vector angle is affected by experimental error, the threshold of vector angle should be predetermined by experiment of standard samples.  In our simulated experiment, its value is 1.02.
}
\references{
Cui, Q., Lewis, I. A., Hegeman, A. D., Anderson, M. E., Li, J., Schulte, C. F., Westler, W. M., Eghbalnia, H. R., Sussman, M. R. and Markley, J. L. (2008) \emph{Nature biotechnology}. \bold{26:2}, 162--164.

Bi, Z., Bruner, A. P., Li, J., Scott, K. N., Liu, Z. S., Stopka, C. B., Kim, H. W., Wilson, D. C. (1999) \emph{Journal of Magnetic Resonance}. \bold{140:1}, 108--119.
}
\author{
Xu Song. China. 
}
\note{
This package is for MAC platforms.
}


\keyword{ function }

