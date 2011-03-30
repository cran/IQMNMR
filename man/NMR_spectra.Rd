\name{NMR_spectra}
\alias{NMR_spectra}

\title{
Show NMR Spectra 
}
\description{
This function is used to show NMR spectra of stimlated NMR experiment.
}
\usage{
NMR_spectra(SWH, AQ, TD, fidfile, SFO1, O1)
}

\arguments{
  \item{SWH}{
The spectral width of NMR experiment. The unit is Hz. It is an NMR experiment parameter. Please refer to BRUKER NMR Guide.
}
  \item{AQ}{
The acquisition time. It is an NMR experiment parameter. The unit is second. Please refer to BRUKER NMR Guide.
}
  \item{TD}{
Time domain data size. It is an NMR experiment parameter. Please refer to BRUKER NMR Guide.
}
  \item{fidfile}{
The file path used to load FID file.
}
  \item{SFO1}{
The transmitter frequency. The unit is MHz. It is an NMR experiment parameter. Please refer to BRUKER NMR Guide.
}
  \item{O1}{
The transmitter frequency offset in Hz. It is an NMR experiment parameter. this value can be obtained from NMR spectra.
}
}
\details{
This function is used to show NMR spectra of simulated NMR experiment.
}
\value{
The value returned show a plot of stimulated NMR spectra, and give the data of time domain and the data of frequency domain.
}
\author{
XU Song. China. 
}

\examples{

data(lists_metabolites)

data(space_x)

stimulated_metabolites<-select_metabolites(n_metabolites=10,
                             s.lists_metabolites=lists_metabolites
                             )



stimulated_experiment_a<-NMR_experiment(
                         t.stimulated_metabolites_1=as.data.frame(stimulated_metabolites),
                         s_data_x=space_x,
                         s_data_y=lists_metabolites,
                         variance_frequency=10,
                         variance_frequency_interval=0.5,
                         SNR= 100, 
                         TD=65536,
                         SWH=6218.905,
                         SFO1=400.1881,
                         O1=2001.6,
                         fid_save_path="fid_stimulate"
                         )



SWH_a<-stimulated_experiment_a[1,3]
AQ_a<-stimulated_experiment_a[1,1]
TD_a<-stimulated_experiment_a[1,2]
fidfile_a<-"fid_stimulate"
O1_a<-stimulated_experiment_a[1,5]
SFO1_a<-stimulated_experiment_a[1,4]

stimulated_spectra_a<-NMR_spectra(
                                  SWH=SWH_a,
                                  AQ=AQ_a, 
                                  TD=TD_a,
                                  fidfile=fidfile_a,
                                  SFO1=SFO1_a,
                                  O1=O1_a
                                  )

}

\keyword{ function }
