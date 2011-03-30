\name{select_metabolites}
\alias{select_metabolites}
\title{
Select Metabolites
}
\description{
This function and the NMR_experiment function are designed to select metabolites and simulate NMR experiment.
}
\usage{
select_metabolites(n_metabolites, s.lists_metabolites)
}

\arguments{
  \item{n_metabolites}{
This argument defines the number of selected metabolites.
}
  \item{s.lists_metabolites}{
This argument defines the space of possible metabolites, i.e. prior knowledge data set of target metabolites. The selected metabolites are sampled from this data set.
}
}
\details{
In order to validate algorithm, the simulated NMR data are used firstly.  This function is designed to select metabolites which generate the simulated NMR signals.
}
\value{
The value give the information on selected metabolites and their concentration. 
}

\author{
Xu Song 
}

\examples{

data(lists_metabolites)

select_metabolites(n_metabolites=10,
                             s.lists_metabolites=lists_metabolites
                             )

}

\keyword{ function }
