\name{check.counts.format} 
\alias{check.counts.format} 
\title{modelling-age-at-death-with-Gamma::check.counts.format} 
\usage{check.counts.format(counts)}
\arguments{
  \item{counts}{A one-row data.frame comprising integer counts, with column names in capital letters}
}
\value{No value is returned}
\description{Helper function to check counts is in the required format}
\details{Performs a few basic checks to ensure the format of the raw counts data is in the correct format. 	Required format is a one-row data.frame comprising column names of age classes (or multi-classes) in capital letters from 'A'. For example, Payne's 1973 caprine age classes are 'A','B','C','D','E','F','G','H',and 'I'. However, it is more likely that the raw counts will include some multi-class assignments, for example 'ABCD','A','AB','B','BC' etc. Order does not matter. Row name is unconstrained, but typically specifies the name of the data, e.g. the site/phase or any other ID. Outputs warnings of any formatting problems}
\examples{
# load the raw data counts for all ten Neolithic sites
load(data)

# select the 'TRA1' subset
counts <- data['TRA1',]

# eyeball data
print(counts)

# check format
check.counts.format(counts)

# introduce an error
counts.problem <- counts
counts[4] <- 3.2
check.counts.format(counts)
}
