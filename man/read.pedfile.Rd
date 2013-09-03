\name{read.pedfile}
\alias{read.pedfile}

\title{
Reading a Ped File
}
\description{
Reads a ped file into R and creates a data frame in ped format, or transform the ped file into a matrix in
genotype format.
}
\usage{
read.pedfile(file, first.row = NA, coded = NULL, naVal = 0, sep = " ", 
   p2g = FALSE, non.rs.IDs = FALSE, cols4ID=FALSE)
}

\arguments{
  \item{file}{the filename (if necessary with path) of a ped file that should be read into R.
}
  \item{first.row}{logical indicating whether the first row of \code{file} also contains data for a subject.
    If \code{FALSE}, the first row is assumed to contain the SNP names. By default, \code{read.pedfile} tries
    to figure out automatically if the first column contains the SNP names or data for a subject.
}
  \item{coded}{a character string stating how the alleles of the SNPs are coded. Possible values are \code{"12"},
  	\code{"AB"}, \code{"1234"}, \code{"ATCG"}. For details, see \code{\link{ped2geno}}. By default, \code{read.pedfile}
  	tries to figure out automatically how the alleles are coded. 
}
  \item{naVal}{value or character string specifying how missing values in the SNP data are coded.
}
  \item{sep}{character string specifying how the SNP names in the first row of \code{file} are separated. Ignored
    if \code{first.row = TRUE}.
}
  \item{p2g}{logical indicating whether the ped file should be transformed into a matrix in genotype format. If \code{FALSE},
    a data frame in ped format is returned. Otherwise, \code{\link{ped2geno}} is called within \code{read.pedfile} to transform
    the data frame into a matrix in genotype format, and the matrix is returned.
}
  \item{non.rs.IDs}{logical indicating whether (some of) the SNP names are specified by other names than rs-IDs.
}
  \item{cols4ID}{logical indicating whether columns should be added to output matrix containing the family ID and the individual ID.
    If \code{FALSE}, the individual IDs are used as the row names of the output matrix.
}
}

\value{A data frame in ped format (if \code{p2g = FALSE}), or a matrix in genotype format (if \code{p2g = TRUE}).
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{ped2geno}}
}

\examples{\dontrun{
# Assuming there is a ped-file called pedfile.ped in the 
# R working directory, this file can be read into R by
ped <- read.pedfile("pedfile.ped")

# The resulting data frame is in the typical ped format
# which needs to be transformed into the genotype format
# for applications of most of the functions in the trio
# package. This transformation can be done by
geno <- ped2geno(ped)

# This transformation can also be done directly when
# reading the ped-file into R by
geno2 <- read.pedfile("pedfile.ped", p2g = TRUE)
}}

\keyword{IO}
\keyword{file}

