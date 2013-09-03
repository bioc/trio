\name{ped2geno}
\alias{ped2geno}

\title{
Transformation of Ped-File
}
\description{
Transforms a ped-file into a genotype file as required by, e.g., the functions for computing the genotypic TDT.
}
\usage{
ped2geno(ped, snpnames = NULL, coded = c("12", "AB", "ATCG", "1234"), 
   naVal = 0, cols4ID = FALSE)
}

\arguments{
  \item{ped}{a data frame in ped format, i.e. the first six columns must contain information on the families as typically
    presenteed in ped files, where
    the column names of these six columns must be "famid", "pid", "fatid", "motid", "sex","affected". The last two
    of these six columns are ignored. The IDs of individuals in the second column must be unique (not only within the family,
    but among all individuals). The columns following the six columns are assumed to contain the alleles of the SNPs, where
    the alleles are coded using the letters/numbers in \code{coded}, and missing values are coded by \code{naVal}. Thus, the
    seventh and the eigth column contain the two alleles for the first SNP, the ninth and tenth the two alleles for the second
    SNP, and so on. Contrary to the names of the first six columns, the names of the columns representing the SNPs are ignored,
    and SNP names can be specified using \code{snpnames}.
}
  \item{snpnames}{a character vector containing the names of the SNPs. If not specified, generic names are assigned (i.e. 
    \code{SNP1}, \code{SNP2}, ...). Ignored if \code{ped} just contains one SNPs.
}
  \item{coded}{the coding used for the alleles of the SNPs. \code{coded = "12"}, e.g., means that one of the alleles is coded
    by 1, and the other by 0. \code{coded = "ATCG"} means that the alleles are coded by the actual base.
}
  \item{naVal}{the value used for specifying missing values.
}
  \item{cols4ID}{logical indicating whether columns should be added to output matrix containing the family ID and the individual ID.
    If \code{FALSE}, the individual IDs are used as the row names of the output matrix.
}
}

\value{A vector (if \code{ped} consists of alleles for one SNP) or matrix (otherwise) 
containing one column for each SNP representing the genotypes of the respective SNP, where the genotypes are coded
by 0, 1, 2 (i.e. the number of minor alleles), and missing values are represented by \code{NA}. The vector or matrix contains
\eqn{3 * t} values for each SNP genotyped at the \eqn{t} trios, where each block of 3 values is composed of the genotypes of the father,
the mother, and the offspring (in this order) of a specific trio. If data for a family with more than one children are available,
each of the children is treated as a separate trio. 
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{tdt}}, \code{\link{tdt2way}}, \code{\link{trio.check}}
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


\keyword{manip}

