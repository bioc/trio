\name{trio.data}
\docType{data}
\alias{trio.data}
\alias{trio.gen1}
\alias{trio.gen2}
\alias{trio.gen.err}
\alias{trio.ped1}
\alias{trio.ped2}
\alias{trio.ped.err}
\alias{freq.hap}
\alias{simuBkMap}
\alias{step3way}
\alias{LDdata}
\alias{mat.test}
\alias{prob.mat.test}
\title{Case-Parent Trio Data}
\description{
\code{trio.data} contains several simulated data sets used in the different examples for the analyses with the functions in the \code{R}
package \code{trio}.

For the applications of genotypic TDTs for individual SNPs and two-way interactions with, for example,
\code{\link{tdt}} and \code{\link{tdt2way}}, respectively, \code{trio.data} contains a 300 x 6 matrix
called \code{mat.test} consisting of genotype data for 100 trios genotyped at 6 SNPs.

For the application of \code{\link{probTDT}} to genotype probabilities, \code{trio.data} contains a 334 x 180 matrix called
\code{prob.mat.test} containing genotype probabilities for 334 SNPs and 20 trios.


For the preparation of the trio data for an application of trio logic regression with \code{\link{trio.check}} and
\code{\link{trio.prepare}}, \code{trio.data} contains different data set containing genotype data for 10 SNPs in 100
trios in different formats. 

\code{trio.gen1}, \code{trio.gen2}, and \code{trio.gen.err} consist of 12 columns and 300 rows, where 
the first two columns contain family identifier and individual identifier. In the columns afterwards, 
each SNPs is encoded in one variable denoting the number of minor alleles.

\code{trio.ped1}, \code{trio.ped2}, and \code{trio.ped.err} consist of 26 columns and 300 rows, where
the first six columns identify the family structure of the data, and the phenotype. Besides the variables
providing information on the family structure and the phenotypes
(columns 1 to 6), each SNPs is encoded in two variables denoting the alleles.

Contrary to the other data sets, \code{trio.gen.err} and \code{trio.ped.err} contain Mendelian errors. 

For the application of the functions \code{\link{getLD}} and \code{\link{findLDblocks}}
for computing the pairwise LD values and for detecting the LD blocks, respectively,
\code{trio.data} contains a 500 x 50 matrix called \code{LDblock} that is composed of genotype data for 
10 LD blocks each consisting of 5 SNPs in strong LD.

Finally, for the simulation of trio data with \code{\link{trio.sim}},
\code{trio.data} contains examples for haplotype frequencies used in these
simulations. Both \code{freq.hap} and \code{simuBkMap} are \code{data.frame}s containing haplotype information, including the
haplotype block identifier, haplotype, and haplotype frequency. While \code{freq.hap} is a data frame consisting of
20 rows and 3 columns, \code{simuBkMap} consists of 66 rows and 3 columns. \code{step3way} is a list internally used
for simulation, containing some indexes and sampling frequencies.  
}

\author{\code{LDdata} and \code{mat.test}: Holger Schwender, \email{holger.schwender@udo.edu};
\code{prob.mat.test}: Margaret Taub, \email{mtaub@jhsph.edu};
all other data sets: Qing Li, \email{mail2qing@yahoo.com}}

\examples{
# Data can be loaded by
data(trio.data)
}

\keyword{datasets}
