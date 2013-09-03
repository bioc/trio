\name{trio.check}
\alias{trio.check}

\title{Check Case-Parent Trio Data for Mendelian Errors}

\description{This function checks case-parent trio data in linkage or
genotype format for Mendelian errors. If no errors are found, the
function returns an object suitable for input to the \code{\link{trio.prepare}}
function. Otherwise, an object identifying the Mendelian errors is
returned.}

\usage{trio.check(dat, is.linkage=TRUE, replace=FALSE)}

\arguments{

  \item{dat}{A matrix or data frame of pedigree data in linkage format,
    or in genotype format.

    If the data are in \strong{linkage format}, the file has to have the
    standard linkage/pedigree format. Each row describes an individual,
    and the columns are \emph{<famid> <pid> <fatid> <motid> <sex>
    <affected> <genotype:1_1> <genotype:1_2> ... <genotype:n_1>
    <genotype:n_2>}.  Here, \emph{<famid>} is a unique identifier for
    each family, \emph{<pid>} is a unique identifier for an individual
    within each family, \emph{<fatid>} and \emph{<motid>} identify the
    father and mother of the individual, \emph{<sex>} denotes the
    gender, using the convention 1=male, 2=female, \emph{<affected>}
    denotes the disease status (0=unknown, 1=unaffected, 2=affected).
    Only one phenotype column is allowed. Each genotype is encoded using
    two columns (\emph{<genotype:k_1>} and \emph{<genotype:k_2>}),
    identifying the alleles (1 for the major allele, 2 for the minor
    allele, 0 if missing).  Other values for the alleles will result in
    an error. Please see the data frames \code{trio.ped1} and
    \code{trio.ped2} contained in this package as examples for trio data
    in linkage file format (complete and with missing records,
    respectively).
    
    If the data are in \strong{genotype format}, each row in the
    object describes an individual, and each block of three consecutive
    rows describes the two parents and the affected child in a trio. The
    columns in the object are \emph{<famid> <pid> <genotype_1>
    ... <genotype_n>}.  Here, \emph{<famid>} is a unique identifier for
    each family, \emph{<pid>} is a unique identifier for an individual
    within each family (with each block of three consecutive rows
    describing the two parents and the affected child in a trio).  Each
    \emph{<genotype>} is encoded as an integer indicating the number of
    variant alleles (e.g. 0=common homozygote, 1=heterozygote, and
    2=rare homozygote, and \code{NA}=missing genotype).  Please see the
    data frames \code{trio.gen1} and \code{trio.gen2} contained in this
    package as examples for trio data in linkage file format (complete
    and with missing records, respectively).
}
    
  \item{is.linkage}{A logical value indicating if the case parent data
    are in linkage file format (\code{TRUE}) or in genotype format
    (\code{FALSE}).}

  \item{replace}{A logical value indicating whether existing Mendelian
    errors should be replaced by missing values. For each Mendelian
    error found (for a particular trio at a particular locus), all three
    genotypes are replaced by \code{NA}, and an object suitable for
    input to the \code{\link{trio.prepare}} function is returned.}

}

\details{ The first function used from this package should always be
  \code{trio.check}.  Unless otherwise specified, this function assumes
  that the data are in linkage format, however, genotype data can also
  be accommodated.  If no Mendelian inconsistencies in the data provided
  are identified, \code{trio.check} creates an object that can be
  processed in the subsequent analysis with the \code{\link{trio.prepare}}
  function. If the data were in linkage format, the genotype information
  for each SNP will be converted into a single variable, denoting the
  number of variant alleles.

  To delineate the genotype information for the pseudo-controls in the
  subsequent analysis, the trio data must not contain any Mendelian
  errors. The function \code{trio.check} returns a warning, and an R
  object with relevant information when Mendelian errors are encountered
  in the supplied trio data.  It is the users responsibility to find the
  cause for the Mendelian errors and correct those, if
  possible. However, Mendelian inconsistencies are often due to
  genotyping errors and thus, it might not be possible to correct those
  in a very straightforward manner. In this instance, the user might
  want to encode the genotypes that cause theses Mendelian errors in
  some of the trios as missing data.  The function \code{trio.check}
  allows for this possibility, using the argument \code{replace=T}.}

\value{

  The function \code{trio.check} returns a list with the following
  elements:
  
  \item{trio}{A data frame with the genotypes of the trios, suitable for
    input to the function \code{\link{trio.prepare}}. This element will be
    \code{NULL} if Mendelian errors are detected.}

  \item{errors}{This element will be \code{NULL} if no Mendelian errors
    are detected. Otherwise, this element will be a data frame with
    five columns, indicating the Mendelian errors detected in the object
    \code{dat}. The five columns of the data frame refer to the trio 
    (\code{trio}), the family id (\code{famid}), the genotype (\code{snp}), the row numbers (\code{r}),
       and the column numbers (\code{c}).}

  \item{trio.err}{This element will be \code{NULL} if no Mendelian
    errors are detected. Otherwise, this element will be a data frame
    with the trio genotype data. If the input was a linkage file, the
    data will be converted from alleles to genotypes. If the input was a
    genotype file, this element will be identical to the input.}

}

\references{
Li, Q., Fallin, M.D., Louis, T.A., Lasseter, V.K., McGrath, J.A., Avramopoulos, D., Wolyniec, P.S., Valle, D., Liang, K.Y., Pulver, A.E., 
and Ruczinski, I. (2010). Detection of SNP-SNP Interactions in Trios of Parents with Schizophrenic Children. Genetic Epidemiology, 34, 396-406. 
}

\author{Qing Li, mail2qing@yahoo.com}

\seealso{\code{\link{trio.prepare}}}

\examples{
data(trio.data)
trio.tmp <- trio.check(dat=trio.ped1)
str(trio.tmp, max=1)
trio.tmp$trio[1:6,]

trio.tmp <- trio.check(dat=trio.ped.err)
str(trio.tmp, max=1)
trio.tmp$errors
trio.tmp$trio.err[1:3, c(1,2, 11:12)]
trio.ped.err[1:3,c(1:2, 23:26)]

trio.tmp <- trio.check(dat=trio.gen.err, is.linkage=FALSE)
trio.tmp$errors
trio.tmp$trio.err[1:6, c(1,2,7), drop=FALSE]

trio.rep <- trio.check(dat=trio.gen.err, is.linkage=FALSE, replace=TRUE)
trio.rep$trio[1:6,c(1,2,7)]
}

\keyword{manip}



  
