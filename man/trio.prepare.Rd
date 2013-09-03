\name{trio.prepare}
\alias{trio.prepare}

\title{Generate Trio Data Format Suitable for Trio Logic Regression}

\description{This function transforms case-parent data into a format
  suitable as input for trio logic regression. The function can also be
  used for the imputation of missing genotypes in case-parent data,
  while taking the existing SNP block structure into account.}

\usage{trio.prepare(trio.dat, freq=NULL, blocks=NULL, logic=TRUE, \ldots)}

\arguments{

  \item{trio.dat}{An object returned from the function
    \code{\link{trio.check}}.}

  \item{freq}{An optional data frame specifying haplotype blocks and
    frequencies. For an example, see the data frame \code{simuBkMap}
    contained in this package.  If provided, the following argument
    \code{blocks} will be ignored.
      
    The object must have three columns in the following order: block
    identifiers (\code{key}), haplotypes (\code{hap}), and haplotype
    frequencies (\code{freq}).  The block identifiers must be unique for
    each block. For each block, the haplotypes must be encoded as a
    string of the integers 1 and 2, where 1 refers to the major allele
    and 2 refers to the minor allele. The respective haplotype
    frequencies will be normalized to sum one.
 }
  \item{blocks}{An optional vector of integers, specifying (in sequence)
    the lengths of the linkage disequilibrium blocks.  The sum of these
    integers must be equal to the total numbers of SNPs in the data set
    used as input.  Using the integer 1 for SNPs not contained in LD
    blocks is required if this argument is used.  If both arguments
    \code{freq} and \code{blocks} are \code{NULL}, complete linkage
    equilibrium is assumed (i.e., no correlation between the
    genotypes).}

  \item{logic}{A logical value indicating whether the trio data are
    returned with genotypes in dominant and recessive coding, suitable
    as input for trio logic regression (\code{TRUE}), or if the imputed
    data should be returned in genotype format, using one variable per
    SNP (\code{FALSE}).}

  \item{...}{Optional arguments that can be passed to function
    \code{\link[haplo.stats]{haplo.em}}.}

}

\details{ To create the genotypes for the pseudo-controls it is
  necessary to take the LD structure of the SNPs into account.  This
  requires information on the LD blocks.  It is assumed that the user
  has already delineated the block structure according to his or her
  method of choice.  The function \code{trio.prepare}, which operates on an
  output object of \code{\link{trio.check}}, accepts the block length
  information as an argument. If this argument is not specified, a
  uniform block length of 1 (i.e., no LD structure) is assumed. If the
  haplotype frequencies are not specified, they are estimated from the
  parents' genotypes using the function \code{\link[haplo.stats]{haplo.em}}.  The
  function then returns a list that contains the genotype information in
  binary format, suitable as input for trio logic regression.  Since
  trio logic regression requires complete data, the function \code{trio.prepare}
  also performs an imputation of the missing genotypes. The imputation
  is based on the estimated or supplied haplotype information.
}

\value{

  \item{bin}{A matrix suitable as input for trio logic regression.  The
    first column specifies the cases and pseudo-controls as required by
    logic regression using conditional logistic regression (the integer
    3 for the probands followed by three zeros indicating the
    pseudo-controls). The following columns specify the (possibly
    imputed) genotypes in dominant and recessive coding, with two binary
    variables for each SNP. This is returned only if \code{logic = TRUE}.}

  \item{trio}{A data frame with imputed SNPs in genotype format derived
    from the input. This is returned only if \code{logic = FALSE}.}

  \item{miss}{A data frame with five columns indicating the missing genotypes in the input object.
     The five columns of the data frame refer to the family id (\code{famid}), the individual id (\code{pid}),
     the genotype (\code{snp}), the row numbers (\code{r}), and the column numbers (\code{c}). This element
      will be \code{NULL} if there are no missing data. }
  
  \item{freq}{The estimated or supplied haplotype information, in the
    same format as described in the \bold{Arguments} above.}

}

\references{
Li, Q., Fallin, M.D., Louis, T.A., Lasseter, V.K., McGrath, J.A., Avramopoulos, D., Wolyniec, P.S., Valle, D., Liang, K.Y., Pulver, A.E., 
and Ruczinski, I. (2010). Detection of SNP-SNP Interactions in Trios of Parents with Schizophrenic Children. Genetic Epidemiology, 34, 396-406. 
}

\section{Acknowledgments}{Support was provided by NIH grants R01 DK061662 and HL090577.}

\author{Qing Li, mail2qing@yahoo.com}

\seealso{\code{\link{trio.check}}, \code{\link[haplo.stats]{haplo.em}}}

\examples{
data(trio.data)
trio.tmp <- trio.check(dat=trio.ped1)
trio.bin <- trio.prepare(trio.dat=trio.tmp, blocks=c(1,4,2,3))
trio.bin$bin[1:8,]
}

\keyword{manip}
