\name{findLDblocks}
\alias{findLDblocks}
\alias{print.LDblock}
\alias{splitBlocks}

\title{
Identifying LD blocks
}
\description{
Finds LD blocks using the procedure proposed by Gabriel et al. (2002).
}
\usage{
findLDblocks(x, alpha = 0.1, ciLD = c(0.7, 0.98), cuRecomb = 0.9, 
   ratio = 9, alsoOthers = FALSE, parentsOnly = FALSE, iter = 50,
   snp.in.col = TRUE)
   
splitBlocks(blocks)
}

\arguments{
  \item{x}{either the output of \code{\link{getLD}} or \code{\link{getLDlarge}}, respectively,
     or a numeric matrix
     consisting of the integers 0, 1, and 2, where these integers are
     assumed to be the number of minor alleles that the respective
     SNPs shows at the respective subject. Missing values are allowed. By default,
     each column of this matrix represents a SNP, and each row a subject
     (for details, see \code{snp.in.col}). The SNPs must be ordered by their
     position on the considered chromosome.
}
  \item{alpha}{numeric value between 0 and 1. For each pair of SNPs, a two-sided 
     100 * (1 - \code{alpha})\% confidence interval of D' is computed, and used
     to specify pairs of SNPs that are either in strong LD, or show historical
     evidence of recombination (see \code{ciLD} and \code{cuRecomb}). All SNP
     pairs not falling into these two categories are specified as 'Others'.  
}
  \item{ciLD}{numeric vector consisting of two values between 0 and 1. If the 
     lower bound of the confidence interval of D' for a SNP pair is larger than
     or equal to the first value in \code{ciLD} and the upper bound is larger
     than or equal to the second value, then this pair of SNP is considered to
     be in strong LD.
}
  \item{cuRecomb}{numeric value between 0 and 1. If the upper bound of the confidence
     interval of D' for a SNP pair is smaller than \code{cuRecomb}, then this pair
     of SNP is considered to show evidence of recombination.
}
  \item{ratio}{numeric value larger than 1. If in a block of SNPs, the ratio of the
     number of SNP pairs being in strong LD to the number of SNPs showing evidence
     of recombination is larger than or equal to \code{ratio}, then this block will
     be identified as an LD-block. (Note that Gabriel et al. (2002) use \code{ratio} = 19
     instead of \code{ratio} = 9.) Overlapping blocks are avoided by employing
     the approach described in Wall and Pritchard (2003).
}
  \item{alsoOthers}{logical value. Following the description of Wall and Pritchard (2003)
     the endmarkers of a LD block must be in strong LD. By default (i.e.\ if \code{alsoOthers = FALSE}),
     this condition is used. If \code{alsoOthers = TRUE}, the endmarkers can also be categorized
     as 'Others'.
}
  \item{parentsOnly}{logical indicating whether only the genotypes of the parents, i.e.\ rows 1, 2,
    4, 5, ... of \code{x}, should be used in the computation of the LD measures when \code{x}
    is in genotype format and contains case-parent trio data (see \code{\link{ped2geno}} and \code{\link{read.pedfile}}).
    If \code{FALSE} (default), all rows are used in the determination of the pairwise LD measure. Ignored
    if \code{x} is the output of \code{getLD} or \code{getLDlarge}.
} 
  \item{iter}{integer specifying the number of iterations used in the computation of D (for details,
      see \code{\link{getLD}}). Ignored if \code{x} is the output of \code{getLD}.
}
  \item{snp.in.col}{logical specifying whether each column of \code{x} represents a SNP (and each row a subject).
     If \code{FALSE}, each row represents a SNP (and each column a subject). Ignored if \code{x} is the
     output of \code{getLD} or \code{getLDlarge}.
}     
  \item{blocks}{output of \code{findLDblocks}. See Details.
}

}
\details{
   The LD-blocks are estimated using the method of Gabriel et al. (2002) as described in Wall and Pritchard (2003),
   where we use the approximate variance estimates of D' proposed by Zabaleta et al. (1997).
   
   Since in \code{\link{trio.prepare}} the LD blocks are restricted to a maximum of 7 SNPs, \code{splitBlocks}
   can be used to split LD blocks composed of more than 7 SNPs into smaller blocks, if the output of \code{findLDblocks}
   should be used in \code{\link{trio.prepare}} to prepare a matrix for a \code{\link{trioLR}} or \code{\link{trioFS}}
   analysis.  
}
\value{An object of class \code{LDblocks} consisting of 
    \item{\code{ld}}{the output of \code{getLD},}
    \item{\code{blocks}}{a vector specifying which SNP belongs to which LD-block,}
    \item{\code{vec.blocks}}{a list in which each entry contains the names of the SNPs belonging to a specific LD-block,}
    \item{\code{param}}{a list of the input parameters.}
}
\references{
  Gabriel, S.B. et al. (2002). The Structure of Haplotype Blocks in the Human Genome. \emph{Science}, 296, 2225-2229.
  
  Wall, J.D. and Pritchard J.K. (2003). Assessing the Performance of the Haplotype Block Model of Linkage Disequilibrium. 
    \emph{American Journal of Human Genetics}, 73, 502-515.
    
  Zapata, C., Alvarez, G., and Carollo, C. (1997). Approximate Variance of the Standardized Measure of Gametic
     Disequilibrium D'. \emph{American Journal of Human Genetics}, 61, 771-774.
}
\author{
  Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
   \code{\link{plot.LDblocks}}, \code{\link{getLD}}
}

\examples{
# Load the simulated data.
data(trio.data)

# Estimate LD blocks.
blocks <- findLDblocks(LDdata)

# Alternatively, the LD blocks can be estimated by
ld.out <- getLD(LDdata, addVarN=TRUE)
blocks2 <- findLDblocks(ld.out)
}

\keyword{multivariate}
