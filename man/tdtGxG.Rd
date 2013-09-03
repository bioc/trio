\name{tdtGxG}
\alias{tdtGxG}
\alias{tdt2way}
\alias{colGxG}
\alias{colTDT2way}
\alias{print.tdtEpi}
\alias{print.colTDTepi}
\alias{colTDTinter2way}
\alias{fastGxG}
\alias{fastGxGrec}

\title{
Genotypic TDT for Two-Way Interactions
}
\description{
\code{tdtGxG} and \code{colGxG} perform the genotypic TDT for the interaction of two SNPs or of each pair of columns 
of a genotype matrix, respectively.

\code{fastGxG} provides a fast implementation for the genotypic TDT for two-way interactions when considering the
simplest conditional logistic regression model only containing one parameter for the interaction effect. It thus leads
to the same results as \code{colGxG} with \code{test = "screen"}. In \code{fastGxGrec},
an analytic solution to the genotypic TDT based on the simplest model for testing a recessive x recessive model is
implemented, which is even faster than \code{fastGxG} with \code{model = "recessive"}. In future versions of this
package, \code{fastGxG} and \code{fastGxGrec} will be joint with \code{colGxG}.

The genotypic TDT for testing two-way interactions makes use of the 16 possible genotypes that can be obtained 
from combining the parents' genotypes of the two considered SNPs. Thus, for each family, genotypes for one case 
(i.e. the affected offspring) and 15 pseudo-controls are used.
}

\usage{
tdtGxG(snp1, snp2, test = c("epistatic", "lrt", "full", "screen"), 
     model = c("additive", "dominant", "recessive"))
     
colGxG(mat.snp, test = c("epistatic", "lrt", "full", "screen"), genes = NULL,
   maf = FALSE, model = c("additive", "dominant", "recessive"))
   
fastGxG(mat.snp, model = c("additive", "dominant", "recessive"), 
   genes = NULL, interval = c(-10, 10), tol = 10^-8, maxiter = 1000, 
   size = 20)
   
fastGxGrec(mat.snp, genes = NULL, size = 20)
}

\arguments{
  \item{snp1, snp2}{numeric vectors of length \eqn{3 * t} representing two SNPs genotyped at \eqn{t} trios. Each of the \eqn{t}
    blocks (i.e. \code{snp1[1:3]}, \code{snp1[4:6]}, ..., and \code{snp2[1:3]}, \code{snp2[4:6]}, ...) must consist of the 
    genotypes of father, mother, and offspring
    (in this order). The genotypes must be coded by 0, 1, and 2. Missing values are allowed and need to be coded by \code{NA}.
    The vectors must thus have the same structure as the output of \code{\link{trio.check}}, or the genotype example data
    sets such as \code{trio.gen1} (see \code{data(trio.gen1)}), and can be generated from a ped-file by, e.g., employing 
    \code{\link{ped2geno}}.
}
  \item{mat.snp}{a numeric matrix in which each column represents a SNP. Each of the SNPs must have
    the same structure as \code{snp}, and can, e.g., be generated from a ped-file by employing \code{\link{ped2geno}}.
}
  \item{test}{character string naming the GxG test that should be performed. If \code{test = "epistatic"}, then a conditional
    logistic regression version of the test proposed by Cordell (2002) is used to test for epistatistical interactions. If
    \code{test = "full"}, a conditional logistic regression model containing one parameter for each SNP and one parameter for
    the interaction of these two SNPs will be fitted and a Wald test for the interaction term will be performed, where a
    genetic model specified by \code{model} is assumed for both SNPs. If \code{test = "lrt"}, a likelihood ratio test is
    performed comparing the fit of this model with the fit of a conditional logistic regression model only containing the
    two parameters for the main effects of the SNPs. If \code{test = "screen"}, a conditional logistic regression model
    only composed of one parameter for the interaction of the two SNPs will be fitted and a Wald test will be performed,
    where the genetic model specified by \code{model} is assumed for both SNPs.
}
  \item{genes}{a character vector containing the names of the genes to which the SNPs belong. If specified, only the two-way interactions
    between SNPs from different genes are tested. If \code{NULL}, all two-way interactions between all possible pairs of SNPs are tested.
}
  \item{maf}{logical indicating whether the minor allele frequency (computed by considering the genotypes of only the parents) should
    be added to the output.
}
  \item{model}{type of model that should be considered. Abbreviations are allowed. Thus, e.g., \code{model = "dom"} will
    consider a dominant model for each of the respective two SNPs, and \code{model = "r"} an recessive model. Ignored if 
    \code{epistatic = TRUE}. 
}
  \item{interval}{the end-points of the interval to be searched for the root. For details, see \code{\link[stats]{uniroot}}.
}
  \item{tol}{the desired accuracy/convergence tolerance. For details, see \code{\link[stats]{uniroot}}.
}
  \item{maxiter}{the maximum number of iterations. For details, see \code{\link[stats]{uniroot}}.
}
  \item{size}{the number of interactions considered simultaneously when computing the parameter estimates. 
}
}

\value{
Depending on \code{test}, the output contains statistics and p-values either of a likelihood ratio test (\code{test = "epistatic"} or
\code{test = "lrt"}) or the Wald statistics and the corresponding p-values for the interaction term in the conditional logistic
regression model (\code{test = "full"} or \code{test = "screen"}). If \code{maf = TRUE}, a vector \code{maf}
containing the minor allele frequencies of each SNP and a matrix \code{mat.maf} with two columns containing the SNP-wise minor allele
frequencies for each tested pair of SNPs are added to the output of \code{colGxG}.
}

\references{
Cordell, H. J. (2002). Epistasis: What it Means, what it Doesn't mean, and Statistical Methods to Detect it in Humans.
Human Molecular Genetics, 11, 2463-2468.

Schwender, H., Taub, M.A., Beaty, T.H., Marazita, M.L., and Ruczinski, I. (2011). 
Rapid Testing of SNPs and Gene-Environment Interactions in Case-Parent Trio Data Based on 
Exact Analytic Parameter Estimation. \emph{Biometrics}, 68, 766-773.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link{tdt}}, \code{\link{ped2geno}}
}

\examples{
# Load the simulated data.
data(trio.data)

# The interaction between a particular pair of SNPs
# (e.g., the ones in the first and second column of 
# mat.test) can be tested by
gxg.out <- tdtGxG(mat.test[,1], mat.test[,2])

# All pairs of SNPs in mat.test can be tested by
gxg.out2 <- colGxG(mat.test)

# By default, Cordell's likelihood ratio test for
# epistatistic interactions is used. This is the
# most sophisticated, but also most time-consuming
# test. If another test, e.g., the one considering
# a conditional logistic regression model only
# containing a term for the interaction, should
# be used, then this can be done by
gxg.out3 <- colGxG(mat.test, test = "screen")

# In this case, different modes of inheritance can
# be considered (by default, the additive mode is
# considered). If a dominant model (for both SNPs)
# should be tested, this can be done by
gxg.out4 <- colGxG(mat.test, test = "screen", model ="dom")

# If just a subset of all pairs of SNPs should be 
# tested, e.g., only pairs of SNPs belonging to different
# genes, then this can be done by first specifying a
# vector specifying which SNP belongs to which genes.
# If we, e.g., assume that the first two SNPs in mat.test 
# belong to gene G1 and the other four SNPs to G2, then
# this vector can be specified by
genes <- paste("G", rep(1:2, c(2,4)), sep="")

# and only the pairs of SNPs in which the two SNPs belong
# to different genes can be tested with Cordell's 
# likelihood ratio test by
gxg.out5 <- colGxG(mat.test, genes = genes)
}


\keyword{htest}
\keyword{models}
\keyword{array}
