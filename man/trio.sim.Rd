\name{trio.sim}
\alias{trio.sim}

\title{Simulate Case-Parent Trios}

\description{\code{trio.sim} generates case-parents trios when the disease
  risk of children is specified by (possibly higher-order) SNP-SNP
  interactions. The SNP minor allele frequencies and/or haplotypes are
  specified by the user, as are the parameters in the logistic model
  that describes the disease risk. If \code{pi.usr} is specified, a
  specific type of model, namely the well-known Risch model, will be employed.}

\usage{
trio.sim(freq, interaction = "1R and 2D", prev = 1e-3, OR = 1, pi.usr = 0, 
   n = 100, rep = 1, step.save = NULL, step.load = NULL, verbose = FALSE)
}

\arguments{

 \item{freq}{A data frame specifying haplotype blocks and
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

  \item{interaction}{ A string that specifies the risk altering genotype
    interaction as a Boolean term, such as \emph{"7D or 19R"}, or
    \emph{"(not 10D) or 45D"}.  Each locus can appear at most once in
    the string, and the the Boolean term \emph{not} can appear at most
    once before each locus, and must be enclosed in parenthesis, e.g.,
    \emph{"(not 3D)"}.  Therefore, strings such as \emph{"not (not 3D)"}
    and \emph{"not 3D or 5R"} are prohibited.  Parenthesis are also used
    to unambiguously define the Boolean expression as a binary tree,
    i.e., every parent node has exact two children.  For example Thus, a
    long string such as \emph{"1R or 3D or 5R"} must be written as
    \emph{"(1R or 3D) or 5R"} or as \emph{"1R or (3D or 5R)"}, even
    though the parenthesis are technically redundant. There is also a
    limit on the size of the interactions, please see \bold{Details}
    below.}

  \item{prev}{The prevalence of the disease in the simulated population
    among non-carriers (the "un-exposed" group).}

  \item{OR}{The odds ratio of disease in the simulated population,
    comparing carriers to non-carriers.}
    
  \item{pi.usr}{probability for an individual without the \code{interaction}
     to be affected}

  \item{n}{The number of case-parent trios simulated. The default is
    100.}

  \item{rep}{The number of data set replicates generated. The default is
    1.}

  \item{step.save}{The name of the binary file (without ".RData"
    extension) in which the object specifying the simulation mating
    tables and probabilities will be saved.  The default value is
    \code{NULL} In that case, the object will not be saved for re-use in
    later run. See \bold{Details}.}

  \item{step.load}{The name of an existing binary file (without ".RData"
    extension) in which the object specifying the simulation mating
    tables and probabilities have been saved (see above).  The default
    value is \code{NULL}. In that case, a new object will be generated.}
  
  \item{verbose}{A logical value indicating whether or not to print
    information about memory and time usage.}
}


\details{The function \code{trio.sim} simulates case-parent trio data
  when the disease risk of children is specified by (possibly
  higher-order) SNP-SNP interactions. The mating tables and the
  respective sampling probabilities depend on the haplotype frequencies
  (or SNP minor allele frequencies when the SNP does not belong to a
  block). This information is specified in the \code{freq} argument of
  the function.  The probability of disease is assumed to be described
  by the logistic term \emph{logit(p) = a + b} I\emph{[Interaction]},
  where \emph{a} = logit (\code{prev}) and \emph{b} = log(\code{OR}),
  with \code{prev} and \code{OR} specified by the user. Note that at
  this point only data for two risk groups (carriers versus
  non-carriers) can be simulated. Since the computational demands for
  generating the mating is dependent on the number of loci involved in
  the interactions and the lengths of the LD blocks that contain these
  disease loci, the interaction term can only consist of up to six loci,
  not more than one of those loci per block, and haplotype (block)
  lengths of at most 5 loci.
  
  Generating the mating tables and the respective sampling probabilities
  necessary to simulate case-parent trios can be very time consuming for
  interaction models involving three or more SNPs. In simulation
  studies, many replicates of similar data are usually required, and
  generating these sampling probabilities in each instance would be a
  large and avoidable computational burden (CPU and memory). The
  sampling probabilities depend foremost on the interaction term and the
  underlying haplotype frequencies, and as long as these remain constant
  in the simulation study, the mating table information and the sampling
  probabilities can be "recycled". This is done by storing the relevant
  information (denoted as "step-stone") as a binary R file in the
  working directory (using the argument \code{step.save}), and loading
  the binary file again in future simulations (using the argument
  \code{step.load}), speeding up the simulation process dramatically.
  It is even possible to change the parameters \code{prev} and \code{OR}
  (corresponding to \emph{a} and \emph{b} in the logistic model) in
  these additional simulations, as the sampling probabilities can be
  adjusted accordingly.
}


\value{A list of matrices, containing the simulated data sets, in
  genotype format (indicating the number of variant alleles), including
  family and subject identifiers.
}  

\references{
Li, Q., Fallin, M.D., Louis, T.A., Lasseter, V.K., McGrath, J.A., Avramopoulos, D., Wolyniec, P.S., Valle, D., Liang, K.Y., Pulver, A.E., 
and Ruczinski, I. (2010). Detection of SNP-SNP Interactions in Trios of Parents with Schizophrenic Children. Genetic Epidemiology, 34, 396-406. 
}

\author{Qing Li, mail2qing@yahoo.com}

\seealso{\code{\link{trio.prepare}}}


\examples{
data(trio.data)
sim <- trio.sim(freq=simuBkMap, interaction="1R and 5R", prev=.001, OR=2, n=20, rep=1)
sim[[1]][1:6, 1:12]
}

\keyword{datagen}
