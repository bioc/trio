\name{plot.LDblocks}
\alias{plot.LDblocks}
\title{
 Plotting a LDblock Object
}
\description{
 Plots either the pairwise D' values or the pairwise LD categorization used in the procedure of Gabriel et al. (2002).
 Additionally, the LD blocks are marked in this plot.
}

\usage{
\method{plot}{LDblocks}(x, y = "gabriel", col = NULL, start = 1, end = NA, xlab = "",
   ylab = "", cexAxis = 0.8, block.col = 2, block.lwd = 3, ...)
}
\arguments{
  \item{x}{the output of \code{findLDblocks}.
}
  \item{y}{either \code{"Dprime"} or \code{"gabriel"} (default) specifying the LD values that should be plotted.
}
  \item{col}{a vector specifying the colors used in plotting of the LD values. If \code{y = "Dprime"},
    different levels of gray will be used by default (the darker, the higher is the LD value).
    If \code{y = "gabriel"}, strong LD is by default marked by blue fields, evidence of recombination by white color,
    and others by yellow.
}
  \item{start}{integer or character string specifying the index or name of the first SNP, respectively,
    that should be plotted, where the index corresponds to the column (or row if \code{snp.in.col = FALSE})
    of the matrix used as input in \code{getLD} or \code{findLDblocks}.
}
  \item{end}{integer or character string specifying the index or name of the last SNP, respectively,
    that should be plotted.
}
  \item{xlab}{character string naming the label of the x-axis.
}
  \item{ylab}{character string naming the label of the y-axis.
} 
  \item{cexAxis}{a numeric value specifying the relative size of the SNP names displayed at the axes of the plot.
}  

  \item{block.col}{the color of the lines used to show the borders of the LD blocks.
}
  \item{block.lwd}{numeric value specifying the size of the lines used to show the borders of the LD blocks 

}
  \item{\dots}{further arguments of \code{image}.
}
}

\references{
  Gabriel, S.B. et al.~(2002). The Structure of Haplotype Blocks in the Human Genome. \emph{Science}, 296, 2225-2229.
}
\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
 \code{\link{findLDblocks}}, \code{\link{plot.getLD}}
}

\examples{
# Load the simulated data.
data(trio.data)

# Estimate LD blocks.
blocks <- findLDblocks(LDdata)

# Alternatively, the LD blocks can be estimated by
ld.out <- getLD(LDdata, addVarN=TRUE)
blocks2 <- findLDblocks(ld.out)

# Plot the LD blocks showing the Gabriel categorization.
plot(blocks)

# Plot the LD blocks showing the Dprime values.
plot(blocks, "Dprime")
}

\keyword{hplot}
