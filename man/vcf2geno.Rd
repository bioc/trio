\name{vcf2geno}
\alias{vcf2geno}

\title{
Transformation of VFC File
}
\description{
Transforms a vcf file into a matrix in genotype format required by, e.g., the functions for computing the genotypic TDT.
}
\usage{
vcf2geno(vcf, ped, none = "0/0", one = c("0/1"), both = "1/1", use.rownames = FALSE, 
    allowDifference = FALSE, removeMonomorphic = TRUE, changeMinor = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf}{a matrix resulting from reading a vcf file into \code{R}, or an object of class \code{collapsedVCF} (i.e. the output
	  of, e.g., the function \code{readVcf} from the \code{VariantAnnotation} package). If \code{use.rownames = FALSE}, 
		the column names of the genotype matrix must correspond to the personal IDs in \code{ped} (i.e. either the column \code{pid}
		of \code{ped}, if the entries in \code{pid} are unique, or otherwise, a combination of the columns \code{famid} and 
		\code{pid} from \code{ped}, combined using an underscore). If \code{use.rownames = TRUE}, the column names of the genotype
		matrix specified by \code{vcf} must correspond to the row names of \code{ped}.
}
  \item{ped}{a data frame containing the family information for the subjects in \code{vcf} (might also contain information for
	  other subjects, see \code{allowDifference}). This data frame must contain the columns \code{famid}, \code{pid}, \code{fatid},
		and \code{motid} comprising the family ID, the personal ID as well as the ID of the father and the mother, respectively.
}
  \item{none}{a character string or vector specifying the coding for the homozygous reference genotype.
}
  \item{one}{a character string or vector specifying the coding for the heterozygous genotype.
}
  \item{both}{a character string or vector specifying the coding for the homozygous variant genotype.
}
  \item{use.rownames}{a logical value specifying whether the row names of \code{ped} correspond to the sample names in 
	  \code{vcf}. For details, see \code{vcf}.
}
  \item{allowDifference}{a logical value specifying whether \code{ped} and \code{vcf} are allowed to also contain samples not
	  available in the respective other object. If \code{FALSE}, all samples in \code{ped} must also be available in \code{vcf},
		and vice versa (matched as described in \code{vcf}). If \code{TRUE}, at least 10\% of the samples must be contained in both
		\code{vcf} and \code{ped}.
}
  \item{removeMonomorphic}{a logical value specifying whether monomorphic SNVs should be removed from the output.
}
  \item{changeMinor}{a logical value specifying whether the coding of the genotypes should be changed for SNVs for which the 
	  default coding leads to a minor allele frequency larger than 0.5. The genotypes are coded by the number of minor alleles,
		i.e. the genotype(s) specified by \code{none} is coded by 0, the genotype(s) specified by \code{one} is coded by 1,
		and the genotype(s) specified by \code{both} is coded by 2. If for an SNV this leads to a minor allele frequency larger 
		than 0.5 and \code{changeMinor = TRUE}, this 0, 1, 2-coding will be changed into a 2, 1, 0-coding.
}
}

\value{
A matrix in genotype format required, e.g., by functions for performing different types of the genotypic TDT, such as 
\code{\link{colTDT}}.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}


\seealso{
\code{\link{colTDT}}, \code{\link{colGxG}}, \code{\link{colGxE}}, \code{\link{ped2geno}}
}

\keyword{manip}

