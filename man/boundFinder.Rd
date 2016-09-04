\name{boundFinder}
\alias{boundFinder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Find a pair of reasonable distances of group means 
for hyper- and hypomethylation based on the quantile of two-group difference.
}
\description{
	This function takes the M-values to produce the distance $D$ to be the maximum 
value satisfies that the proportion of absolute values of two-group difference larger
than $D$ is at certain level. Due to the precision limitation, the $D$'s for 
hyper- and hypomethylation are not necessarily the same. 
If the samples are not totally paired, than user should
first call 'reformData' to process M-values.
}
\usage{
boundFinder(mv, prop = 0.1)
}

\arguments{
  \item{mv}{The input M-values matrix. If the samples are not totally 
  paired, than user MUST first call "reformData" to process M-values.
}
  \item{prop}{The proportion that absolute values of two-group difference larger than $D$
	must be satisfied. Default value is 0.1
}
}
\details{
	The choices of 'prop' are not going to be too extreme or stringent, which will produce
	to dominated prior. This value should depend on the belief that around 'prop' proportion
	of loci are differentially methylated. In general 0.1 to 0.2 should be reasonable and
	and well-performed. Users may also choose different $D$'s for two differentially
	methylation status freely, in this situation, values around 1.5 to 3 are recommended.
	
	Users must ensure the M-values come from paired samples or has been processed by 
	'reformData' according to the experiment design.
}
\value{ A two-value vector contains $D_1$ and $D_2$ for the group-mean difference of
hypermethylation and hypomethylation respectively. Due to the precision limitation, 
$D_1$ may not necessarily equal to $D_2$.
}
\author{
	Linghao SHEN <sl013@ie.cuhk.edu.hk>
}

\seealso{
	\code{\link{reformData}} to tackle unpaired data.
}
\examples{

	# Finding the 5% and 95% quantile of normal samples
	set.seed(0)
	mv <- cbind(rep(0,100000),rnorm(100000))
	boundFinder(mv)
	
	# Output matched the normal p-values
	#     5.0%     94.9% 
	#-1.639578  1.639691
}
