\name{DMRMark-package}

\alias{DMRMark-package}

\docType{package}
\title{
DMR Detection by Non-Homogeneous Hidden Markov Model from Methylation Array Data
}

\description{

Perform differential analysis for
methylation array data. DMRMark detects differentially
methylated regions (DMRs) from array M-values. 
The core of DSS is a Non-homogeneous Hidden Markov Model
for estimating spatial correlation and a novel Constrained 
Gaussian Mixture Model for modeling the M-value pairs of each individual locus.

DMRMark only works for two-group comparisons currently. We have the plan to
extend the transition and response model that make then suitable for complex
experimental designs in the future. 
}

\author{
	Linghao SHEN <sl013@ie.cuhk.edu.hk>
}

\keyword{ package }

\examples{
	# DMR detection performed on chr18 of a small BLCA dataset from TCGA
	data(BLCA)
	
	# Use a small subset
	nprobe <- 500 
	# M-values
	mv <- BLCA$mv[1:nprobe,]
	
	# Distance between probes, L<0 indicates acorssing chromosomes 
	L = BLCA$distance[1:nprobe]
	
	# Initialize new chain when probe distance too long 
	# or across different chromosomes
	newChains <- which((L > 100000) | L < 0)
	# The starting positions of new chains
	starting <- c(1, newChains[-length(newChains)]+1)
	
	# Run DMRMark with default options
	set.seed(0)
	par <- DMRMark(mv, L, starting)
	
	# Get the posterior of being certain states
	# Return the result of DML for plotting by setting 'region=FALSE'
	results <- DMRViterbi(mv, par, L, starting, region=FALSE)
	
	# The MAP states being 3 or 4 indicate DMLs
	isDML <- (results$states > 2) + 0
	mvScatter(mv, isDML, nPlot=10000)
}
