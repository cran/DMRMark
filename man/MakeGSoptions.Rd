\name{MakeGSoptions}
\alias{MakeGSoptions}

\title{
	Encapsulate prior parameters and Gibbs Sampler (GS) control parameters
}
\description{
	This function encapsulate prior parameters and Gibbs Sampler 
	control parameters. All parameters with initial values.
	The encapsulation is for easy initiating, managing and 
	passing of parameters.
}
\usage{
MakeGSoptions(pi0 = c(100, 100, 5, 5), 
			  cmu0 = c(11.5, 11.5, 8, 8), 
			  theta0 = c(-3, 2), 
			  mu0 = matrix(c(-2, 2, 2, -2), 2, byrow = TRUE), 
			  kappa0 = c(50, 50, 5, 5), 
			  nu0 = rep(4, 2), 
			  A0 = array(rep(c(2, 0.8, 0.8, 4), 2), 
			  dim = c(2, 2, 2)), 
			  alpha12N = rep(40, 3), 
			  beta12N = rep(60, 3), 
			  D_mu = rep(-2, 2), 
			  chi_alpha = 0.2, #This and above for priors 
			  burnin = 500, #This and below for Gibbs Sampler Control
			  nsamples = 100, 
			  sampleSep = 10, 
			  onHMM = TRUE, 
			  track = FALSE, 
			  verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{pi0}{Length-4 vector, the concentration of Dirichlet distribution. 
  Prior of initial states.
}
  \item{cmu0}{Single value, the mean of Normal distribution.
  Prior of characteristic length.
}
  \item{theta0}{Length-2 vector, each value is the mean of a Normal distributions.
  Priors for means of control groups of two non-differentially 
  methylated CpG sites (non-DMCs) responses.
}
  \item{mu0}{2-by-2 matrix, each row is the means of a bivariate Normal distributions.
  Priors for means of two DMCs responses
}
  \item{kappa0}{Length-4 vector, each value is the prior observation number of 
  Normal-Inverse-Gamma (NIG) or Normal-Inverse-Wishart (NIW) depends on the corresponding
  state. 
}
  \item{nu0}{Length-2 vector, each value is the degree 
  of freedom of an IW distribution.
  Priors for covariance of DMC responses.
}
  \item{A0}{2-by-2-by-2 array, each 2-by-2 matrix along the 
  third dimension is the scale matrix of an IW distribution.
  Priors for covariance of DMC responses.
}
  \item{alpha12N}{Length-3 vector, each value is the shape of an IG distribution.
  Priors for variance of non-DMC responses.
}
  \item{beta12N}{Length-3 vector, each value is the rate of an IG distribution.
  Priors for variance of non-DMC responses.
}
  \item{D_mu}{Length-2 vector, each value is the minimum 
  distance between two group means of DMCs.
  Prior for truncating the means of bivariate normals of DMC's responses.
}
  \item{chi_alpha}{p-value of chi-square distribution with 2 degrees of freedom.
  Prior for truncating the covariant matrices 
  of bivariate normals of DMC's responses.
}

  \item{burnin}{Number of iterations for burn-in.
  Gibbs Sampler control parameter. Default is 500.
}
  \item{nsamples}{Number of samples to compute the point estimators.
  Gibbs Sampler control parameter. Default is 100.
}
  \item{sampleSep}{Only keep every 'sampleSep'-th samples to estimate point 
  estimators. Gibbs Sampler control parameter. Default is 10.
}
  \item{onHMM}{Set to FALSE will disable HMM, and reduce to simple clustering
  of Mixture Model. Gibbs Sampler control parameter. Default is TRUE.
}
  \item{track}{Set to TRUE will make DMRMark return all samples 
  from the beginning of burn-in to the end of sampling
  instead of point estimators. Useful for inspecting convergence. 
  Please know well about this issue before you decide to set it to TRUE.
  Gibbs Sampler control parameter. Default is TRUE.
}
  \item{verbose}{Set to TRUE to show the details when running the 
  Gibbs Sampler. Gibbs Sampler control parameter. Default is FALSE.
}
}

\value{
	Simply a list with all items are the same with input.
	Just an encapsulation.
}

\author{
	Linghao SHEN <sl013@ie.cuhk.edu.hk>
}

\seealso{
	\code{\link{DMRMark}}
}
\examples{

	# MakeGSoptions
	opts <- MakeGSoptions()
}

