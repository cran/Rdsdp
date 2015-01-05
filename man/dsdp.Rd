\name{Rdsdp::dsdp}
\alias{dsdp}
\title{
Solve semidefinite programm with DSDP
}
\description{
Interface to DSDP semidefinite programming library.
}


\usage{
dsdp(A,b,C,K,OPTIONS=NULL)
}

\arguments{
  \item{A}{An object of class \code{"matrix"} with \eqn{m} rows defining the block diagonal constraint matrices
    \eqn{A_i}. Each constraint matrix \eqn{A_i} is specified by a row
	of \eqn{A} as explained in the Details section.}
  \item{b}{A numeric vector of length \eqn{m} containg the right hand side of the constraints.}
  \item{C}{An object of class \code{"matrix"} with one row or a valid
	class from the class hierarchy in the \code{"Matrix"} package. 
	It defines the objective coefficient matrix \eqn{C} with the same structure of \code{A} as explained above.}
  \item{K}{Describes the sizes of each block of the sdp problem. It is a list with the following elements:
    \describe{
      \item{\code{"s"}:}{A vector of integers listing the dimension of positive semidefinite cone blocks.}
      \item{\code{"l"}:}{A scaler integer indicating the dimension of the linear nonnegative cone block.}      
    }
  }	
  \item{OPTIONS}{A list of OPTIONS parameters passed to dsdp. It may contain any of the following fields:}
    \describe{
      \item{print:}{\code{= k} to display output at each k iteration, else = 0 [default 10].}
      \item{logsummary:}{\code{= 1}  print timing information if set to 1.}
      \item{save:}{to set the filename to save solution file in SDPA format.}
      \item{outputstats:}{\code{= 1} to output full information about the solution statistics in STATS.}
	  \item{gaptol:}{tolerance for duality gap as a fraction of the value of the objective functions [default 1e-6].}
      \item{maxit:}{maximum number of iterations allowed [default 1000].}
	  }
	Please refer to DSDP User Guide for additional OPTIONS parameters available. 
}
\details{
  All problem matrices are assumed to be of block diagonal structure, the input matrix \code{A} must be specified as follows:
  \enumerate{
    \item The coefficients for nonnegative cone block are put in the first \code{K$l} columns of \code{A}.
	\item The coefficients for positive semidefinite cone blocks are put after nonnegative cone block in the the same order as those in \code{K$s}.
	The \code{i}th positive semidefinite cone block takes \code{(K$s)[i]} times \code{(K$s)[[i]]} columns, with each row defining 
     a symmetric matrix of size \code{(K$s)[[i]]}.
  }
  

  
  This function does not check for symmetry in the problem data.
}

\value{
  Returns a list of three objects:
  \item{X}{Optimal primal solution \eqn{X}. A vector containing blocks in the
	same structure as explained above.}
  \item{y}{Optimal dual solution \eqn{y}. A vector of the same length as
    argument \code{b}}
  \item{STATS}{A list with three to eight fields that describe the solution of the problem:
    \describe{
      \item{stype:}{\code{PDFeasible} if the solutions to both (D) and (P) are feasible, 
	                \code{Infeasible} if (D) is infeasible, and \code{Unbounded} if (D) is unbounded.}
      \item{dobj:}{objective value of (D).}
      \item{pobj:}{objective value of (P).}
      \item{r:}{the multiple of the identity element added to \eqn{C-\mathcal{A}^{*}y} in the final solution to make \eqn{S} positive definite.}
      \item{mu:}{the final barrier parameter \eqn{\nu}.}
      \item{pstep:}{the final step length in (P)}
      \item{dstep:}{the final step length in (D)}
      \item{pnorm:}{the final value \eqn{\| P(\nu)\|}.}
    }
	The last five fields are optional, and only available when \code{OPTIONS$outputstats} is set to \eqn{1}.
    }
}

\references{
  \itemize{
	\item Steven J. Benson and Yinyu Ye:\cr
      \emph{{DSDP5} User Guide --- Software for Semidefinite Programming} Technical Report ANL/MCS-TM-277, 2005\cr
      \url{http://www.mcs.anl.gov/hs/software/DSDP/DSDP5-Matlab-UserGuide.pdf}
	}
}

\examples{
	K=NULL
	K$s=c(2,3)
	K$l=2

	C=matrix(c(0,0,2,1,1,2,c(3,0,1,
                       0,2,0,
                       1,0,3)),1,15,byrow=TRUE)
	A=matrix(c(0,1,0,0,0,0,c(3,0,1,
                         0,4,0,
                         1,0,5),
          	1,0,3,1,1,3,rep(0,9)), 2,15,byrow=TRUE)
	b <- c(1,2)
	
    OPTIONS=NULL
    OPTIONS$gaptol=0.000001
    OPTIONS$logsummary=0
    OPTIONS$outputstats=1
	
    result = dsdp(A,b,C,K,OPTIONS)
}
