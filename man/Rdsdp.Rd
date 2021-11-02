\name{Rdsdp}
\alias{Rdsdp}
\docType{package}
\title{
R interface to DSDP semidefinite programming library
}
\description{
Rdsdp is the R package providing a R interface to DSDP semidefinite programming library. The DSDP package implements a dual-scaling algorithm to
find solutions (\eqn{X}, \eqn{y}) to linear and semidefinite optimization problems of the form
	\deqn{\mbox{(P)} \ \inf\, \mbox{tr}(CX)}{inf tr(CX)}
	\deqn{\mbox{subject to}\; \mathcal{A}X = b}{subject to AX = b}
	\deqn{X \succeq 0}{X >= 0}
	with
	\eqn{(\mathcal{A}X)_i = \mbox{tr}(A_iX)}{(AX)_i = tr(A_iX)}
	where \eqn{X \succeq 0}{X >= 0} means X is positive
	semidefinite, \eqn{C} and all \eqn{A_i} are symmetric matrices of the same
	size and \eqn{b} is a
	vector of length \eqn{m}. 
	
	The dual of the problem is
	\deqn{\mbox{(D)} \ \sup\, b^{T}y}{sup b'y}
	\deqn{\mbox{subject to}\; \mathcal{A}^{*}y + S = C}{subject to A'y + S = C}
	\deqn{S \succeq 0}{S >= 0}

	where
	\eqn{\mathcal{A}y = \sum_{i=1}^m y_i A_i.}{A'y = \sum_{i=1}^m y_i A_i.}

	Matrices \eqn{C} and \eqn{A_i} are assumed to be block diagonal
	structured, and must be specified that way (see Details). 
}

\references{
  \itemize{
	\item \url{https://www.mcs.anl.gov/hs/software/DSDP/}
	\item Steven J. Benson and Yinyu Ye:\cr
      \emph{Algorithm 875: DSDP5 software for semidefinite programming} ACM Transactions on Mathematical Software (TOMS)  34(3), 2008\cr
      \url{http://web.stanford.edu/~yyye/DSDP5-Paper.pdf}
	\item Steven J. Benson and Yinyu Ye and Xiong Zhang:\cr
      \emph{Solving Large-Scale Sparse Semidefinite Programs for Combinatorial Optimization} SIAM Journal on Optimization  10(2):443-461, 2000\cr
      \url{http://web.stanford.edu/~yyye/yyye/largesdp.ps.gz}
  }
}

