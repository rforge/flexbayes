\name{posterior.subset}
\alias{posterior.subset}
\alias{[.posterior}
\alias{[[.posterior}
\title{Extract a subset of an object of class posterior}

\usage{
x[i,j]
x[[k]]
}

\arguments{
   \item{x}{An object of class \code{posterior}}
   \item{i}{Row to extract}
   \item{j}{Column to extract}
   \item{k}{Chain to extract}
}

\description{
   Methods for subsetting objects of class \code{posterior}.  
   
   \code{x[[k]]}: Returns a \code{posterior}
   object with a single chain, equal to the \code{k}th chain of \code{x}.  
   
   \code{x[,j]}, where \code{j} is an index, parameter name, or a vector of either of these:
   Returns a \code{posterior} object with only the parameters specified by \code{j}.

   \code{x[i,]}, where \code{i} is an index or a vector of indices:
   Returns a \code{posterior} object with only the iterations of \code{x} specified by \code{i}.

   \code{x[i,j]}: Returns a list of length equal to the number
   of chains of \code{x}; each element of the list is a scalar, vector, or matrix,
   depending on the form of \code{i} and \code{j}, containing the 
   \code{i}th iterations of the \code{j}th parameter(s) in the corresponding chain.
   
}
