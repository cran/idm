\name{enron}
\alias{enron}
\docType{data}
\title{
enron data set
}
\description{The data set is a subset of the Enron e-mail corpus from the UCI Machine Learning Repository (Lichman, 2013). The original data is a collection of 39,861 email messages with roughly 6 million tokens and a 28,102 term vocabulary. The subset is a binary (presence/absence) data set containing the 80 most frequent words which appear in the original corpus.
}
\usage{data("enron")}
\format{
  A binary data frame with 39,861 observations (e-mail messages) on 80 variables (words).
  }

\references{
Lichman, M. (2013). UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. Irvine, CA: University of California, School of Information and Computer Science.
}
\examples{
data(enron)
}
\keyword{datasets}
