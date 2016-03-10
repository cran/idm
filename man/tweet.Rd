\name{tweet}
\alias{tweet}
\docType{data}
\title{
twitter data set
}
\description{
The data set refers to a small corpus of messages or tweets mentioning seven
major hotel brands. It was gathered by continuously querying and archiving
the Twitter Streaming API service, using the \code{twitteR} package in \code{R}. A total of 7,296 tweets were extracted within a time period of 6 days, from June 23th to June 28th 2013. Only tweets in the English language were considered. A sentiment polarity variable was calculated, indicating the sentiment value of each message and a third variable, user visibility or popularity, as measured by
the number of followers each user had, was also included in the dataset
}
\usage{data("tweet")}
\format{
  A data frame with the following variables:
  \describe{
    \item{\code{Brand}}{The hotel brand mentioned in the tweet: 1=\kbd{Hilton}, 2=\kbd{Intercontinental}, 3=\kbd{Marriott}, 4=\kbd{Bestwestern}, 5=\kbd{Starwood}, 6=\kbd{Hyatt}, 7=\kbd{Choice}}
    \item{\code{Sentiment}}{Sentiment  for each tweet: 1=\kbd{negative (-)}, 2=\kbd{mixed (+/-)}, 3=\kbd{positive (+)}, 4=\kbd{very positive (++)}}
    \item{\code{UserVis}}{User popularity/visibility in Twitter: 1=\kbd{low}, 2=\kbd{medium}, 3=\kbd{high}}
  }
}

\references{
Iodice D' Enza, A., & Markos, A. (2015). Low-dimensional tracking of association structures in categorical data, \emph{Statistics and Computing}, \emph{25}(5), 1009-1022.
}
\examples{
data(tweet)
}
\keyword{datasets}
