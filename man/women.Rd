\name{women}
\alias{women}
\docType{data}
\title{
women data set
}
\description{
The data are from the third Family and Changing Gender Roles survey conducted in 2002. The questions retained are those related to working women in Spain and the effect on the family. A total of 2,107 respondents answered eight questions on a 5-point Likert scale, as well as four demographic variables (gender, martial status, education and age). There are no cases with missing data.
}
\usage{data("women")}
\format{
  A data frame with the following variables:
  \describe{
    \item{\code{A}}{"a working mother can establish a warm relationship with her child" \cr 1=\kbd{strongly agree}, 2=\kbd{agree}, 3=\kbd{neither agree or disagree}, 4=\kbd{disagree}, 5=\kbd{strongly disagree}}
    \item{\code{B}}{"a pre-school child suffers if his or her mother works" \cr 1=\kbd{strongly agree}, 2=\kbd{agree}, 3=\kbd{neither agree or disagree}, 4=\kbd{disagree}, 5=\kbd{strongly disagree}}
    \item{\code{C}}{"when a woman works the family life suffers" \cr 1=\kbd{strongly agree}, 2=\kbd{agree}, 3=\kbd{neither agree or disagree}, 4=\kbd{disagree}, 5=\kbd{strongly disagree}}
    \item{\code{D}}{"what women really want is a home and kids" \cr 1=\kbd{strongly agree}, 2=\kbd{agree}, 3=\kbd{neither agree or disagree}, 4=\kbd{disagree}, 5=\kbd{strongly agree}}
    \item{\code{E}}{"running a household is just as satisfying as a paid job" \cr 1=\kbd{strongly agree}, 2=\kbd{agree}, 3=\kbd{neither agree or disagree}, 4=\kbd{disagree}, 5=\kbd{strongly disagree}}
  \item{\code{F}}{"work is best for a woman's independence" \cr 1=\kbd{strongly agree}, 2=\kbd{agree}, 3=\kbd{neither agree or disagree}, 4=\kbd{disagree}, 5=\kbd{strongly disagree}}
    \item{\code{G}}{"a man's job is to work; a woman's job is the household" \cr 1=\kbd{strongly agree}, 2=\kbd{agree}, 3=\kbd{neither agree or disagree}, 4=\kbd{disagree}, 5=\kbd{strongly disagree}}
    \item{\code{H}}{"working women should get paid maternity leave" \cr 1=\kbd{strongly agree}, 2=\kbd{agree}, 3=\kbd{neither agree or disagree}, 4=\kbd{disagree}, 5=\kbd{strongly disagree}}
    \item{\code{g}}{gender: 1=\kbd{male}, 2=\kbd{female}}
    \item{\code{m}}{marital status: 1=\kbd{married/living as married}, 2=\kbd{widowed}, 3=\kbd{divorced}, 4=\kbd{separated, but married}, 5=\kbd{single, never married}}
    \item{\code{e}}{education: 1=\kbd{no formal education}, 2=\kbd{lowest education}, 3=\kbd{above lowest education}, 4=\kbd{highest secondary completed}, 5=\kbd{above higher secondary level, below full university}, 6=\kbd{university degree completed}}
    \item{\code{a}}{age: 1=\kbd{16-25 years}, 2=\kbd{26-35}, 3=\kbd{36-45}, 4=\kbd{46-55}, 5=\kbd{56-65}, 6=\kbd{66 and older}}
    }
}
\source{
\url{http://www.econ.upf.edu/~michael/women_Spain2002_original.xls}
}
\references{
Greenacre, M. J. (2010). \emph{Biplots in practice}. Fundacion BBVA.
}
\examples{
data(women)
}
\keyword{datasets}
