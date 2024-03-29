% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccmean.R
\name{ccmean}
\alias{ccmean}
\title{Calculates estimates of the mean cost with censored data}
\usage{
ccmean(x, L = max(x$surv), addInterPol = 0)
}
\arguments{
\item{x}{A dataframe with columns: id, cost, delta and surv. If Cost history is available it can be specified by: start and stop,}

\item{L}{Limit. Mean cost is calculated up till L, if not specified L = max(surv)}

\item{addInterPol}{This parameter affects the interpolation of cost between two observed times. Defaults to zero.}
}
\value{
An object of class "ccobject".
}
\description{
This function calculates the mean cost for right-censored cost 
data over a period of L time units (days, months, years,...)
}
\details{
The function returns four estimates. The first two are simple and biased
downwards, and included for comparison. The estimates are:

- AS: "Available Sample estimator" - The simple sample mean

- CC: "Complete Case estimator" - The mean of fully observed cases

- BT: "Weighted Complete Case estimator" - Bang and Tsiatis's estimator

- ZT: "Weighted Available estimator" - Zhao and Tian's estimator

The function needs the following in a dataframe:

- id: The id separating each individual

- cost: The total cost, or if start and stop provided the specific cost

- start: Start of cost

- stop: End of cost, if one time cost then start = stop

- delta: Event variable, 1 = event, 0 = no event

- surv: Survival
}
\examples{
hcost
ccmean(hcost, L = 1461, addInterPol = 1)

}
\references{
\insertRef{Bang2000}{ccostr}

\insertRef{Zhao2001}{ccostr}
}
