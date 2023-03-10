% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networkmap_ggnetwork.R
\name{networkmap_ggnetwork}
\alias{networkmap_ggnetwork}
\title{Network map: make ggraph}
\usage{
networkmap_ggnetwork(
  tg,
  layout = "nicely",
  node_size_range = NULL,
  edge_alpha = 0.25,
  seed = 2023
)
}
\arguments{
\item{layout}{The type of layout to create. Either a valid string, a
function, a matrix, or a data.frame (see Details)}

\item{node_size_range}{The minimum / maximum size of each node.}
}
\value{
\link[ggnetwork]{ggnetwork}
}
\description{
Subfunction of \link[ThreeWayTest]{networkmap}.
}
\keyword{internal}
