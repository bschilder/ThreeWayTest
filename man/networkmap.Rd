% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networkmap.R
\name{networkmap}
\alias{networkmap}
\title{Network map}
\source{
\href{https://ggraph.data-imaginist.com/articles/Layouts.html}{
 ggraph layouts}

\href{https://github.com/thomasp85/ggraph/issues/75}{
ggraph issue with scale functions}
}
\usage{
networkmap(
  dat,
  cols = grep("^w", names(dat), value = TRUE),
  node_vars = c("dataset", "GENE", "SNP"),
  i = seq_len(50),
  snps = unique(dat$SNP)[i],
  agg_var = NULL,
  agg_fun = mean,
  as_cor = FALSE,
  k_row = 3,
  annot_vars = c("TYPE", "GENE"),
  show_plot = TRUE,
  layout = "nicely",
  node_size_range = NULL
)
}
\arguments{
\item{dat}{Summary statistics data.}

\item{cols}{Numeric columns to plot in the heatmap.}

\item{node_vars}{Columns within \code{dat} to set as nodes. 
Each node variable will link to the next node in the character vector.
You can repeat column names to create more connections between nodes.}

\item{i}{Indices of rows to include. Set to \code{NULL} to include all rows, 
but be warned that this can become very computationally expensive.}

\item{snps}{A character vector of SNP RSIDs to subset \code{dat} by.}

\item{agg_var}{Variable to aggregate data by. 
Set to \code{NULL} to skip this step.}

\item{agg_fun}{Function to aggregate \code{cols} with.}

\item{as_cor}{Show the heatmap as a correlation matrix
instead of a feature x sample matrix.}

\item{k_row}{an integer scalar with the desired number of groups by which to
color the dendrogram's branches in the rows (uses \link[dendextend]{color_branches})
If NA then \link[dendextend]{find_k} is used to deduce the optimal number of clusters.}

\item{annot_vars}{Variables in \code{dat} to include as row-wise annotations.}

\item{show_plot}{Print the plot.}

\item{layout}{The type of layout to create. Either a valid string, a
function, a matrix, or a data.frame (see Details)}

\item{node_size_range}{The minimum / maximum size of each node.}
}
\value{
A named list containing 
a network plot and the data used to create it.
}
\description{
Generate an interactive network map of the summary statistics data.
}
\examples{
dat <- ThreeWayTest::data_matrix_final
nm <- networkmap(dat = dat)
}
