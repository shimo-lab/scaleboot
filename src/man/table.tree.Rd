\name{table.tree}
\alias{table.tree}
\alias{prepare.table.tree}
\title{Tables for phylogenetic inference}
\description{
Creating tables of p-values and tree/edge associaitons for phylogenetic inference.
Trees are sorted by the likelihood value.
}
\usage{
table.tree(x, k = 2)

prepare.table.tree(trees,edges,edge2tree,
                   treename=NULL,edgename=NULL,taxaname=NULL,mt=NULL) 
}
\arguments{
  \item{x}{output of \code{prepare.table.tree}.}
  \item{k}{integer of \eqn{k} for calculating p-values.}
  \item{trees}{\code{relltest} output for trees.}
  \item{edges}{\code{relltest} output for edges.}
  \item{edge2tree}{\code{read.ass} output for "edge to tree" association.}
  \item{treename}{character vector for tree descriptions.}
  \item{edgename}{character vector for edge descriptions.}
  \item{taxaname}{character vector for taxa names.}
  \item{mt}{\code{read.mt} output for the site-wise log-likelhiood values.}
}
\details{
First, apply \code{prepare.table.tree} to prepare a list of information used by \code{table.tree}.
Output tables are suitable for publication.
}
\value{
\code{prepare.table.tree} returns a list of several information, including p-values.
Then \code{table.tree} returns a list of two character matrices; 
one for tree table and the other for edge table.
}

\author{Hidetoshi Shimodaira}

\seealso{\code{\link{relltest}}, \code{\link{read.ass}}, \code{\link{read.mt}}}

\examples{
data(mam15)
mam15.trees <- mam15.relltest[1:15]  # 15 trees
mam15.edges <- mam15.relltest[16:25]  # 10 edges
mam15.edge2tree <- mam15.ass[16:25] # 10 edges
mam15 <- prepare.table.tree(mam15.trees,mam15.edges,mam15.edge2tree)
table.tree(mam15)
}

\keyword{ phylogenetics }

