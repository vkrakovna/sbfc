##' @useDynLib sbfc
##' @importFrom Rcpp evalCpp
##' @importFrom DiagrammeR grViz
##' @importFrom Matrix sparseMatrix
##' @importFrom discretization mdlp
##' @importFrom graphics axis barplot plot
##' @importFrom stats acf cor

##' @title Data set discretization and formatting
##' @description Removes rows containing missing data, and discretizes the data set using Minimum Description Length Partitioning (MDLP).
##' @param data Data frame, where the last column must be the class variable.
##' @param n_train Number of data frame rows to use as the training set - the rest are used for the test set. If NULL, all rows are used for training, and there is no test set (default=NULL).
##' @param missing Label that denotes missing values in your data frame (default='?').
##' @return A discretized data set:
##' \describe{     
##' \item{\code{TrainX}}{Matrix containing the training data.}
##' \item{\code{TrainY}}{Vector containing the class labels for the training data.}
##' \item{\code{TestX}}{Matrix containing the test data (optional).}
##' \item{\code{TestY}}{Vector containing the class labels for the test data (optional).}
##' }
##' @examples data(iris)
##' iris_disc = data_disc(iris)
##' @export
data_disc = function(data, n_train = NULL, missing = '?') {
  for (i in 1:ncol(data)) { 
    mrows = grep(missing, data[,i], fixed=T)
    if (length(mrows) > 0) data = data[-mrows,] # remove missing rows
    if (!is.null(n_train)) n_train = n_train - sum(mrows <= n_train)
  }
  for (i in 1:ncol(data)) data[,i] = as.numeric(data[,i])
  data = mdlp(data)$Disc.data
  if (is.null(n_train)) n_train = nrow(data)

  data_disc = list()
  data_disc$TrainY = data[1:n_train,ncol(data)]
  data_disc$TrainX = data[1:n_train,-ncol(data)]
  if (n_train < nrow(data)) {
    data_disc$TestY = data[(n_train+1):nrow(data),ncol(data)]
    data_disc$TestX = data[(n_train+1):nrow(data),-ncol(data)]
  }
  data_disc
}

##' @title
##' Selective Bayesian Forest Classifier (SBFC) algorithm
##' @description
##' Runs the SBFC algorithm on a discretized data set. To discretize your data, use the \code{\link{data_disc}} command.
##' 
##' @param data Discretized data set:
##' \describe{     
##' \item{\code{TrainX}}{Matrix containing the training data.}
##' \item{\code{TrainY}}{Vector containing the class labels for the training data.}
##' \item{\code{TestX}}{Matrix containing the test data (optional).}
##' \item{\code{TestY}}{Vector containing the class labels for the test data (optional).}
##' }
##' @param nstep Number of MCMC steps, default max(10000, 10 * ncol(TrainX)).
##' @param thin Thinning factor for the MCMC. 
##' @param burnin_denom Denominator of the fraction of total MCMC steps discarded as burnin (default=5).
##' @param cv Do cross-validation on the training set (if test set is not provided).
##' @param thinoutputs Return thinned MCMC outputs (parents, groups, trees, logposterior), rather than all outputs (default=FALSE).
##' @param alpha Dirichlet hyperparameter(default=1)
##' @param y_penalty Prior coefficient for y-edges, which penalizes signal group size (default=1)
##' @param x_penalty Prior coefficient for x-edges, which penalizes tree size (default=4)
##' @details
##' Data needs to be discretized before running SBFC. \cr
##' If the test data matrix TestX is provided, SBFC runs on the entire training set TrainX, and provides predicted class labels for the test data. 
##' If the test data class vector TestY is provided, the accuracy is computed. 
##' If the test data matrix TestX is not provided, and cv is set to TRUE, SBFC performs cross-validation on the training data set TrainX, 
##' and returns predicted classes and accuracy for the training data. \cr
##' @return An object of class \code{sbfc}:
##' \describe{     
  ##' \item{\code{accuracy}}{Classification accuracy (on the test set if provided, otherwise cross-validation accuracy on training set).}
  ##' \item{\code{predictions}}{Vector of class label predictions (for the test set if provided, otherwise for the training set).}
  ##' \item{\code{probabilities}}{Matrix of class label probabilities (for the test set if provided, otherwise for the training set).}
  ##' \item{\code{runtime}}{Total runtime of the algorithm in seconds.}
  ##' \item{\code{parents}}{Matrix representing the structures sampled by MCMC, where parents[i,j] is the index of the parent of node j at iteration i (0 if node is a root).}
  ##' \item{\code{groups}}{Matrix representing the structures sampled by MCMC, where groups[i,j] indicates which group node j belongs to at iteration j (0 is noise, 1 is signal).}
  ##' \item{\code{trees}}{Matrix representing the structures sampled by MCMC, where trees[i,j] indicates which tree node j belongs to at iteration j.}
  ##' \item{\code{logposterior}}{Vector representing the log posterior at each iteration of the MCMC.}
  ##' \item{Parameters}{\code{nstep}, \code{thin}, \code{burnin_denom}, \code{cv}, \code{thinoutputs}, \code{alpha}, \code{y_penalty}, \code{x_penalty}.}
  ##' }
  ##' If \code{cv=TRUE}, the MCMC samples from the first fold are returned (\code{parents}, \code{groups}, \code{trees}, \code{logposterior}).
##' @examples
##' data(madelon)
##' madelon_result = sbfc(madelon)
##' data(heart)
##' heart_result = sbfc(heart, cv=FALSE)
##' @export
sbfc = function(data, nstep = NULL, thin = 50, burnin_denom = 5, cv = T, thinoutputs = F, 
                alpha = 5, y_penalty = 1, x_penalty = 4) {
  sbfc_cpp(if (is.null(data$TrainX)) NULL else apply(as.matrix(data$TrainX), 2, as.integer), 
       if (is.null(data$TrainY)) NULL else as.integer(data$TrainY), 
       if (is.null(data$TestX)) NULL else apply(as.matrix(data$TestX), 2, as.integer), 
       if (is.null(data$TestY)) NULL else as.integer(data$TestY),
       nstep, thin, burnin_denom, cv, thinoutputs, alpha, y_penalty, x_penalty)
}

#### SBFC graphs ####

# produces Graphviz code for a graph for a single MCMC sample
single_sbfc_graph = function(groups, parents, i, single_noise_nodes=F, 
                             names=paste0("X", 1:ncol(parents)), colorscheme="blues", ncolors=7) {
  s = paste0('digraph G { subgraph cluster_g1 {
  node [colorscheme=', colorscheme, ncolors, ' color=6, fontcolor=white, style=filled, fontname=Arial]; label="Group 1";')
  for (j in ncol(parents):1) {
    if(groups[i, j] == 1)
      s = paste(s, " \"", names[j], "\";")
    if(groups[i, j] == 1 && parents[i, j] != 0)
      s = paste(s,  " \"", names[parents[i, j]], "\" -> \"", names[j], "\";")
  }
  s = paste0(s, '} subgraph cluster_g0 {
            node [colorscheme=', colorscheme, ncolors, ' color=2, style=filled, fontname=Arial]; label="Group 0";')
  for (j in ncol(parents):1) {
    if((groups[i, j] == 0) && single_noise_nodes)
      s = paste(s, " \"", names[j], "\";")
    if(groups[i, j] == 0 && parents[i, j] != 0)
      s = paste(s,  " \"", names[parents[i, j]], "\" -> \"", names[j], "\";")
  }
  s = paste(s, "}}")
  s
}

# determines a set of edges to include in the average graph
average_sbfc_graph_edges = function(parents, cutoff=0.1, names=paste0("X", 1:ncol(parents))) {
  j1=c(); j2=c(); freq=c();
  for (j in 1:ncol(parents)) {
    par = sort(unique(parents[,j]))
    freq_edge = table(parents[,j])/nrow(parents)
    stopifnot(length(par) == length(freq_edge))
    if (par[1] == 0) {
      if (length(par) == 1) next
      par = par[-1]
      freq_edge = freq_edge[-1]
    }
    j1 = c(j1, rep(j, length(par)), par)
    j2 = c(j2, par, rep(j, length(par)))
    freq = c(freq, freq_edge, freq_edge)
  }
  stopifnot(length(j1) == length(j2), length(j1) == length(freq))
  if (length(j1) == 0) return(list(s="", edge_nodes=c()))
  freq_mat = sparseMatrix(i=j1, j=j2, x=freq)
  mat = cbind(freq_mat@i+1, rep(seq_along(diff(freq_mat@p)), diff(freq_mat@p)), freq_mat@x)
  s =""
  edge_nodes = c()
  for (i in 1:nrow(mat)) {
    if ((mat[i, 1] < mat[i, 2]) && (mat[i, 3] >= cutoff)) {
      edge_nodes = c(edge_nodes, mat[i, 1], mat[i, 2])
      s = paste0(s, " \"", names[mat[i, 1]], "\" -- \"", names[mat[i, 2]], "\" [penwidth=", 5*mat[i, 3], "];")
    }
  }
  list(s = s, edge_nodes = unique(edge_nodes))
}

# produces Graphviz code for an average graph over a set of MCMC sample graphs
average_sbfc_graph = function(groups, parents, edge_cutoff=0.1, single_noise_nodes=F,
                              names = paste0("X", 1:ncol(parents)), colorscheme="blues", ncolors=7) {
  freq.group1 = apply((groups == 1), 2, mean)
  freq.group0 = apply((groups == 0), 2, mean)
  ae = average_sbfc_graph_edges(parents, edge_cutoff, names = names)
  vars = 1:ncol(parents)
  #if (edges_only) vars = ae$edge_nodes
  if (!single_noise_nodes) vars = unique(c(ae$edge_nodes, which(freq.group1 >= .2))) #add parameter
  fontcolor=c(rep("black", floor(ncolors/2)), rep("white", ceiling(ncolors/2)))
  
  s = "strict graph G { node[fontname=Arial style=filled] "
  for (i in ncolors:1) {
    for (j in vars) {
      if ((freq.group1[j] >= (i-1)*1.0/ncolors) && (freq.group1[j] <= i*1.0/ncolors))
        s = paste0(s, " \"", names[j], "\" [colorscheme=", colorscheme, ncolors,
                   " fontcolor=", fontcolor[i], " color=", i, "];")
    }
  }
  s = paste(s, ae$s, "}")
  s
}

##' @title SBFC graph
##' @description Plots a sampled MCMC graph or an average of sampled graphs using Graphviz. \cr
##' In average graphs, nodes are color-coded according to importance - the proportion of samples where the node appeared in Group 1 (dark-shaded nodes appear more often).
##' In average graphs, thickness of edges also corresponds to importance: the proportion of samples where the edge appeared.
##' @param sbfc_result An object of class \code{sbfc}.
##' @param average Plot an average of sampled MCMC graphs (default=TRUE).
##' @param iter MCMC iteration of the sampled graph to plot, if \code{average=F} (default=10000).
##' @param edge_cutoff The average graph includes edges that appear in at least this fraction of the sampled graphs, if \code{average=T} (default=0.1).
##' @param single_noise_nodes Plot single-node trees that appear in the noise group (Group 0) in at least 80 percent of the samples, which can be numerous for high-dimensional data sets (default=FALSE).
##' @param labels A vector of node labels (default=\code{c("X1","X2",...)}).
##' @param save_graphviz_code Save the Graphviz source code in a .gv file (default=FALSE).
##' @param colorscheme \href{http://www.graphviz.org/doc/info/colors.html}{Graphviz color scheme} for the nodes (default="blues").
##' @param ncolors number of colors in the palette (default=7).
##' @param width An optional parameter for specifying the width of the resulting graphic in pixels.
##' @param height An optional parameter for specifying the height of the resulting graphic in pixels.
##' @examples
##' data(madelon)
##' madelon_result = sbfc(madelon)
##' sbfc_graph(madelon_result) 
##' sbfc_graph(madelon_result, average=FALSE, iter=5000) # graph for 5000th iteration
##' sbfc_graph(madelon_result, single_noise_nodes=TRUE) # wide graph with 480 single nodes
##' 
##' data(heart)
##' heart_result = sbfc(heart)
##' heart_labels = c("Age", "Sex", "Chest Pain", "Rest Blood Pressure", "Cholesterol", 
##' "Blood Sugar", "Rest ECG", "Max Heart Rate", "Angina", "ST Depression", "ST Slope",
##' "Fluoroscopy Colored Vessels", "Thalassemia")
##' sbfc_graph(heart_result, labels=heart_labels, width=700)
##' @export
sbfc_graph = function(sbfc_result, iter=10000, average=T, edge_cutoff=0.1, single_noise_nodes=F,
                      labels=paste0("X", 1:ncol(sbfc_result$parents)), save_graphviz_code = F,
                      colorscheme="blues", ncolors=7, width=NULL, height=NULL) {
  parents = sbfc_result$parents
  groups = sbfc_result$groups
  if (length(labels) != ncol(parents)) stop("Size of label vector must be equal to number of variables.")
  labels = gsub(" ", "\n", labels)
  if (average) {
    rows = seq(sbfc_result$nstep/sbfc_result$burnin_denom + 1, sbfc_result$nstep, by=sbfc_result$thin)
    gv_source = average_sbfc_graph(groups[rows,], parents[rows,], edge_cutoff, single_noise_nodes, labels, colorscheme, ncolors)
  } else 
    gv_source = single_sbfc_graph(groups, parents, iter, single_noise_nodes, labels, colorscheme, ncolors)
  if (save_graphviz_code) writeLines(gv_source, "sbfc_graph_code.gv")
  if (is.null(width)) 
    width = 500 + (40 + 10*average) * max(0, length(gregexpr(';', gv_source)[[1]])-2*length(gregexpr(' -', gv_source)[[1]])-30)
  grViz(gv_source, width=width, height=height)
}

#### Graph counts and plots ####

# helper function for trace plots over MCMC iterations and autocorrelation plots
iteration_plot = function(vector, ylabel, start=0, end=1, type="trace", acf_window=100, ...) {
  x = (start*length(vector)+1):(end*length(vector))
  if (type=="trace") plot(x, vector[x], type='l', ylab=ylabel, xlab='Iterations')
  if (type=="acf") {
    acf(vector[x], acf_window, main="")
    axis(4, at=seq(0,1,by=.1))
  }
}

##' @title Log posterior  plot
##' @param sbfc_result An object of class \code{sbfc}.
##' @param start The start of the included range of MCMC iterations (default=0, i.e. starting with the first iteration).
##' @param end The end of the included range of MCMC iterations (default=1, i.e. ending with the last iteration).
##' @param type Type of plot (either \code{trace} or \code{acf}, default=\code{trace}).
##' @description Plots the log posterior for a range of the MCMC iterations (indicated by \code{start} and \code{end}).
##' @export
logposterior_plot = function(sbfc_result, start=0, end=1, type="trace") {
  iteration_plot(sbfc_result$logposterior, "Log posterior", start, end, type)
}

# frequency matrix for graph edges
freq_matrix = function(sbfc_result, ...) {
  parents = sbfc_result$parents
  nvar = ncol(parents)
  corr = matrix(0, nvar, nvar)
  for (j in 1:ncol(parents)) {
    anc = sort(unique(parents[,j])) + 1
    freq_edge = table(parents[,j])/nrow(parents)
    for (k in 1:length(freq_edge)) {
      i = anc[k]
      if (i > 0) {
        corr[j, i] = corr[j, i] + freq_edge[k]
        corr[i, j] = corr[i, j] + freq_edge[k]
      }
    }
  }
  corr
}

##' @title Trace plot of Group 1 size
##' @param sbfc_result An object of class \code{sbfc}.
##' @param start The start of the included range of MCMC iterations (default=0, i.e. starting with the first iteration).
##' @param end The end of the included range of MCMC iterations (default=1, i.e. ending with the last iteration).
##' @param samples Calculate signal group size based on sampled MCMC graphs after burn-in and thinning,
##' rather than graphs from all iterations (default=FALSE).
##' @description Plots the Group 1 size for a range of the MCMC iterations (indicated by \code{start} and \code{end}).
##' @export
signal_size_plot = function(sbfc_result, start=0, end=1, samples=F) {
  n = nrow(sbfc_result$groups)
  if (samples) rows = seq(n/sbfc_result$burnin_denom + 1, n, by=sbfc_result$thin)
  else rows = 1:n
  ss = apply(sbfc_result$groups[rows,], 1, sum)
  iteration_plot(ss, "signal group size", start, end)
}

##' @title Signal variable proportion
##' @param sbfc_result An object of class \code{sbfc}.
##' @param nvars Number of top signal variables to include in the plot (default=10).
##' @param samples Calculate signal variable proportion based on sampled MCMC graphs after burn-in and thinning,
##' rather than graphs from all iterations (default=FALSE).
##' @param labels A vector of node labels (default=\code{c("X1","X2",...)}).
##' @param label_size Size of variable labels on the X-axis (default=1).
##' @param rotate_labels Rotate x-axis labels by 90 degrees to make them vertical (default=FALSE)
##' @description For each variable, computes the proportion of the samples in which this variable is in the signal group (Group 1). 
##' Plots the top \code{nvars} variables in decreasing order of signal proportion.
##' @return Signal proportion for the top \code{nvars} variables in decreasing order.
##' @export
signal_var_proportion = function(sbfc_result, nvars=10, samples=F, 
                                 labels=paste0("X", 1:ncol(sbfc_result$parents)), label_size=1, rotate_labels=F) {
  n = nrow(sbfc_result$groups)
  if (samples) rows = seq(n/sbfc_result$burnin_denom + 1, n, by=sbfc_result$thin)
  else rows = 1:n
  sig_prop = apply(sbfc_result$groups[rows,], 2, mean)
  names(sig_prop) = labels
  sort_prop = sort(sig_prop, decreasing=T)
  barplot(sort_prop[1:nvars], cex.names = label_size, cex.lab=1, ylab="Group 1 proportion", las=ifelse(rotate_labels, 3, 1))
  sort_prop[1:nvars]
} 

# plots the total number of trees of different sizes that appeared in the sample graphs
tree_size_plot = function(sbfc_result, samples=F, ...) {
  n = nrow(sbfc_result$trees)
  if (samples) rows = seq(n/sbfc_result$burnin_denom + 1, n, by=sbfc_result$thin)
  else rows = 1:n
  sizes = table(unlist(apply(sbfc_result$trees[rows,], 1, table)))
  barplot(sizes, log="y", ylab="number of occurrences", xlab="tree size")
}

##' @title Plots the density of edges in a given group over the MCMC iterations
##' @param sbfc_result An object of class \code{sbfc}.
##' @param group Which group (0 or 1) to plot edge density for.
##' @param start The start of the included range of MCMC iterations (default=0, i.e. starting with the first iteration).
##' @param end The end of the included range of MCMC iterations (default=1, i.e. ending with the last iteration).
##' @description Plots the edge density for the given group for a range of the MCMC iterations (indicated by \code{start} and \code{end}).
##' @export
edge_density_plot = function(sbfc_result, group, start=0, end=1) {
  n = nrow(sbfc_result$groups)
  edge_den = seq(0, 0, length=n)
  for (i in 1:n) {
    subset = which(sbfc_result$groups[i,] == group)
    treeset = unique(sbfc_result$trees[i, subset])
    if (length(subset) > 1) edge_den[i] = 1 - (length(treeset)-1)/(length(subset)-1)
  }
  iteration_plot(edge_den, paste("Edge density in group", group), start, end)
}

# scatterplot of edge frequency between pairs of variables vs correlation between those variables in the data set
edge_freq_plot = function(sbfc_result, data, ...) {
  corr = cor(data$TrainX)-diag(ncol(data$TrainX))
  freq = freq_matrix(sbfc_result, ...)
  plot(abs(corr), freq)
}
