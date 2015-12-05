##' @useDynLib sbfc
##' @importFrom Rcpp evalCpp

#### SBFC graphs ####

# produces GraphViz code for a graph for a single MCMC sample
single_sbfc_graph = function(groups, parents, i, single_noise_nodes=F, names=paste0("X", 1:ncol(parents))) {
  s = 'strict graph G { node[fontname=default shape=circle]'
  s = paste(s, "node[fontcolor=white fillcolor=dodgerblue3]")
  for (j in 1:ncol(parents)) {
    if(groups[i, j] == 1)
      s = paste(s, names[j], ";")
  }
  s = paste(s, "node[fontcolor=black fillcolor=cadetblue1]")
  for (j in 1:ncol(parents)) {
    if((groups[i, j] == 0) && single_noise_nodes)
      s = paste(s, names[j], ";")
  }
  for (j in 1:ncol(parents)) {
    if(parents[i, j] != 0)
      s = paste(s, names[parents[i, j]], "--", names[j], ";")
  }
  s = paste(s, "}")
  s
}

# determines a set of edges to include in the average graph
average_sbfc_graph_edges = function(parents, cutoff=0.2, names=paste0("X", 1:ncol(parents))) {
  s =""
  edge_nodes = c()
  for (j in 1:ncol(parents)) {
    par = sort(unique(parents[,j]))
    freq.edge = table(parents[,j])/nrow(parents)
    for (k in 1:length(freq.edge)) {
      if (par[k] > 0 && freq.edge[k] >= cutoff) {
        s = paste(s, names[par[k]], "--", names[j], ";")
        edge_nodes = c(edge_nodes, j, par[k])
      }
    }
  }
  list(s = s, edge_nodes = unique(edge_nodes))
}

# produces GraphViz code for an average graph over a set of MCMC sample graphs
average_sbfc_graph = function(groups, parents, cutoff=0.2, single_noise_nodes=F,
                              names=paste0("X", 1:ncol(parents))) {
  ncolors=5
  freq.group1 = apply((groups >= 1 & groups < 3), 2, mean)
  freq.group0 = apply((groups == 0 | groups == 3), 2, mean)
  ae = average_sbfc_graph_edges(parents, cutoff, names = names)
  vars = 1:ncol(parents)
  #if (edges_only) vars = ae$edge_nodes
  if (!single_noise_nodes) vars = unique(c(ae$edge_nodes, which(freq.group1 >= .2)))
  fontcolor=c(rep("black", floor(ncolors/2)), rep("white", ceiling(ncolors/2)))
  
  col = c("aliceblue", "cadetblue1", "deepskyblue", "dodgerblue3", "dodgerblue4")
  s = "strict graph G { node[fontname=default shape=circle] "
  for (i in ncolors:1) {
    for (j in vars) {
      if ((freq.group1[j] >= (i-1)*1.0/ncolors) && (freq.group1[j] <= i*1.0/ncolors)) 
        s = paste(s, names[j], "[fontcolor=", fontcolor[i], "fillcolor=", col[i], "];")
    }
  }
  s = paste(s, ae$s, "}")
  s
}

graphviz_plot = function(gv_source) {
  require('Rgraphviz')
  requireNamespace('Rgraphviz')
  file = tempfile()
  write(gv_source, file)
  plot(agread(file))
}

##' @title SBFC graph
##' @description Plots a sampled MCMC graph or an average of sampled graphs.
##' @param sbfc_result An object of class \code{sbfc}.
##' @param average Plot an average of sampled MCMC graphs (default=TRUE).
##' @param iter MCMC iteration of the sampled graph, if \code{average=F} (default=10000).
##' @param edge_cutoff The average graph includes edges that appear in at least this fraction of the sampled graphs, if \code{average=T} (default=0.2).
##' @param edges_only Omit single-node trees (default=FALSE).
##' @param single_noise_nodes Plot single-node trees in the noise group (Group 0), which can be numerous for high-dimensional data sets (default=FALSE).
##' @param labels Node labels (default=\code{c("X1","X2",...)}).
##' @details If the graph renders poorly in the R plot window, try changing the aspect ratio of the plot window.
##' @export
sbfc_graph = function(sbfc_result, iter=10000, average=T, edge_cutoff=0.2, single_noise_nodes=F, labels=paste0("X", 1:ncol(sbfc_result$parents))) {
  parents = sbfc_result$parents
  groups = sbfc_result$groups
  if (average) gv_source = average_sbfc_graph(groups, parents, edge_cutoff, single_noise_nodes, labels)
  else gv_source = single_sbfc_graph(groups, parents, iter, single_noise_nodes, labels)
  graphviz_plot(gv_source)
}

#### Graph counts and plots ####

# helper function for trace plots over MCMC iterations and autocorrelation plots
iteration_plot = function(vector, ylabel, start = 0, end = 1, type="trace", acf_window=100, ...) {
  x = (start*length(vector)+1):(end*length(vector))
  if (type=="trace") plot(x, vector[x], type='l', ylab=ylabel)
  if (type=="acf") {
    acf(vector[x], acf_window, main="")
    axis(4, at=seq(0,1,by=.1))
  }
}

##' @title Log posterior trace plot
##' @param sbfc_result An object of class \code{sbfc}.
##' @param start The start of the included range of MCMC iterations (default=0, i.e. starting with the first iteration).
##' @param end The end of the included range of MCMC iterations (default=1, i.e. ending with the last iteration).
##' @description Plots the log posterior for a range of the MCMC iterations (indicated by \code{start} and \code{end}).
##' @export
logposterior_plot = function(sbfc_result, start=0, end=1, ...) {
  iteration_plot(sbfc_result$logposterior, "Log posterior", start, end, ...)
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
signal_size_plot = function(sbfc_result, start=0, end=1, samples=F, ...) {
  n = nrow(sbfc_result$groups)
  if (samples) rows = seq(n/sbfc_result$burnin_denom + 1, n, by=sbfc_result$thin)
  else rows = 1:n
  ss = apply(sbfc_result$groups[rows,], 1, sum)
  iteration_plot(ss, "signal group size", start, end, ...)
}

##' @title Signal variable proportion
##' @param sbfc_result An object of class \code{sbfc}.
##' @param nvars Number of top signal variables to include in the plot (default=10).
##' @param samples Calculate signal variable proportion based on sampled MCMC graphs after burn-in and thinning,
##' rather than graphs from all iterations (default=FALSE).
##' @param label_size Size of variable labels on the X-axis (default=1).
##' @description For each variable, computes the proportion of the samples in which this variable is in the signal group (Group 1). 
##' Plots the top \code{nvars} variables in decreasing order of signal proportion.
##' @return Signal proportion for the top \code{nvars} variables in decreasing order.
##' @export
signal_var_proportion = function(sbfc_result, nvars=10, samples=F, label_size=1) {
  n = nrow(sbfc_result$groups)
  if (samples) rows = seq(n/sbfc_result$burnin_denom + 1, n, by=sbfc_result$thin)
  else rows = 1:n
  sig_prop = apply(sbfc_result$groups[rows,], 2, mean)
  names(sig_prop) = paste0("X", 1:ncol(sbfc_result$groups))
  sort_prop = sort(sig_prop, decreasing=T)
  barplot(sort_prop[1:nvars], cex.names = label_size, cex.lab=1.5, ylab="Group 1 proportion", xlab="Variable")
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
edge_density_plot = function(sbfc_result, group, start=0, end=1, ...) {
  n = nrow(sbfc_result$groups)
  edge_den = seq(0, 0, length=n)
  for (i in 1:n) {
    subset = which(sbfc_result$groups[i,] == group)
    treeset = unique(sbfc_result$trees[i, subset])
    if (length(subset) > 1) edge_den[i] = 1 - (length(treeset)-1)/(length(subset)-1)
  }
  iteration_plot(edge_den, paste("Edge density in group", group), start, end, ...)
}

# scatterplot of edge frequency between pairs of variables vs correlation between those variables in the data set
edge_freq_plot = function(sbfc_result, data, ...) {
  corr = cor(data$TrainX)-diag(ncol(data$TrainX))
  freq = freq_matrix(sbfc_result, ...)
  plot(abs(corr), freq)
}
