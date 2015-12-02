##' @useDynLib sbfc
##' @importFrom Rcpp evalCpp
#### SBFC graphs ####

# produces GraphViz code for a graph for a single MCMC sample
single_sbfc_graph = function(groups, parents, i, samples=F, thin=50, noise_singletons=F, names=paste0("X", 1:ncol(parents))) {
  if (samples) i = i/thin
  s = 'strict graph G { node[fontname=default shape=circle]'
  for (j in 1:ncol(parents)) {
    if(groups[i, j] == 1)
      s = paste(s, names[j], "[fontcolor=white fillcolor=dodgerblue3];")
  }
  s = paste(s, "node[fontcolor=black fillcolor=cadetblue1]")
  for (j in 1:ncol(parents)) {
    if((groups[i, j] == 0) && noise_singletons)
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
average_sbfc_graph = function(groups, parents, cutoff=0.2, edges_only=F, noise_singletons=F,
                              names=paste0("X", 1:ncol(parents)), ...) {
  ncolors=5
  freq.group1 = apply((groups >= 1 & groups < 3), 2, mean)
  freq.group0 = apply((groups == 0 | groups == 3), 2, mean)
  ae = average_sbfc_graph_edges(parents, cutoff, names = names)
  vars = 1:ncol(parents)
  if (edges_only) vars = ae$edge_nodes
  else if (!noise_singletons) vars = unique(c(ae$edge_nodes, which(freq.group1 >= .2)))
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
  requireNamespace(Rgraphviz)
  file = tempfile()
  write(gv_source, file)
  plot(agread(file))
}

##' @title SBFC graph
##' @description Plots a sampled MCMC graph or an average of sampled graphs.
##' @param sbfc_result An object of class \code{sbfc}
##' @param average Whether to plot an average of sampled MCMC graphs (default=TRUE)
##' @param iter MCMC iteration of the sampled graph, if \code{average=F} (default=10000)
##' @param edge_cutoff The average graph includes edges that appear in at least this fraction of the sampled graphs, if \code{average=T} (default=0.2)
##' @param noise_singletons Whether to plot single-node trees in the noise group (Group 0), which can be numerous for high-dimensional data sets (default=FALSE).
##' @param names Node names for the graph labels (default=\code{c("X1","X2",...)})
##' @details If the graph renders poorly in the R plot window, try changing the aspect ratio of the plot window.
##' @export
sbfc_graph = function(sbfc_result, iter=10000, average=T, edge_cutoff=0.2, edges_only=F, ...) {
  parents = sbfc_result$parents
  groups = sbfc_result$groups
  if (average) gv_source = average_sbfc_graph(groups, parents, edge_cutoff, edges_only, ...)
  else gv_source = single_sbfc_graph(groups, parents, iter, ...)
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
##' @description Plots the log posterior for a range of the MCMC iterations (indicated by \code{start} and \code{end}).
##' @export
logpost_plot = function(sbfc_result, start=0, end=1, ...) {
  iteration_plot(sbfc_result$logpost, "logpost", start, end, ...)
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

# trace plot of Group 1 size
signal_size_plot = function(sbfc_result, start=0, end=1, subset=F, thin=50, ...) {
  n = nrow(sbfc_result$groups)
  if (subset) rows = seq(n/5 + 1, n, by=thin)
  else rows = 1:n
  ss = apply(sbfc_result$groups[rows,], 1, sum)
  iteration_plot(ss, "signal group size", start, end, ...)
}

##' @title Signal variable proportion
##' @description Computes proportion of the samples in which each variable is in the signal group (Group 1). Plots the first (nvars) variables in decreasing order of proportion.
##' @export
signal_var_prop = function(sbfc_result, subset=F, thin=50, nvars=10, label_size=1) {
  n = nrow(sbfc_result$groups)
  if (subset) rows = seq(n/5 + 1, n, by=thin)
  else rows = 1:n
  sig_prop = apply(sbfc_result$groups[rows,], 2, mean)
  names(sig_prop) = paste0("X", 1:ncol(sbfc_result$groups))
  sort_prop = sort(sig_prop, decreasing=T)
  barplot(sort_prop[1:nvars], cex.names = label_size, cex.lab=1.5, ylab="Group 1 proportion", xlab="variable")
  sort_prop[1:nvars]
} 

# plots the total number of trees of different sizes that appeared in the sample graphs
tree_size_plot = function(sbfc_result, subset=F, thin=50, ...) {
  n = nrow(sbfc_result$trees)
  if (subset) rows = seq(n/5 + 1, n, by=thin) # without burnin
  else rows = 1:n
  sizes = table(unlist(apply(sbfc_result$trees[rows,], 1, table)))
  barplot(sizes, log="y", ylab="number of occurrences", xlab="tree size")
}

# plots density of edges in a given group over the MCMC iterations
edge_density = function(sbfc_result, group, start=0, end=1, subset=F, thin=50, ...) {
  n = nrow(sbfc_result$groups)
  edge_den = seq(0, 0, length=n)
  for (i in 1:n) {
    subset = which(sbfc_result$groups[i,] == group)
    treeset = unique(t(sbfc_result$trees[i, subset]))
    if (length(subset) > 1) edge_den[i] = 1 - (length(treeset)-1)/(length(subset)-1)
  }
  iteration_plot(edge_den, paste("edge density in group", group), start, end, ...)
}

# scatterplot of edge frequency between pairs of variables vs correlation between those variables in the data set
edge_freq_plot = function(sbfc_result, data, ...) {
  corr = cor(data$TrainX)-diag(ncol(data$TrainX))
  freq = freq_matrix(sbfc_result, ...)
  plot(abs(corr), freq)
}
