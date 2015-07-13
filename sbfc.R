########## I/O helper functions ###########

paste0 = function(...) {
  paste(..., sep="")
}

gen_path = function(output_type, id, dir, pref, samples=F, fold=0, filetype="txt", ...) {
  path = paste0(pref, id, "/", dir, "/")
  if (samples) suff = ""
  else suff = "_all"
  if (fold > 0) fold_str = paste0("_", fold)
  else fold_str=""
  paste(path, id, paste("_", dir, sep=""), fold_str, "_", output_type, suff, ".", filetype, sep = "")
}

data_import = function(id, pref, nvar=0, test=T, factor=F, ...) {
  path = paste0(pref, id, "/")
  data = list()
  data$TrainX = read.table(paste0(path, "TrainX"))
  if (nvar > 0) data$TrainX = data$TrainX[,1:nvar]
  data$TrainY = as.factor(read.table(paste0(path, "TrainY"))[,1])
  if (test) {
    data$TestX = read.table(paste0(path, "TestX"))
    if (nvar > 0) data$TestX = data$TestX[,1:nvar]
    data$TestY = factor(read.table(paste0(path, "TestY"))[,1], levels(data$TrainY))
  }
  if (factor) {
    if (test) {
      joined = as.data.frame(lapply(rbind(data$TestX,data$TrainX), factor))
      data$TrainX = joined[(1+nrow(data$TestX)):nrow(joined),]
      data$TestX = joined[1:nrow(data$TestX),]
    } else {
      data$TrainX = as.data.frame(lapply(data$TrainX, factor))
    }
  }
  data$Train = cbind(data$TrainY, data$TrainX)
  colnames(data$Train)[1]='y'
  data
}

output_data = function(data, id, pref, test=T) {
  path = paste(pref, id, sep="") 
  dir.create(file.path(path), showWarnings = F)
  write.table(data$TrainX, file=paste(path, "/TrainX",sep=""), quote=F, col.names=F, row.names=F)
  write.table(data$TrainY, file=paste(path, "/TrainY",sep=""), quote=F, col.names=F, row.names=F)
  if (test) {
    write.table(data$TestX, file=paste(path, "/TestX",sep=""), quote=F, col.names=F, row.names=F)
    write.table(data$TestY, file=paste(path, "/TestY",sep=""), quote=F, col.names=F, row.names=F)
  }
  if (substr(id, 1, 1) < 2)   
    write.table(data$tm, file=paste(path, "/TrueModel",sep=""), quote=F, col.names=F, row.names=F)
}

################ Data structure summaries ######

## required parameters: 
# data set name/id (id), SBFC output directory (dir), directory path prefix (pref)
## required parameters for small data sets that require cross-validation: 
# which cross-validation fold to examine (fold), whether there is a test set (test)

#### Tree graphs ####

# produces GraphViz code for a graph for a single MCMC sample
tree_graph = function(groups, ancestors, i, samples=F, thin=50, ...) {
  if (samples) i = i/thin
  s = 'digraph G { subgraph cluster_g1 { node [color=blue]; label="Group 1";'
  for (j in ncol(ancestors):1) {
    if(groups[i, j] == 1)
      s = paste0(s, " X", j, ";")
    if(groups[i, j] == 1 && ancestors[i, j] != -1)
      s = paste0(s, " X", ancestors[i, j]+1, "-> X", j, ";")
  }
  s = paste(s, '} subgraph cluster_g0 { node [color=red]; label="Group 0";')
  for (j in ncol(ancestors):1) {
    if(groups[i, j] == 0)
      s = paste0(s, " X", j, ";")
    if(groups[i, j] == 0 && ancestors[i, j] != -1)
      s = paste0(s, " X", ancestors[i, j]+1, "-> X", j, ";")
  }
  s = paste(s, "}}")
  s
}

# determines a set of edges to include in the average graph
average_tree_graph_edges = function(ancestors, cutoff=0.2, names = paste0("X", 1:ncol(ancestors))) {
  s =""
  edge_nodes = c()
  for (j in 1:ncol(ancestors)) {
    anc = sort(unique(ancestors[,j])) + 1
    freq.edge = table(ancestors[,j])/nrow(ancestors)
    for (k in 1:length(freq.edge)) {
      if (anc[k] > 0 && freq.edge[k] >= cutoff) {
        s = paste0(s, " \"", names[anc[k]], "\"-> \"", names[j], "\";")
        edge_nodes = c(edge_nodes, j, anc[k])
      }
    }
  }
  list(s = s, edge_nodes = unique(edge_nodes))
}

# produces GraphViz code for an average graph over a set of MCMC sample graphs
average_tree_graph = function(groups, ancestors, cutoff=0.2, edges_only = F, noise_singletons=F,
                              names = paste0("X", 1:ncol(ancestors)), colorscheme="rdylbu", ncolors=7, ...) {
  freq.group1 = apply((groups >= 1 & groups < 3), 2, mean)
  freq.group0 = apply((groups == 0 | groups == 3), 2, mean)
  ae = average_tree_graph_edges(ancestors, cutoff, names = names)
  vars = 1:ncol(ancestors)
  if (edges_only) vars = ae$edge_nodes
  else if (!noise_singletons) vars = unique(c(ae$edge_nodes, which(freq.group1 >= .5)))
  
  s = "digraph G { "  
  
  if (colorscheme!="discrete") {
    for (i in ncolors:1) {
      for (j in vars) {
        if ((freq.group1[j] >= (i-1)*1.0/ncolors) && (freq.group1[j] <= i*1.0/ncolors)) 
          s = paste0(s, " \"", names[j], "\"[colorscheme=", colorscheme, ncolors, " color=", i, "];")
      }
    }
  } else {
    for (j in vars) {
      if ((freq.group1[j] >= cutoff) & (freq.group0[j] < cutoff)) 
        s = paste0(s, " ", names[j], "[color=blue];")
    }
    for (j in vars) {
      if ((freq.group1[j] >= cutoff) & (freq.group0[j] >= cutoff)) 
        s = paste0(s, " ", names[j], "[color=green];")
    }
    for (j in vars) {
      if ((freq.group1[j] < cutoff) & (freq.group0[j] >= cutoff)) 
        s = paste0(s, " ", names[j], "[color=red];")
    }
  }
  s = paste(s, ae$s, "}")
  s
}

# outputs GraphViz code for a graph
# required parameters: data set id (id), SBFC output directory (dir), directory path prefix (pref)
output_tree_graph = function(row=10000, average=F, cutoff=0.2, edges_only=F, ...) {
  ancestors = read.table(gen_path("Ancestors", ...), sep="")
  groups = read.table(gen_path("Groups", ...), sep="")
  if (average) {
    gv_source = average_tree_graph(groups, ancestors, cutoff, edges_only, ...)
    row = "average"
  }
  else gv_source = tree_graph(groups, ancestors, row, ...)
  
  s = paste0("graph_", row, "_cutoff_", cutoff)
  if (edges_only) s = paste0(s, "_edges_only")
  path = gen_path(s, filetype="gv", ...)
  write(gv_source, file = path)
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

# plots the log posterior 
logpost_plot = function(start=0, end=1, save=F, ...) {
  path = gen_path("LogPost", ...)
  logpost = read.table(path)[,1]
  iteration_plot(logpost, "logpost", start, end, ...)
  if (save) {
    path = gen_path(paste(type, "_plot(", start, ",", end,")",sep=""), filetype="png", ...)
    png(path)
    iteration_plot(logpost, "logpost", start, end, ...)
    dev.off()
  }
}

# frequency matrix for graph edges
freq_matrix = function(...) {
  ancestors = read.table(gen_path("Ancestors", ...), sep="")
  nvar = ncol(ancestors)
  corr = matrix(0, nvar, nvar)
  for (j in 1:ncol(ancestors)) {
    anc = sort(unique(ancestors[,j])) + 1
    freq_edge = table(ancestors[,j])/nrow(ancestors)
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
signal_size_plot = function(start=0, end=1, save=T, subset=F, thin=50, ...) {
  groups = read.table(gen_path("Groups", ...))
  n = nrow(groups)
  if (subset) rows = seq(n/5 + 1, n, by=thin)
  else rows = 1:n
  ss = apply(groups[rows,], 1, sum)
  iteration_plot(ss, "signal group size", start, end, ...)
  if (save) {
    path = gen_path(paste("signal_size(", start, ",", end,")", sep=""), filetype="png", ...)
    png(path)
    iteration_plot(ss, "signal group size", start, end, ...)
    dev.off()
  }
}

# acceptance rate for the switch move
switch_acc_rate = function(k=10, ...) {
  switch_acc = read.table(gen_path("SwitchAcc", ...))[,1]
  mean(switch_acc)/k
}

# computes proportion of the samples in which each node/variable is in Group 1
# plots the first (nvars) variables in decreasing order of proportion
signal_node_prop = function(save=T, subset=F, thin=50, nvars=20, ...) {
  groups = read.table(gen_path("Groups", ...))
  n = nrow(groups)
  if (subset) rows = seq(n/5 + 1, n, by=thin)
  else rows = 1:n
  sig_prop = apply(groups[rows,], 2, mean)
  sort_prop = sort(sig_prop, decreasing=T)
  barplot(sort_prop[1:nvars], cex.names = .5, ylab="Group 1 proportion", xlab="variable")
  
  if (save) {
    path = gen_path("signal_prop", filetype="png", ...)
    png(path)
    barplot(sort_prop[1:nvars], cex.names = .5, ylab = "Group 1 proportion", xlab="variable")
    dev.off()
  }
  sort_prop[1:nvars]
} 

# plots the total number of trees of different sizes that appeared in the sample graphs
tree_size_plot = function(save=T, subset=F, thin=50, ...) {
  trees = read.table(gen_path("Trees", ...))
  n = nrow(trees)
  if (subset) rows = seq(n/5 + 1, n, by=thin) # without burnin
  else rows = 1:n
  sizes = table(unlist(apply(trees[rows,], 1, table)))
  barplot(sizes, log="y", ylab="number of occurrences", xlab="tree size")
  if (save) {
    path = gen_path("tree_sizes", filetype="png", ...)
    png(path)
    barplot(sizes, log="y", ylab="number of occurrences", xlab="tree size")
    dev.off()
  }
}

# plots density of edges in a given group over the MCMC iterations
edge_density = function(group, start=0, end=1, save=T, subset=F, thin=50, ...) {
  groups = read.table(gen_path("Groups", ...))
  trees = read.table(gen_path("Trees", ...))
  n = nrow(groups)
  edge_den = seq(0, 0, length=n)
  for (i in 1:n) {
    subset = which(groups[i,] == group)
    treeset = unique(t(trees[i, subset]))
    if (length(subset) > 1) edge_den[i] = 1 - (length(treeset)-1)/(length(subset)-1)
  }
  iteration_plot(edge_den, paste("edge density in group", group), start, end, ...)
  if (save) {
    path = gen_path(paste0("edge_density", group, "(", start, ",", end,")"), filetype="png", ...)
    png(path)
    iteration_plot(edge_den, paste("edge density in group", group), start, end, ...)
    dev.off()
  } 
}

# scatterplot of edge frequency between pairs of variables vs correlation between those variables in the data set
edge_freq_plot = function(save=T, ...) {
  data = data_import(...)
  corr = cor(data$TrainX)-diag(ncol(data$TrainX))
  freq = freq_matrix(...)
  plot(abs(corr), freq)
  if (save) {
    path = gen_path("edge_freq_plot", filetype="png", ...)
    png(path)
    plot(abs(corr), freq)
    dev.off()
  }
}

################ Other classification methods for comparison ###############################

# naive Bayes
nb = function(data, e1071=T, smooth=F) {  
  t = proc.time()
  testcategory = if(e1071) {
    require("e1071")
    laplace=0
    if (smooth) laplace = mean(1/sapply(data$TrainX,nlevels))/nlevels(data$TrainY)
    predict(naiveBayes(y~.,data=data$Train,laplace=laplace), data$TestX)
  } else {
    base.rates = as.numeric(1+table(data$TrainY))
    clas = matrix(log(base.rates),nrow(data$TestX),nlevels(data$TrainY), byrow=T)
    for(i in 1:ncol(data$TrainX)) {
      clas = clas + log((table(data$TrainX[,i], data$TrainY) + 
                           1/nlevels(data$TrainX[,i])/nlevels(data$TrainY))[data$TestX[,i],] /
                        base.rates[col(clas)])
    }
    levels(data$TrainY)[apply(clas,1,which.max)]
  }
  accuracy = mean(testcategory==data$TestY)
  runtime = proc.time() - t
  list(accuracy = accuracy, runtime = runtime, testcategory=testcategory)
}

# Random Forest
rf = function(data) {  
  require(randomForest)
  t = proc.time()
  rf = randomForest(y~.,data=data$Train, importance=T)
  imp_mat = rf$importance
  imp = imp_mat[,ncol(imp_mat)-1]
  testcategory = predict(rf, data$TestX)
  #levels(testcategory) = levels(y)
  accuracy = mean(testcategory==data$TestY)
  runtime = proc.time() - t
  list(accuracy = accuracy, runtime = runtime, testcategory=testcategory, importance=sort(imp, decreasing=T))
  
}

# CART
ct = function(data) { 
  require(tree)
  t = proc.time()
  testcategory = predict(tree(y~.,data=data$Train),data$TestX, type="class")
  accuracy = mean(testcategory==data$TestY)
  runtime = proc.time() - t
  list(accuracy = accuracy, runtime = runtime, testcategory=testcategory)
}

# BART
bt = function(data) { 
  require(BayesTree)
  t = proc.time()
  bart.res = bart(x.train = data$TrainX * 1.0, y.train = data$TrainY, x.test = data$TestX, verbose=F)
  testcategory = as.factor(round(apply(pnorm(bart.res$yhat.test), 2, mean)))
  levels(testcategory) = levels(data$TestY)
  accuracy = mean(testcategory==data$TestY)
  runtime = proc.time() - t
  list(accuracy = accuracy, runtime = runtime, testcategory=testcategory)
}

# logistic regression
lr = function(data) { 
  t = proc.time()
  testcategory = as.factor(round(predict(glm(y~.,family=binomial,data$Train),data$TestX,type="response")))
  levels(testcategory) = levels(data$TestY)
  accuracy = mean(testcategory==data$TestY)
  runtime = proc.time() - t
  list(accuracy = accuracy, runtime = runtime, testcategory=testcategory)
}

# Lasso (glmnet)
ls = function(data) { 
  require('glmnet')
  t = proc.time()
  model = cv.glmnet(x=as.matrix(data$TrainX), y=data$TrainY, family="multinomial")
  testcategory = as.numeric(predict(model, as.matrix(data$TestX), type="class"))
  accuracy = mean(testcategory==data$TestY)
  runtime = proc.time() - t
  list(accuracy = accuracy, runtime = runtime, testcategory=testcategory)
}

# Support Vector Machines
sv = function(data){  
  require("e1071")
  t = proc.time()
  testcategory = predict(svm(y~.,data=data$Train), data$TestX)
  accuracy = mean(testcategory==data$TestY)
  runtime = proc.time() - t
  list(accuracy = accuracy, runtime = runtime, testcategory=testcategory)
}

# C.50 (decision tree classifier)
c5 = function(data) { 
  require("C50")
  t = proc.time()
  testcategory = predict(C5.0(y~.,data=data$Train), newdata=data$TestX)
  accuracy = mean(testcategory==data$TestY)
  runtime = proc.time() - t
  list(accuracy = accuracy, runtime = runtime, testcategory=testcategory)
}

# run a classification method with cross-validation (for small data sets)
run_method_cv = function(data, method, n_folds=5) {
  t = proc.time()
  n_units = length(data$TrainY)
  row_shuffle = sample(1:n_units, replace=F)
  testcategory = as.factor(rep(NA, n_units))
  levels(testcategory) = levels(data$TrainY)
  fold_start = 1
  fold_size = floor(n_units / n_folds)
  mod = n_units %% n_folds
  
  for (i in 1:n_folds) {
    fold_end = fold_start + fold_size
    if (i > mod) fold_end = fold_end - 1
    test_subset = row_shuffle[fold_start:fold_end]
    train_subset = row_shuffle[-(fold_start:fold_end)]
    
    datafold = list()
    datafold$TrainX = data$TrainX[train_subset,]
    datafold$TestX = data$TrainX[test_subset,] 
    datafold$TrainY = data$TrainY[train_subset]
    datafold$TestY = data$TrainY[test_subset]
    datafold$Train = cbind(datafold$TrainY, datafold$TrainX)
    colnames(datafold$Train)[1]='y'
    
    testcategory[test_subset] = method(datafold)$testcategory
    fold_start = fold_end + 1
  }
  stopifnot(fold_end == n_units)
  levels(testcategory) = levels(data$TrainY)
  accuracy = mean(testcategory==data$TrainY)
  runtime = proc.time() - t
  list(accuracy = accuracy, testcategory = testcategory, runtime = runtime)
}

################# Results ##########################################

#### Helper functions #####

list_dirs <- function(path=".", pattern=NULL, all.dirs=FALSE, full.names=FALSE, ignore.case=FALSE, dirs_only = TRUE) {
  # code from http://stackoverflow.com/questions/4749783/how-to-obtain-a-list-of-directories-within-a-directory-like-list-files-but-i
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
                    full.names=TRUE, recursive=FALSE, ignore.case)
  if (dirs_only) {
    dirs <- all[file.info(all)$isdir]
  } else {
    dirs <- all
  }
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}

strip = function(out) {
  as.numeric(substr(out, 1, nchar(out)-1))
}

info_import = function(type="Accuracy", ...) {
  output = scan(gen_path("Output", ...), "raw", quiet=T)
  acc_rows = grep(type, output)+1
  if (type=="Time") {
    strip(output[acc_rows])
  } else {
    as.numeric(output[acc_rows])
  }
}

# adding extra noise variables
data_augment = function(id, pref, nvar=10000, test=F) {
  alldata = data_import(id, pref, test=test)
  data = alldata$TrainX
  for (i in (ncol(data)+1):nvar) {
    index = sample(1:ncol(data), 1)
    rows = sample(1:nrow(data), nrow(data), replace=F)
    data[,i] = data[rows, index]
  }
  alldata$TrainX = data
  output_data(alldata, paste0(id, "_extra"), pref, test)
}

#### Summarization functions ####

# summarize SBFC classification results
info_sbfc = function(id, pref, suff="", nosuff="", ...) {
  dirs = list_dirs(paste0(pref, id))
  rows = grep("sbfc", dirs)
  if (suff != "") {
    suff_rows = grep(paste0("_", suff), dirs)
    rows = intersect(rows, suff_rows)
  }
  if (nosuff != "") {
    nosuff_rows = grep(paste0("_", nosuff), dirs, invert=T)
    rows = intersect(rows, nosuff_rows)
  }
  x = c()
  for (i in rows) {
    x = c(x, info_import(id=id, dir=dirs[i], samples=T, pref=pref, ...))
  }
  x
}

# summarize classification results for a given method
info_method = function(id, method_name, pref, type="Accuracy") {
  files = list_dirs(paste0(pref, id), dirs_only = F)
  rows = intersect(grep(method_name, files), grep(type, files))
  info = c()
  for (i in rows) {
    info = c(info, read.table(paste0(path1, "/", files[i]), header=F)[1,1])
  }
  info
}

# output a confidence interval for an info vector
summ = function(info) {
  paste0(signif(mean(info),3), "+-", signif(sd(info),1))
}

#### Result tables for UCI data sets ####

# for small data sets with 2 classes
method_name_list = c("bart", "rf", "nb", "cart", "c5.0", "lr", "svm", "lasso")
method_list = c(bt, rf, nb, ct, c5, lr, sv, ls)
data_list = c("australian","breast","chess","cleve","corral","crx","diabetes","flare","german","glass2","heart","hepatitis","mofn","pima","vote")
cv_list = c(T, T, F, T, T, T, T, T, T, T, T, T, F, T, T) # whether a data set requires cross-validation
# for small multiclass data sets (with >2 classes)
method_multi_list = c(rf, nb, ct, c5, sv) # methods that can handle multi-class data
method_name_multi_list = c("rf", "nb", "cart", "c5.0", "svm")
data_multi_list = c("glass", "iris", "letter", "lymphography", "satimage", "segment", "shuttle", "soybean", "vehicle", "waveform")
cv_multi_list = c(T, T, F, T, F, F, F, T, T, F)
# large data sets with >300 variables
large_data_list = c("microsoft", "madelon", "isolet", "ad", "gisette", "arcene")

# produce result table for a set of other methods on a collection of data sets
results_method_data = function(data.list, cv.list, method.list, method.name.list, pref, ...) {
  results = matrix(0, length(data.list), length(method.name.list))
  dimnames(results) = list(data.list, method.name.list)
  for (i in 1:length(data.list)) {
    id = data.list[i]
    data = data_import(id, pref, test = !cv.list[i])
    data_nb = data_import(id, pref, factor=T, test = !cv.list[i])
    for (j in 1:length(method.list)) {
      method = method.list[[j]]
      if (method.name.list[j] == "nb") { 
        d = data_nb
      } else {
        d = data
      }
      if (cv.list[i]) {
        results[i,j] = run_method_cv(d, method, ...)$accuracy
      } else {
        results[i,j] = run_method(d)$accuracy
      }
    }
  }
  results
}

# produce result table for a set of other methods on one data set with different numbers of variables
# requires the data set to include a large number of extra noise variables
results_method_nvar = function(id, cv, method.list, method.name.list, nvar= c(100,1000,10000), pref, ...) {
  results = matrix(0, length(nvar), length(method.list))
  dimnames(results) = list(nvar, method.name.list)
  for (i in 1:length(nvar)) {
    data = data_import(id, pref, nvar=nvar[i], test = !cv)
    data_nb = data_import(id, pref, nvar=nvar[i], factor=T, test = !cv)
    for (j in 1:length(method.list)) {
      method = method.list[[j]]
      if (method.name.list[j] == "nb") { 
        d = data_nb
      } else {
        d = data
      }
      if (cv) {
        results[i,j] = run_method_cv(d, method, ...)$accuracy
      } else {
        results[i,j] = method(d)$accuracy
      }
    }
  }
  results
}

