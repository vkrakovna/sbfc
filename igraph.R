require('igraph')


g <- make_empty_graph(n = 5) 
g=  add_edges(g, c(1,2, 2,3, 3,4, 4,5), color="red") 
g=  add_edges(g, c(5,1), color = "green")
plot(g)



single_sbfc_graph = function(groups, parents, i, noise_singletons=F, names=paste0("X", 1:ncol(parents))) {
  #g <- make_empty_graph(n = ncol(parents), with_vertex_(name = names))
  
  g = make_empty_graph()
  for (j in 1:ncol(parents)) {
    if(groups[i, j] == 1)
      g = g + vertex(names[j], fontcolor="white", color="darkblue")
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