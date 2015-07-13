# sbfc
Bayesian algorithm for simultaneous variable selection and classification. See paper: http://arxiv.org/pdf/1506.02371v1

Variable selection and classification are fundamental tasks in machine learning that are related yet usually achieved separately.
We propose a Bayesian method that strikes a balance between predictive power and interpretability by simultaneously performing 
classification, variable selection and variable interaction detection. We build a correlation structure on top of naive Bayes,
and introduce a latent variable to partition the variables into two groups according to their relationships with the outcome of 
interest. In order to achieve both model flexibility and parsimony, we use trees to approximate the dependence relationships 
among the variables, and set a complexity-penalizing prior on the tree structure parameters. We use Markov chain Monte Carlo to 
explore the partition and forest structure space, and combine the predictions using Bayesian model averaging. Our method performs
competitively with state-of-the-art classifiers on low- and high-dimensional data sets, and provides insight into relevant 
variables and variable interactions, complete with a visualization tool.

Files:
sbfc.cpp: C++ code for running the SBFC algorithm
sbfc.R: R code for result summarization and visualization
