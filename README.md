# sbfc
Bayesian algorithm for simultaneous feature selection and classification. See paper: http://arxiv.org/abs/1506.02371

Feature selection and classification are fundamental tasks in machine learning that are related yet usually achieved separately. We propose a Bayesian method that strikes a balance between predictive power and interpretability by simultaneously performing classification, feature selection and feature interaction detection. We build a correlation structure on top of Naive Bayes, and introduce a latent feature to partition the features into two groups according to their relationships with the outcome of interest. In order to achieve both model flexibility and parsimony, we use trees to approximate the dependence relationships among the features, and set a complexity-penalizing prior on the tree structure parameters. We use Markov chain Monte Carlo to explore the partition and forest structure space, and combine the predictions using Bayesian model averaging. Our method performs competitively with state-of-the-art classifiers on low- and high-dimensional data sets, and provides insight into relevant features and feature interactions, complete with a visualization tool.

Files:

- sbfc.cpp: C++ code for running the SBFC algorithm

- sbfc.R: R code for result summarization and visualization
