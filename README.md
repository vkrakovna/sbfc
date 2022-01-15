# sbfc

Selective Bayesian Forest Classifier - R package. 
Paper: Interpretable Selection and Visualization of Features and Interactions Using Bayesian Forests (http://arxiv.org/abs/1506.02371)

It is becoming increasingly important for machine learning methods to make predictions that are interpretable as well as accurate. In many practical applications, it is of interest which features and feature interactions are relevant to the prediction task. We present a novel method, Selective Bayesian Forest Classifier, that strikes a balance between predictive power and interpretability by simultaneously performing classification, feature selection, feature interaction detection and visualization. It builds parsimonious yet flexible models using tree-structured Bayesian networks, and samples an ensemble of such models using Markov chain Monte Carlo. We build in feature selection by dividing the trees into two groups according to their relevance to the outcome of interest. Our method performs competitively on classification and feature selection benchmarks in low and high dimensions, and includes a visualization tool that provides insight into relevant features and interactions. 

To install SBFC, run install.packages() on CRAN, or download the R source package sbfc-*.tar.gz, and run
install.packages(path_to_source_package, repos = NULL, type="source")
