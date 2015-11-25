// [[Rcpp::depends(RcppArmadillo)]]
//#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <armadillo>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <cassert>
#include <sys/timeb.h>
#include <sstream>
#include <cmath>

using namespace arma;
using namespace std;
using namespace Rcpp;
#define to_vec conv_to<vec>::from
#define to_svec conv_to<svec>::from
#define to_uvec conv_to<uvec>::from
#define to_smat conv_to<smat>::from

const unsigned null_value = 65535;
const double cutoff_equal = 1e-6;

typedef std::vector<unsigned> stdvec;
typedef Cube<unsigned short> scube;
typedef Mat<unsigned short> smat;
typedef Col<unsigned short> svec;

///// Data structures

enum group {
	noise = 0,
	sig = 1
};

struct graph {
	svec Group;
	svec Tree;
	svec Parent;
	graph(): Group(), Tree(), Parent() {};
	graph(int n_var): Group(n_var), Tree(n_var), Parent(n_var) {};
};

struct data {
	smat X_train;
	svec Y_train;
	smat X_test;
	svec Y_test;
	smat X;
	svec Y;
	imat true_model;
};

struct counts {
	field<scube> var_var_y; // Counts for triples of var, parent of var, and y
	field<smat> var_var; // Counts for pairs of variables
	field<smat> var_y; // Counts for pairs of var and y
	field<svec> var; // Counts for var
	svec y; // Counts for y

	counts(unsigned n_var): var_var_y(n_var, n_var), var_var(n_var, n_var), var_y(n_var), var(n_var) {};
};

struct nlevels {
	svec x;
	unsigned y;
};

struct parameters {
	// model parameters that are kept constant
	double alpha; // Dirichlet hyperparameter
	double scaling; // scaling factor
	double edge_mult; // multiplier for the x-edge penalty
	double yedge_mult; // multiplier for the y-edge penalty
	unsigned k; // number of repeats for the Switch move
	unsigned n_var; // number of variables
	unsigned n_units; // number of units
	unsigned n_step; // total number of MCMC steps including burnin
	unsigned n_rows; // number of rows in output matrices
	unsigned burnin_frac; // (denominator of) fraction of total MCMC steps discarded as burnin
	unsigned thin; // thinning factor for MCMC samples
	unsigned n_samples; // number of MCMC samples
	bool thin_output; // whether to output only the thinned samples or all samples from the MCMC
	unsigned n_folds; // number of cross-validation folds
	string start; // starting graph for the algorithm (single noise nodes, random groups with single nodes, or random trees)
	string output_id; // name of output directory
	bool thread; // whether to run cross-validation with parallel threads

	parameters(): alpha(5), edge_mult(4), yedge_mult(1), k(10), n_var(0), n_step(0), burnin_frac(5),
				 thin(50), thin_output(false), n_folds(5), start("noise"), thread(false) {};
};

struct outputs {
	smat Groups;
	smat Trees;
	smat Parents;
	vec switch_acc;
	vec logpost;
	vec move_times;
	mat probs;
	svec testclass;

	outputs(): move_times(5) {};
	outputs(unsigned n_var, unsigned n_rows, unsigned n_step): 
		Groups(n_rows, n_var), Trees(n_rows, n_var), Parents(n_rows, n_var), switch_acc(n_step),
		logpost(n_rows), move_times(5) {};
};

///// Functions for computing the log posterior

field<svec> Categories(const smat &X, const unsigned n_var) {
	// returns the possible values for each variable (in increasing order)
	field<svec> cat(n_var);
	for(unsigned i=0; i < n_var; i++) {
		cat(i) = unique(X.col(i));
	}
	return cat;
}

uvec Level(const svec &vals, const svec &cat) {
	// assumes the categories (possible values) for each variable are nonnegative integers in increasing order
	uvec indices(vals.n_elem, fill::zeros);
	for(uword i=0; i < cat.n_elem - 1; i++)
		indices += (vals > cat(i));
	assert (all(vals == cat(indices)));
	return indices;
}

nlevels ComputeLevels(const data &Data, const unsigned n_var) {
	nlevels n_levels;
	svec x_nlevels(n_var);
	svec levels;
	for (unsigned i = 0; i < n_var; i++) {
		levels = unique(Data.X.col(i));
		x_nlevels(i) = levels.n_elem;
	}
	n_levels.x = x_nlevels;
	levels = unique(Data.Y);
	n_levels.y = levels.n_elem;
	return n_levels;
}

counts ComputeCounts(const data &Data, const parameters &Parameters, const nlevels &n_levels,
	const field<svec> &cat, const svec &cat_y) {
	const smat &X = Data.X_train;
	const svec &Y = Data.Y_train;
	unsigned n_units = Parameters.n_units;
	unsigned n_var = Parameters.n_var;
	counts Counts(n_var);

	svec count_y(n_levels.y);
	count_y.fill(0);
	uvec y_level = Level(Y, cat_y);
	for (unsigned i = 0; i < n_units; i++) {
		count_y(y_level(i))++;
	}
	Counts.y = count_y;
	assert(accu(Counts.y) == n_units); 

	for (unsigned var1 = 0; var1 < n_var; var1++) {
		uvec x1_level = Level(X.col(var1), cat(var1));

		svec count_var(n_levels.x(var1));
		smat count_var_y(n_levels.x(var1), n_levels.y);
		count_var.fill(0);
		count_var_y.fill(0);
		for (unsigned i = 0; i < n_units; i++) {
			count_var(x1_level(i))++;
			count_var_y(x1_level(i), y_level(i))++;
		}
		Counts.var(var1) = count_var;
		Counts.var_y(var1) = count_var_y;
		assert(accu(Counts.var(var1)) == n_units);
		assert(accu(Counts.var_y(var1)) == n_units);

		for (unsigned var2 = var1+1; var2 < n_var; var2++) {
			uvec x2_level = Level(X.col(var2), cat(var2));

			smat count_var_var(n_levels.x(var1), n_levels.x(var2));          
			scube count_var_var_y(n_levels.x(var1), n_levels.x(var2), n_levels.y);
			count_var_var.fill(0);           
			count_var_var_y.fill(0);
			for (unsigned i = 0; i < n_units; i++) {
				count_var_var(x1_level(i), x2_level(i))++;
				count_var_var_y(x1_level(i), x2_level(i), y_level(i))++;
			}
			Counts.var_var(var1, var2) = count_var_var;
			Counts.var_var_y(var1, var2) = count_var_var_y;
			assert(accu(Counts.var_var(var1, var2)) == n_units);
			assert(accu(Counts.var_var_y(var1, var2)) == n_units);

		}
	}
	return Counts;
}

cube LogLikTermMatrix(const counts &Counts, const nlevels &n_levels, const parameters &Parameters) {
	unsigned n_var = Parameters.n_var;
	cube term_matrix(n_var, n_var, 2);
	term_matrix.zeros();
	unsigned m = n_levels.y;
	double alpha = Parameters.alpha;
	for (unsigned var1 = 0; var1 < n_var; var1++) {

		unsigned v1 = n_levels.x(var1);
		for (unsigned i = 0; i < v1; i++) {
			term_matrix(var1, var1, noise) += lgamma(Counts.var(var1)(i) + alpha/v1) - lgamma(alpha/v1);
			for (unsigned l = 0; l < m; l++) {
				term_matrix(var1, var1, sig) +=
				lgamma(Counts.var_y(var1)(i, l) + alpha/(v1*m)) - lgamma(alpha/(v1*m));
			}
		}

		for (unsigned var2 = var1+1; var2 < n_var; var2++) {
			unsigned v2 = n_levels.x(var2);
			for (unsigned i = 0; i < v1; i++) {
				for (unsigned j = 0; j < v2; j++)	{
					term_matrix(var1, var2, noise) += lgamma(Counts.var_var(var1, var2)(i, j) +
						                              alpha/(v1*v2)) - lgamma(alpha/(v1*v2));
					for (unsigned l = 0; l < m; l++) {
						term_matrix(var1, var2, sig) += lgamma(Counts.var_var_y(var1, var2)(i, j, l) +
							                            alpha/(v1*v2*m)) - lgamma(alpha/(v1*v2*m));
					}
				}
			}
			term_matrix(var2, var1, noise) = term_matrix(var1, var2, noise);
			term_matrix(var2, var1, sig) = term_matrix(var1, var2, sig);
		}

	}
	return term_matrix;
}

vec LogLikTerms(const counts &Counts, const nlevels &n_levels, const parameters &Parameters) {
	// term 0 is over all data (corresponds to Group 0 roots)
	// term 1 is broken down by Y (corresponds to Group 1 roots)
	unsigned m = n_levels.y;
	double alpha = Parameters.alpha;
	vec terms(2);
	terms.zeros();
	terms(noise) = lgamma(Parameters.n_units + alpha) - lgamma(alpha);
	for (unsigned l = 0; l < m; l++) {
		terms(sig) += lgamma(Counts.y(l) + alpha/m) - lgamma(alpha/m);
	}
	return terms;
}

cube LogLik(const parameters &Parameters, const counts &Counts, const nlevels &n_levels) {
	unsigned n_var = Parameters.n_var;

	vec terms = LogLikTerms(Counts, n_levels, Parameters);
	cube term_matrix = LogLikTermMatrix(Counts, n_levels, Parameters);
	cube loglik_matrix(n_var, n_var, 2); 

	for (unsigned var1 = 0; var1 < n_var; var1++) {
		for (unsigned var2 = 0; var2 < n_var; var2++) {
			if (var1 == var2) {
				loglik_matrix(var1, var1, noise) = term_matrix(var1, var1, noise) - terms(noise);
				loglik_matrix(var1, var1, sig) = term_matrix(var1, var1, sig) - terms(sig);
			} else {				
				loglik_matrix(var1, var2, noise) = term_matrix(var1, var2, noise) - term_matrix(var2, var2, noise);
				loglik_matrix(var1, var2, sig) = term_matrix(var1, var2, sig) - term_matrix(var2, var2, sig);
			}
		}
	}
	for (unsigned var1 = 0; var1 < n_var; var1++) {
		for (unsigned var2 = 0; var2 < n_var; var2++) {
			for (unsigned group = 0; group <= 1; group++) {
				assert(abs(loglik_matrix(var1, var2, group) - loglik_matrix(var2, var1, group) -
				loglik_matrix(var1, var1, group) + loglik_matrix(var2, var2, group)) < cutoff_equal);
			}
		}
	}
	return loglik_matrix;
}

void LogPost(cube &logpost_matrix, const parameters &Parameters, const nlevels &n_levels) {
	unsigned n_var = Parameters.n_var;
	for (unsigned var1 = 0; var1 < n_var; var1++) {
		for (unsigned var2 = 0; var2 < n_var; var2++) {
			if (var1 != var2) {
				logpost_matrix(var1, var2, sig) -= (Parameters.edge_mult + Parameters.yedge_mult) *
				Parameters.scaling / n_levels.y;
				logpost_matrix(var1, var2, noise) -= Parameters.edge_mult * Parameters.scaling;
			} else {
				logpost_matrix(var1, var2, sig) -= Parameters.yedge_mult * Parameters.scaling / n_levels.y;
			}
		}
	}
}

///// Helper functions for performing the MCMC updates

double logsumexp(const vec &logpost) {
	double mlogpost = max(logpost);
	return log(sum(exp(logpost - mlogpost))) + mlogpost;
}

double RandUnif() {
	return runif(1)[0];
}

unsigned RandSample(const unsigned n_max) {
	assert(n_max > 0);
	return floor(runif(1, 0, n_max)[0]);
}

uvec RandShuffle(unsigned size) {
	IntegerVector seq = seq_len(size) - 1;
	IntegerVector shuffle = Rcpp::RcppArmadillo::sample(seq, size, false, NumericVector::create());
	return as<uvec>(shuffle); 
}

uvec SampleWithoutReplacement(const svec &orig, const unsigned k) {
	IntegerVector sample = Rcpp::RcppArmadillo::sample(as<IntegerVector>(wrap(orig)), k, false, NumericVector::create());
	return as<uvec>(sample); 
}

unsigned opp(const unsigned g) {
	assert(g == noise || g == sig);
	if (g == noise) return sig;
	return noise;
}

vec normalize(const vec &logpost) {
	double mlogpost = max(logpost);
	return (logpost - mlogpost) - logsumexp(logpost - mlogpost);
}

unsigned Choose(const vec &logpost) {
	assert(logpost.n_elem > 0);
	vec Ratio = exp(normalize(logpost));
	double unif = RandUnif();
	double s = 0;
	unsigned index = 0;
	for(unsigned i = 0; i < Ratio.n_elem; i++)
	{
		s = s + Ratio(i);
		if(s >= unif)
		{
			index = i;
			break;
		}
	}
	return index;
}

void FindRootNode(const svec &Parent, const unsigned &node) {
	unsigned pos = node;
	int steps = 0;

	while(pos != null_value) {
		unsigned next = Parent(pos);
		if(next == null_value) break;
		pos = next;
		steps++;
		assert(steps<1000); // Super-deep tree means we actually made a cycle
	}
}

unsigned FindRootTree(const graph &Graph, const unsigned &tree_label) {
	uvec tree_index_set = find(Graph.Tree == tree_label);
	uvec root_index_set = find(Graph.Parent.elem(tree_index_set) == null_value);
	assert(root_index_set.n_elem == 1);
	unsigned root = tree_index_set(root_index_set(0));
	assert(Graph.Tree(root) == tree_label);
	assert(Graph.Parent(root) == null_value);
	return root;
}

double LogPostDiffTree(const graph &Graph, const cube &logpost_matrix, const unsigned &tree_label) {
	uvec tree_index_set = find(Graph.Tree == tree_label);
	unsigned tree_group = Graph.Group(tree_index_set(0));
	double diff = 0;

	for (unsigned i = 0; i < tree_index_set.n_elem; i++) {
		unsigned j = tree_index_set(i);
		unsigned par = Graph.Parent(j);
		if (par == null_value) par = j;
		diff += logpost_matrix(j, par, opp(tree_group)) - 
			logpost_matrix(j, par, tree_group);
	}
	return diff;
}

double LogPostTree(const graph &Graph, const cube &logpost_matrix, const unsigned &tree_label, const unsigned &group) {
	uvec tree_index_set = find(Graph.Tree == tree_label);
	double logpost = 0;

	for (unsigned i = 0; i < tree_index_set.n_elem; i++) {
		unsigned j = tree_index_set(i);
		unsigned par = Graph.Parent(j);
		if (par == null_value) par = j;
		logpost += logpost_matrix(j, par, group);
	}
	return logpost;
}

void MergeTreeLabels(svec &Tree, unsigned chosen_tree_label, unsigned top_tree_label) {
	unsigned max_tree_label_old = max(Tree);
	assert(top_tree_label != chosen_tree_label);
	Tree.elem(find(Tree == chosen_tree_label)).fill(top_tree_label);
	uvec max_tree = find(Tree == max_tree_label_old);
	Tree.elem(max_tree).fill(chosen_tree_label);
	assert(max(Tree) == max_tree_label_old - 1);
}

void SplitTreeLabels(svec &Tree, const svec &Parent, unsigned node, unsigned tree_label) {
	uvec nodes(1);
	nodes(0) = node;
	uvec children;
	while (nodes.n_elem != 0) {
		Tree.elem(nodes).fill(tree_label);
		children.reset();
		for (unsigned i = 0; i < nodes.n_elem; i++) {
			children.insert_rows(children.n_elem, find(Parent == nodes(i)));
		}
		nodes = children;
	}
}

void SplitSubtree(graph &Graph, const unsigned &chosen_node) {
	svec &Tree = Graph.Tree;
	svec &Parent = Graph.Parent;
	if (Parent(chosen_node) != null_value) {
		SplitTreeLabels(Tree, Parent, chosen_node, max(Tree) + 1);
		Parent(chosen_node) = null_value;
	}
}

void MergeSubtree(graph &Graph, const unsigned &chosen_node, const unsigned &parent) {
	svec &Tree = Graph.Tree;
	Graph.Parent(chosen_node) = parent;
	if (parent != null_value) {
		MergeTreeLabels(Tree, Tree(chosen_node), Tree(parent));
	}
}

double LogPostProb(const graph &Graph, const cube &logpost_matrix, const parameters &Parameters) {
	double logpost = 0;
	for(unsigned j = 0; j < Graph.Group.n_elem; j++) {
		unsigned par = Graph.Parent(j);
		if(par == null_value) par = j;
		logpost += logpost_matrix(j, par, Graph.Group(j));
	}
	return logpost;
}

///// Debugging functions

void SanityCheck(const graph &Graph) {
	#ifdef DEBUG
		const svec &Group = Graph.Group;
		const svec &Tree = Graph.Tree;
		const svec &Parent = Graph.Parent;
		assert(Group.n_elem == Tree.n_elem);
		assert(Parent.n_elem == Tree.n_elem);
		assert(max(Group) <= 1);

		for (unsigned i = 0; i < Tree.n_elem; i++) {
			if (Parent(i) != null_value) {
				assert(Tree(i) == Tree(Parent(i)));
				assert(Group(i) == Group(Parent(i)));
			} 
		}

		for(unsigned i = 0; i < Parent.n_elem; i++)
			FindRootNode(Parent, i);

		svec treelabels = unique(Tree);
		for(unsigned i = 0; i <= max(Tree); i++) 
		{
			assert(treelabels(i) == i);
			FindRootTree(Graph, treelabels(i));
		}
	#endif
}

void DetailedBalanceCheck(const double &logpost1, const double &transprob1,
	const double &logpost2, const double &transprob2) {
	assert(abs(logpost1 + transprob1 - logpost2 - transprob2) < cutoff_equal);
}

///// MCMC updates

unsigned Switch(graph &Graph, const unsigned tree_label, const cube &logpost_matrix, const parameters &Parameters) {
	// switches a tree to the opposite group
	// returns whether the move was accepted
	uvec chosen_tree = find(Graph.Tree == tree_label);
	unsigned tree_size = chosen_tree.n_elem;
	assert(tree_size > 0);
	unsigned tree_group = Graph.Group(chosen_tree(0));

	double log_accept1 = min(0.0, LogPostDiffTree(Graph, logpost_matrix, tree_label));
	double log_unif = log(RandUnif());

	if(log_unif < log_accept1) {
		#ifdef DEBUG
			double logpost1 = LogPostProb(Graph, logpost_matrix, Parameters);
		#endif

		Graph.Group.elem(chosen_tree).fill(opp(tree_group));

		#ifdef DEBUG
			double logpost2 = LogPostProb(Graph, logpost_matrix, Parameters);
			double log_accept2 = min(0.0, LogPostDiffTree(Graph, logpost_matrix, tree_label));
			DetailedBalanceCheck(logpost1, log_accept1, logpost2, log_accept2);
		#endif
			return 1;
	}
	SanityCheck(Graph);
	return 0;
}

unsigned SwitchRepeat(graph &Graph, const cube &logpost_matrix, const parameters &Parameters) {
	// proposes to switch each of a set of trees to the opposite group
	unsigned num_trees = max(Graph.Tree)+1;
	uvec tree_index_set = SampleWithoutReplacement(Graph.Tree, min(Parameters.k, num_trees));
	unsigned count = 0;
	for(unsigned i=0; i<tree_index_set.n_elem; i++) {
		count += Switch(Graph, tree_index_set(i), logpost_matrix, Parameters);
	}
	SanityCheck(Graph);
	return count;
}

void Pivot(graph &Graph, const cube &logpost_matrix, const unsigned tree_label) {
	// pivot a tree to a randomly chosen new root
	svec &Parent = Graph.Parent;
	svec newParent = Parent;
	uvec tree_index_set = find(Graph.Tree == tree_label);
	unsigned node = tree_index_set(RandSample(tree_index_set.n_elem));
	if (Parent(node) == null_value) return;

	newParent(node) = null_value;
	unsigned new_parent = node;
	unsigned new_child = Parent(new_parent);

	while(new_child != null_value) {
		newParent(new_child) = new_parent;
		new_parent = new_child;
		new_child = Parent(new_parent);
	}

	Graph.Parent = newParent;
	SanityCheck(Graph);
}

void ReassignSubtreeChoose(const graph &Graph, const unsigned &chosen_node, const uvec &chosen_subtree, 
	svec &parent_subset, vec &logpost, const cube &logpost_matrix, const parameters &Parameters) {
	// chooses a node to attach the bottom tree to according to the conditional distribution
	unsigned n_var = Parameters.n_var;
	const svec &Tree = Graph.Tree;
	const svec &Group = Graph.Group;
	assert(Tree(chosen_node) == Tree(chosen_subtree(0)));
	unsigned chosen_tree_label = Tree(chosen_node);
	unsigned chosen_tree_size = chosen_subtree.n_elem;

	vec logpost_tree(2);
	logpost_tree(0) = LogPostTree(Graph, logpost_matrix, chosen_tree_label, 0); // chosen tree is in group 0
	logpost_tree(1) = LogPostTree(Graph, logpost_matrix, chosen_tree_label, 1); // chosen tree is in group 1

	unsigned subset_size = n_var - chosen_tree_size;

	parent_subset.reset();
	parent_subset.set_size(subset_size + 2);
	logpost.reset();
	logpost.set_size(subset_size + 2);

	unsigned count = 0;
	if (subset_size > 0) {
		for(unsigned par = 0; par < n_var; par++) {
			if (Tree(par) != chosen_tree_label) {
				parent_subset(count) = par;
				logpost(count) = logpost_matrix(chosen_node, par, Group(par)) - 
					logpost_matrix(chosen_node, chosen_node, Group(par)) + logpost_tree(Group(par));
				count++;
			}
		}
	}
	assert(count == subset_size);

	for(unsigned group=0; group<2; group++) {
		logpost(count) = logpost_tree(group);
		parent_subset(count) = null_value;
		count++;
	}
}

void ReassignSubtree(graph &Graph, const cube &logpost_matrix, const parameters &Parameters) {
	// reassigns a subtree to a different parent (can be null)
	svec &Group = Graph.Group;
	svec &Tree = Graph.Tree;

	unsigned chosen_node = RandSample(Parameters.n_var);
	Pivot(Graph, logpost_matrix, Tree(chosen_node));
	unsigned group1 = Group(chosen_node);

	#ifdef DEBUG
		unsigned par1 = Graph.Parent(chosen_node);
		double logpost1 = LogPostProb(Graph, logpost_matrix, Parameters);
	#endif

	SplitSubtree(Graph, chosen_node);
	uvec chosen_subtree = find(Tree == Tree(chosen_node));

	svec parent_subset;
	vec logpost;
	ReassignSubtreeChoose(Graph, chosen_node, chosen_subtree, parent_subset, logpost, logpost_matrix, Parameters);

	unsigned index = Choose(logpost);

	#ifdef DEBUG
		double tp1 = normalize(logpost)(index);
	#endif

	unsigned par2 = parent_subset(index);
	unsigned group2;
	if (par2 == null_value) group2 = index - logpost.n_elem + 2;
	else group2 = Group(par2);

	MergeSubtree(Graph, chosen_node, par2);
	if (group2 != group1) Group.elem(chosen_subtree).fill(group2);

	#ifdef DEBUG
		double logpost2 = LogPostProb(Graph, logpost_matrix, Parameters);
		SanityCheck(Graph);

		graph Graph2 = Graph;
		SplitSubtree(Graph2, chosen_node);
		ReassignSubtreeChoose(Graph2, chosen_node, chosen_subtree, parent_subset, logpost, logpost_matrix, Parameters);

		uvec par_index_set = find(parent_subset == par1);
		assert(par_index_set.n_elem > 0);
		unsigned par_index = 0;
		if (par1 == null_value) par_index = group1;
		double tp2 = normalize(logpost)(par_index_set(par_index));

		DetailedBalanceCheck(logpost1, tp1, logpost2, tp2);
	#endif

	SanityCheck(Graph);
}

///// Constructing starting graphs for MCMC

graph InitGraph(const parameters &Parameters) {
	unsigned n_var = Parameters.n_var;
	graph Graph(n_var);

	if (Parameters.start == "random_trees") {
		svec nodes = shuffle(linspace<svec>(0, n_var - 1, n_var));
		svec shuffle_groups(n_var);

		bool is_root;
		unsigned node, group;
		unsigned tree_count = 0, group_size = 0;
		uvec group_index_set;
		for (unsigned i = 0; i < n_var; i++) {
			node = nodes(i);
			group = (unsigned)(RandUnif() < 0.5);
			Graph.Group(node) = group;
			shuffle_groups(i) = group;
			if (i > 0) {
				group_index_set = find(shuffle_groups.subvec(0, i-1) == group);
				group_size = group_index_set.n_elem;
			}
			is_root = (group_size == 0) || (-log(RandUnif()) < Parameters.edge_mult * Parameters.scaling);

			if (is_root) {
				Graph.Parent(node) = null_value;
				Graph.Tree(node) = tree_count;
				tree_count++;
			} else {
				unsigned parent = nodes(group_index_set(RandSample(group_size)));
				Graph.Parent(node) = parent;
				Graph.Tree(node) = Graph.Tree(parent);
				assert(Graph.Group(parent) == group);
			}
		}

	} else {
		Graph.Tree = linspace<svec>(0, n_var - 1, n_var);
		Graph.Parent.fill(null_value);
		Graph.Group.zeros();
		if (Parameters.start == "random_groups") {
			uvec g1_indices = SampleWithoutReplacement(Graph.Tree, floor(n_var*0.5));
			Graph.Group.elem(g1_indices).ones();
		}
	}

	SanityCheck(Graph);
	return Graph;
}

graph TrueModelGraph(const imat &true_model, const unsigned n_var) {
	unsigned n_kernel = true_model.n_cols;
	svec trueGroup = to_svec(true_model.row(0).t());
	svec trueTree = to_svec(true_model.row(1).t());
	unsigned max_tree = max(trueTree);
	ivec Parent = true_model.row(2).t();
	Parent.elem(find(Parent == -1)).fill(null_value);
	svec trueParent = to_svec(Parent);

	graph trueGraph(n_var);
	trueGraph.Group.fill(0); 
	trueGraph.Parent.fill(null_value);
	trueGraph.Group.rows(0, n_kernel-1) = trueGroup;
	trueGraph.Parent.rows(0, n_kernel-1) = trueParent;
	trueGraph.Tree.rows(0, n_kernel-1) = trueTree;
	if (n_kernel < n_var) trueGraph.Tree.rows(n_kernel, n_var-1) =
		linspace<svec>(max_tree+1, n_var-n_kernel+max_tree, n_var-n_kernel);
	SanityCheck(trueGraph);
	return trueGraph;
}

///// Running MCMC

void MCMC(field<graph> &Graphs, vec &logpost,
	const data &Data, const cube &logpost_matrix, const parameters &Parameters, outputs &Outputs) {
	timeb t0, t1, t2, t3;
	vec move_times(5);
	move_times.zeros();
	//string output_id = Parameters.output_id;

	//unsigned n_step = Parameters.n_step;
	unsigned n_burnin = Parameters.n_step / Parameters.burnin_frac;
	
	assert(Graphs.n_elem == Parameters.n_samples);
	assert(logpost.n_elem == Parameters.n_samples);
	
	//unsigned n_var = Parameters.n_var;
	/*
	smat Groups(n_rows, n_var);
	smat Trees(n_rows, n_var);
	smat Parents(n_rows, n_var);

	vec switch_acc(n_step);
	vec logpost_all(n_step);*/

	graph Graph = InitGraph(Parameters);
	if (!Data.true_model.is_empty()) {
		graph TrueGraph = TrueModelGraph(Data.true_model, Parameters.n_var);
		vec true_logpost(1);
		true_logpost(0) = LogPostProb(TrueGraph, logpost_matrix, Parameters);
		//true_logpost.save(output_id + "_TrueModelLogPost.txt", csv_ascii);

		if (Parameters.start == "true") Graph = TrueGraph;
	}

	unsigned count = 0, count1 = 0;
	for(unsigned i = 0; i < Parameters.n_step; i++) {
		ftime(&t0);
		Outputs.switch_acc(i) = SwitchRepeat(Graph, logpost_matrix, Parameters);
		ftime(&t1);
		ReassignSubtree(Graph, logpost_matrix, Parameters);
		ftime(&t2);
		//Outputs.logpost_all(i) = LogPostProb(Graph, logpost_matrix, Parameters);

		if (i % Parameters.thin == 0) {
			if (Parameters.thin_output) {
				Outputs.Groups.row(count1) = Graph.Group.t();
				Outputs.Trees.row(count1) = Graph.Tree.t();
				Outputs.Parents.row(count1) = Graph.Parent.t();
				Outputs.logpost(count1) = LogPostProb(Graph, logpost_matrix, Parameters);
				count1++;
			}
			if (i >= n_burnin) {
				Graphs(count) = Graph;
				//logpost(count) = Outputs.logpost_all(i);
				logpost(count) = LogPostProb(Graph, logpost_matrix, Parameters);
				count++;
			}
		}
		if (!Parameters.thin_output) {
			Outputs.Groups.row(i) = Graph.Group.t();
			Outputs.Trees.row(i) = Graph.Tree.t();
			Outputs.Parents.row(i) = Graph.Parent.t();
		}
		ftime(&t3);

		move_times(0) += t1.time - t0.time + .001*(t1.millitm - t0.millitm);
		move_times(1) += t2.time - t1.time + .001*(t2.millitm - t1.millitm);
		move_times(2) += t3.time - t2.time + .001*(t3.millitm - t2.millitm);
	}

	assert(count == Parameters.n_samples);
	string suff = (Parameters.thin_output)?"":"_all";
	Outputs.Parents += 1;
	//Outputs.Parents.elem(find(Outputs.Parents == null_value + 1)).fill(0);
	/*
	imat printParents = conv_to<imat>::from(Parents);
	printParents.elem(find(Parents == null_value)).fill(-1);
	Groups.save(output_id + "_Groups" + suff + ".txt", raw_ascii);
	Trees.save(output_id + "_Trees" + suff + ".txt", raw_ascii);
	printParents.save(output_id + "_Parents" + suff + ".txt", raw_ascii);
	logpost_all.save(output_id + "_LogPost_all.txt", raw_ascii);
	logpost.save(output_id + "_LogPost.txt", raw_ascii);
	switch_acc.save(output_id + "_SwitchAcc_all.txt", raw_ascii);
	move_times.save(output_id + "_Move_Runtimes.txt", csv_ascii);*/
}

///// Classification functions

vec LogProbY(const graph &Graph, const counts &Counts, const data &Data, const nlevels &n_levels,
	const field<svec> &cat, const unsigned &test_row, const parameters &Parameters) {
	uvec group1 = find(Graph.Group == sig);
	unsigned n_units = Parameters.n_units;

	assert(n_units > 0);
	assert(min(Counts.y) > 0);

	vec logpost = log(to_vec(Counts.y)) - log(n_units);
	assert(logpost.n_elem == n_levels.y);

	for (unsigned y_level = 0; y_level < n_levels.y; y_level++) {
		for (unsigned j = 0; j < group1.n_elem; j++) {
			unsigned var = group1(j);
			unsigned var_level = Level(Data.X_test(span(test_row, test_row), var), cat(var))(0);
			double numer, denom;
			if (Graph.Parent(var) == null_value) {
				numer = Counts.var_y(var)(var_level, y_level) + Counts.var(var)(var_level) * 1.0 / n_units;
				denom = 1.0 + Counts.y(y_level);
			} else {
				unsigned par = Graph.Parent(var);
				unsigned par_level = Level(Data.X_test(span(test_row, test_row), par), cat(par))(0);
				unsigned count_edge;
				if (var < par) count_edge = Counts.var_var_y(var, par)(var_level, par_level, y_level);
				else count_edge = Counts.var_var_y(par, var)(par_level, var_level, y_level);
				numer = count_edge + Counts.var(var)(var_level) * 1.0 / n_units;
				denom = 1.0 + Counts.var_y(par)(par_level, y_level);
				assert (denom > 0); 
			}
			if (numer > 0) logpost(y_level) += log(numer) - log(denom);
		}
	}
	return logpost - logsumexp(logpost);
}

void Classify(const field<graph> &Graphs, const counts &Counts, const nlevels &n_levels, const field<svec> &cat, 
	const svec &cat_y, const vec &logpost, const data &Data, const parameters &Parameters, outputs &Outputs) {
	assert (Graphs.n_elem == logpost.n_elem);
	unsigned n_test = Data.X_test.n_rows;
	mat probs(n_test, n_levels.y);
	probs.fill(0); 
	double mlogpost = max(logpost);

	for (unsigned i = 0; i < Graphs.n_elem; i++) {
		for (unsigned j = 0; j < n_test; j++) {
			probs.row(j) += exp(LogProbY(Graphs(i), Counts, Data, n_levels, cat, j, Parameters) +
				logpost(i) - mlogpost).t();
		}
	}
	svec testclass(n_test);
	uvec max_index_set;
	for (unsigned j = 0; j < n_test; j++) {
		max_index_set = find(probs.row(j) == max(probs.row(j)));
		assert(max_index_set.n_elem > 0); // is usually violated in case of infinity errors somewhere
		testclass(j) = cat_y(max_index_set(0));
	}
	//probs.save(Parameters.output_id + "_Probabilities.txt", csv_ascii);
	Outputs.probs = probs;
	Outputs.testclass = testclass;
	//return testclass;
}

///// Running the SBFC algorithm (with or without cross-validation)
void SBFC(const data &Data, const parameters &Parameters, outputs &Outputs) {
	timeb t1, t2, t3, t4, t5;
	ftime(&t1);
	field<svec> cat = Categories(Data.X, Parameters.n_var);
	svec cat_y = unique(Data.Y);
	nlevels n_levels = ComputeLevels(Data, Parameters.n_var);
	counts Counts = ComputeCounts(Data, Parameters, n_levels, cat, cat_y);
	ftime(&t2);
	cube logpost_matrix = LogLik(Parameters, Counts, n_levels);
	LogPost(logpost_matrix, Parameters, n_levels);
	ftime(&t3);
	field<graph> Graphs(Parameters.n_samples);
	vec logpost(Parameters.n_samples);
	MCMC(Graphs, logpost, Data, logpost_matrix, Parameters, Outputs);
	ftime(&t4);
	Classify(Graphs, Counts, n_levels, cat, cat_y, logpost, Data, Parameters, Outputs);
	ftime(&t5);
	vec times(5);
	times(0) = t2.time - t1.time;
	times(1) = t3.time - t2.time;
	times(2) = t4.time - t3.time;
	times(3) = t5.time - t4.time;
	times(4) = t5.time - t1.time;
	Outputs.move_times = times;
	//times.save(Parameters.output_id + "_Runtimes.txt", csv_ascii);
	//return Outputs.testclass;
}

struct cv_fold {
	uvec test_subset;
	data Data;
	parameters Parameters;
	outputs Outputs;
};

void *CV_SBFC_fold(void *_fold) {
	cv_fold &fold = *(cv_fold *) _fold;
	SBFC(fold.Data, fold.Parameters, fold.Outputs);
	return NULL;
}

double CV_SBFC(const data &Data, const parameters &Parameters, outputs &Outputs) {
	unsigned n_units = Parameters.n_units;
	unsigned n_folds = Parameters.n_folds;
	unsigned fold_size = n_units / n_folds;
	unsigned mod = n_units % n_folds;
	unsigned fold_start = 0, fold_end = 0;
	uvec row_shuffle = RandShuffle(n_units);
	svec testclass(n_units);

	vector<cv_fold> cv_folds(n_folds);
	/*
	vector<pthread_t> threads(n_folds);
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int return_code;
	*/

	for (unsigned i=0; i < n_folds; i++) {
		uvec train_subset;
		fold_end = fold_start + fold_size;
		if (i >= mod) fold_end--;
		cv_folds[i].test_subset = row_shuffle.rows(fold_start, fold_end);
		if (i == 0) {
			train_subset = row_shuffle.rows(fold_end+1, n_units-1);
		} else if (i+1 == n_folds) {
			train_subset = row_shuffle.rows(0, fold_start-1);
		} else {
			train_subset = join_cols(row_shuffle.rows(0, fold_start-1), row_shuffle.rows(fold_end+1, n_units-1));
		}
		assert(cv_folds[i].test_subset.n_elem + train_subset.n_elem == n_units);

		cv_folds[i].Data.X_train = Data.X_train.rows(train_subset);
		cv_folds[i].Data.X_test = Data.X_train.rows(cv_folds[i].test_subset);
		cv_folds[i].Data.Y_train = Data.Y_train.rows(train_subset);
		cv_folds[i].Data.Y_test = Data.Y_train.rows(cv_folds[i].test_subset);
		cv_folds[i].Data.X = Data.X_train;
		cv_folds[i].Data.Y = Data.Y_train;

		cv_folds[i].Parameters = Parameters;
		cv_folds[i].Parameters.n_units = train_subset.n_elem;
		cv_folds[i].Outputs = Outputs;
		
		ostringstream fold_id_stream;
		fold_id_stream << Parameters.output_id << "_" << (i+1);
		cv_folds[i].Parameters.output_id = fold_id_stream.str();

		/*
		if (Parameters.thread) {
			return_code = pthread_create(&threads[i], &attr, CV_SBFC_fold, (void *) &cv_folds[i]);
			if(return_code) {
				cerr << "pthread_create failed with return code " << return_code << endl;
				exit(-1);
			}
		} else {
		*/
		CV_SBFC_fold((void *) &cv_folds[i]);
		//}

		fold_start = fold_end + 1;
	}

	// pthread_attr_destroy(&attr);
	for (unsigned i=0; i < n_folds; i++) {
		/*if (Parameters.thread) {
			return_code = pthread_join(threads[i], NULL);
			if(return_code) {
				cerr << "pthread_join failed with return code " << return_code << endl;
				exit(-1);
			}
		}*/

		testclass.rows(cv_folds[i].test_subset) = cv_folds[i].Outputs.testclass;
	}

	assert(fold_start == n_units);
	Outputs = cv_folds[0].Outputs;
	Outputs.testclass = testclass;
	//testclass.save(Parameters.output_id + "_Predictions.txt", csv_ascii);
	uvec correct = find(testclass == Data.Y_train);
	double accuracy = correct.n_elem * 1.0 / Data.Y_train.n_elem;
	return accuracy;
}

double RunSBFC(const data &Data, parameters &Parameters, outputs &Outputs) {
	double accuracy = 0;
	if (Data.X_test.is_empty()) { // use k-fold cross-validation to compute accuracy
		accuracy = CV_SBFC(Data, Parameters, Outputs); 
	} else { // use test data set to compute accuracy
		SBFC(Data, Parameters, Outputs);
		//testclass.save(Parameters.output_id + "_Predictions.txt", csv_ascii);	
		if (!Data.Y_test.is_empty()) {
			uvec correct = find(Outputs.testclass == Data.Y_test);
			accuracy = correct.n_elem * 1.0 / Data.Y_test.n_elem;
		}
	} 
	return accuracy;
}

///// Initializing the SBFC algorithm

bool FileExists(string &path, string &file_name) {
	string fullpath = path + file_name;
	const char* file_path = fullpath.c_str();
	ifstream infile(file_path);
	return infile;
}

void DataLoad(data &Data, string &path, string (&filenames)[4]) {
	if (FileExists(path, filenames[0]) && FileExists(path, filenames[1])) {
		Data.X_train.load(path + filenames[0]);
		Data.Y_train.load(path + filenames[1]);
		assert(Data.X_train.n_rows == Data.Y_train.n_rows);
		cout << "training data loaded" << endl;
	} else {
		cout << "Please input the correct training data file names" << endl; 
		exit(0);
	}

	if (FileExists(path, filenames[2])) {
		Data.X_test.load(path + filenames[2]);
		cout << "test data loaded" << endl;
		assert(Data.X_test.n_cols == Data.X_train.n_cols);
		Data.X = join_cols(Data.X_train, Data.X_test);
		if (FileExists(path, filenames[3])) {
			Data.Y_test.load(path + filenames[3]);
			assert(Data.X_test.n_rows == Data.Y_test.n_rows);
			Data.Y = join_cols(Data.Y_train, Data.Y_test);
		} else {
			Data.Y = Data.Y_train;
		}
	} else {
		Data.X = Data.X_train;
		Data.Y = Data.Y_train;
	}
}

void SetParam(parameters &Parameters, unsigned n_var, unsigned n_units) {
	Parameters.n_units = n_units;
	if (Parameters.n_var == 0) Parameters.n_var = n_var;
	assert(Parameters.n_var <= n_var);
	if (Parameters.n_var >=1000) Parameters.thin_output = true;
	if (Parameters.n_step == 0) Parameters.n_step = max((unsigned)10000, 10 * Parameters.n_var);
	Parameters.n_rows = Parameters.thin_output?(Parameters.n_step/Parameters.thin):Parameters.n_step;
	Parameters.n_samples = (Parameters.n_step - Parameters.n_step/Parameters.burnin_frac) / Parameters.thin;
	Parameters.scaling = log(Parameters.n_var);
}

void InitParam(int argc, char* argv[], data &Data, parameters &Parameters, string &output_id, 
	ofstream &out_stream, bool &multichain) {
	string filenames[4] = {"TrainX", "TrainY", "TestX", "TestY"};
	string output_dir_suff, root_dir, id;
	bool true_model = false;

	// initializing parameters from user input
	for (int i = 1; i < argc; i++) {
		if(strcmp(argv[i], "-TrainX") == 0) {
			filenames[0] = string(argv[i+1]);
			i++;
		} else if(strcmp(argv[i], "-TrainY") == 0) {
			filenames[1] = string(argv[i+1]);
			i++;
		} else if(strcmp(argv[i], "-TestX") == 0) {
			filenames[2] = string(argv[i+1]);
			i++;
		} else if(strcmp(argv[i], "-TestY") == 0) {
			filenames[3] = string(argv[i+1]);
			i++;            
		} else if(strcmp(argv[i], "-dir") == 0) { 
			root_dir = string(argv[i+1]);
			i++;
		} else if(strcmp(argv[i], "-outdirsuff") == 0) { 
			// custom suffix for the name of output directory (will be created within the above directory)
			output_dir_suff = string(argv[i+1]);
			i++;
		} else if(strcmp(argv[i], "-id") == 0) { // data set id
			id = string(argv[i+1]);
			i++;
		} else if(strcmp(argv[i], "-alpha") == 0) { // Dirichlet hyperparameter value
			Parameters.alpha = atof(argv[i+1]);
			i++;   
		} else if(strcmp(argv[i], "-edge") == 0) { // multiplier for x-edge penalty
			Parameters.edge_mult = atof(argv[i+1]);
			i++;        
		} else if(strcmp(argv[i], "-yedge") == 0) { // multiplier for y-edge penalty
			Parameters.yedge_mult = atof(argv[i+1]);
			i++;               
		} else if(strcmp(argv[i], "-k") == 0) { // number of repeats for the Switch move
			Parameters.k = atof(argv[i+1]);
			i++;            
		} else if(strcmp(argv[i], "-nstep") == 0) { // total number of MCMC steps, including burnin
			Parameters.n_step = atof(argv[i+1]);
			i++;
		} else if(strcmp(argv[i], "-burninfrac") == 0) {
			// (denominator of) fraction of total MCMC steps discarded as burnin 
			Parameters.burnin_frac = atof(argv[i+1]);
			i++;
		} else if(strcmp(argv[i], "-thin") == 0) { // thinning factor for MCMC samples
			Parameters.thin = atof(argv[i+1]);
			i++;
		} else if(strcmp(argv[i], "-start") == 0) { 
			// starting graph for the algorithm (single noise nodes, random groups with single nodes, or random trees) 
			Parameters.start = string(argv[i+1]);
			i++;
		} else if(strcmp(argv[i], "-thinoutput") == 0) { // whether to thin the output files
			Parameters.thin_output = true;
		} else if(strcmp(argv[i], "-truemodel") == 0) { 
			// whether to compute the log posterior for true model (requires TrueModel file)
			true_model = true;
		} else if(strcmp(argv[i], "-multichain") == 0) { // whether to run multiple chains
			multichain = true;
		} else if(strcmp(argv[i], "-nvar") == 0) { // number of variables
			Parameters.n_var = atof(argv[i+1]);
			i++;        
		} else if(strcmp(argv[i], "-thread") == 0) { // whether to run cross-validation with parallel threads
			Parameters.thread = true;
		}
	}
	if (id.empty()) {
		cout << "Please input the data set id" << endl; 
		exit(0);
	}  
	if (root_dir.empty()) root_dir = "";
	string path = root_dir + id + "/";

	DataLoad(Data, path, filenames);
	if (true_model) Data.true_model.load(path + "TrueModel"); 
	SetParam(Parameters, Data.X_train.n_cols, Data.Y_train.n_elem);

	// setting target directory
	if (output_dir_suff.empty()) output_dir_suff = "";
	ostringstream outdir_stream;
	outdir_stream << "sbfc_nvar" << Parameters.n_var << '_' << Parameters.n_step/1000 << "k" << output_dir_suff;
	string output_dir = outdir_stream.str();
	output_id = path + output_dir + "/" + id + "_" + output_dir;
	string command = "mkdir \"" + path + output_dir + "\" ";
	system(command.c_str());

	// outputting parameter values
	string outfile = output_id + "_Output.txt";
	out_stream.open(outfile.c_str(), ofstream::app);

	out_stream  << output_dir << ": data set " << id << " with " << Parameters.n_var << " variables and " 
				<< Parameters.n_units << " units" << endl;
	cout        << output_dir << ": data set " << id << " with " << Parameters.n_var << " variables and " 
				<< Parameters.n_units << " units" << endl;
	out_stream  << Parameters.n_step << " total iterations, with 1/" << Parameters.burnin_frac
				<< " burnin, thinning every " << Parameters.thin << endl;
	cout        << Parameters.n_step << " total iterations, with 1/" << Parameters.burnin_frac
				<< " burnin, thinning every " << Parameters.thin << endl;
	out_stream  << "Edge penalty parameters: y-edge penalty multiplier=" << Parameters.yedge_mult 
				<< ", x-edge penalty multiplier=" << Parameters.edge_mult << ", scaling factor = "
				<< Parameters.scaling << endl;
	cout        << "Edge penalty parameters: y-edge penalty multiplier=" << Parameters.yedge_mult
				<< ", x-edge penalty multiplier=" << Parameters.edge_mult << ", scaling factor = "
				<< Parameters.scaling << endl;
	out_stream  << "Other parameters: alpha=" << Parameters.alpha << ", k=" << Parameters.k 
				<< ", start=" << Parameters.start << endl;
	cout        << "Other parameters: alpha=" << Parameters.alpha << ", k=" << Parameters.k 
				<< ", start=" << Parameters.start << endl;
}

#ifndef TEST
int main(int argc, char* argv[]) {
	data Data;
	parameters Parameters;
	outputs Outputs;
	string output_id;
	ofstream out_stream;
	bool multichain = false;

	InitParam(argc, argv, Data, Parameters, output_id, out_stream, multichain);

	timeb start, end;
	ftime(&start);
	vec accuracy(1);
	if (multichain) {
		for (int k = 0; k <= 4; k++) {
			ostringstream chain_id_stream;
			chain_id_stream << output_id << "_chain" << k+1;
			Parameters.output_id = chain_id_stream.str();
			accuracy.set_size(5);
			accuracy(k) = RunSBFC(Data, Parameters, Outputs);

			out_stream 	<< "Chain " << k+1 << " accuracy is " << accuracy(k) << endl;
			cout 		<< "Chain " << k+1 << " accuracy is " << accuracy(k) << endl;
		}
	} else {
		Parameters.output_id = output_id;
		accuracy(0) = RunSBFC(Data, Parameters, Outputs);

		out_stream 	<< "Accuracy: " << accuracy(0) << endl;
		cout 		<< "Accuracy: " << accuracy(0) << endl;
	}
	accuracy.save(output_id + "_Accuracy.txt", csv_ascii);

	ftime(&end);
	int t = floor(end.time - start.time);
	int hour, min, sec;
	hour = t/3600; 
	t = t%3600; 
	min = t/60; 
	t = t%60; 
	sec = t; 

	out_stream 	<< "Time: " << end.time - start.time << "s = (" << hour << ":" << min << ":" << sec << ")" << endl;
	cout 		<< "Time: " << end.time - start.time << "s = (" << hour << ":" << min << ":" << sec << ")" << endl;

	out_stream.close();
	return 0;
}
#endif

void DataImportR(data &Data, SEXP &TrainX, SEXP &TrainY, SEXP &TestX, SEXP &TestY) {
	if ((TrainX != R_NilValue) && (TrainY != R_NilValue)) {
		Data.X_train = as<smat>(TrainX);
		Data.Y_train = as<svec>(TrainY);
		assert(Data.X_train.n_rows == Data.Y_train.n_rows);
	} else {
		Rf_error("Training data missing");
	}
	if (TestX != R_NilValue) {
		Data.X_test = as<smat>(TestX);
		assert(Data.X_test.n_cols == Data.X_train.n_cols);
		Data.X = join_cols(Data.X_train, Data.X_test);
		if (TestY != R_NilValue) {
			Data.Y_test = as<svec>(TestY);
			assert(Data.X_test.n_rows == Data.Y_test.n_rows);
			Data.Y = join_cols(Data.Y_train, Data.Y_test);
		} else {
			Data.Y = Data.Y_train;
		}
	} else {
		Data.X = Data.X_train;
		Data.Y = Data.Y_train;
	}
}


//' @title
//' Selective Bayesian Forest Classifier algorithm
//' @description
//' Runs the SBFC algorithm on a discretized data set.
//' 
//' @param TrainX matrix containing the training data
//' @param TrainY vector containing the class labels for the training data
//' @param TestX matrix containing the test data, if applicable
//' @param TestY vector containing the class labels for the test data
//' @details
//' Data needs to be discretized before running SBFC.
//' If the test data matrix TestX is provided, SBFC runs on the entire training set TrainX, and provides predicted class labels for the test data. If the test data class vector TestY is provided, the accuracy is computed. If the test data matrix TestX is not provided, SBFC performs cross-validation on the training data set TrainX, and returns predicted classes and accuracy for the training data.  
//' The number of MCMC iterations is \code{max(10000, 10 * n_var)}. 
//' For data sets with 1000 or more variables, the output matrices are thinned by default, and contain only the thinned samples used for classification.
//' @return An object of class \code{sbfc}:
//' \describe{     
//' \item{\code{accuracy}}{classification accuracy (on the test set if provided, otherwise cross-validation accuracy on training set)}
//' \item{\code{predictions}}{vector of class label predictions (for the test set if provided, otherwise for the training set)}
//' \item{\code{probabilities}}{matrix of class label probabilities (for the test set if provided, otherwise for the training set)}
//' \item{\code{runtime}}{total runtime of the algorithm in seconds}
//' \item{\code{parents}}{matrix representing the structures sampled by MCMC, where parents[i,j] is the index of the parent of node j at iteration i (0 if node is a root)}
//' \item{\code{groups}}{matrix representing the structures sampled by MCMC, where groups[i,j] indicates which group node j belongs to at iteration j (0 is noise, 1 is signal)}
//' \item{\code{trees}}{matrix representing the structures sampled by MCMC, where trees[i,j] indicates which tree node j belongs to at iteration j}
//' \item{\code{logposterior}}{vector representing the log posterior at each iteration of the MCMC}
//' }
//' @examples
//' data(chess)
//' chess_result = sbfc(as.matrix(chess$TrainX), as.numeric(chess$TrainY), as.matrix(chess$TestX), as.numeric(chess$TestY))
//' data(corral)
//' corral_result = sbfc(as.matrix(corral$TrainX), as.numeric(corral$TrainY)) # uses cross-validation
//' @export
// [[Rcpp::export]]
List sbfc(SEXP TrainX = R_NilValue, SEXP TrainY = R_NilValue, SEXP TestX = R_NilValue, SEXP TestY = R_NilValue) {
  timeb start, end;
  ftime(&start);
	data Data;
	DataImportR(Data, TrainX, TrainY, TestX, TestY);
	parameters Parameters;
	SetParam(Parameters, Data.X_train.n_cols, Data.Y_train.n_elem);
	outputs Outputs(Parameters.n_var, Parameters.n_rows, Parameters.n_step);
	double accuracy = RunSBFC(Data, Parameters, Outputs);
	ftime(&end);
	List results = List::create(
	  _["accuracy"] = accuracy,
	  _["predictions"] = as<IntegerVector>(wrap(Outputs.testclass)),
	  _["probabilities"] = as<NumericMatrix>(wrap(Outputs.probs)),
	  _["runtime"] = floor(end.time - start.time),
	  _["parents"] = as<IntegerMatrix>(wrap(Outputs.Parents)),
	  _["groups"] = as<IntegerMatrix>(wrap(Outputs.Groups)),
	  _["trees"] = as<IntegerMatrix>(wrap(Outputs.Trees)),
	  _["logposterior"] = as<NumericVector>(wrap(Outputs.logpost))
	);
	Rf_setAttrib(results, wrap("class"), wrap("sbfc"));
	return results;
}
