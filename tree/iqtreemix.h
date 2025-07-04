//
//  iqtreemix.hpp
//  tree
//
//  Created by Thomas Wong on 14/12/20.
//

#ifndef iqtreemix_h
#define iqtreemix_h

#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <cmath>
#include "iqtree.h"
#include "utils/MPIHelper.h"
#include "model/modelmarkov.h"
#include "model/ratefree.h"
// #include "phylokernel.h"
//#include "phylokernelnew.h"
#include <vectorclass/vectorclass.h>

// by default, the minimum branch length for MAST model is 1e-6
#define MAST_MIN_BRANCH_LEN 1e-6;

// for checking the scaling for the likelihood values
#define TINY_SCALE_DIFF 0.5
#define ONE_LOG_SCALE_DIFF 178.0

// change to 0.04 for tree mixture model as 0.02 and 0.03 cause numerical problems
const double MIN_GAMMA_SHAPE_TREEMIX = 0.04;

class IQTreeMix : public IQTree, public vector<IQTree*> {
public:

    /**
     default constructor
     */
    IQTreeMix();
    
    IQTreeMix(Params &params, Alignment *aln);
    
    /**
     destructor
     */
    virtual ~IQTreeMix() override;

    /**
     initialization
     */
    // void init(Alignment* aln, Params* params, Checkpoint* chkpt);
    
    virtual void initializeModel(Params &params, string model_name, ModelsBlock *models_block) override;

    /**
     * @return number of elements per site lhl entry, used in conjunction with computePatternLhCat
     */
    virtual int getNumLhCat(SiteLoglType wsl) override;

    // compute the overall likelihood value by combining all the existing likelihood values of the trees
    double computeLikelihood_combine(double *pattern_lh = nullptr, bool save_log_value = true);

    // compute the log-likelihood values for every site and tree
    // updated array: _ptn_like_cat
    // update_which_tree: only that tree has been updated
    void computeSiteTreeLogLike(int update_which_tree);

    virtual double computeLikelihood(double *pattern_lh = nullptr, bool save_log_value = true) override;

    virtual double computePatternLhCat(SiteLoglType wsl) override;

    /**
     compute pattern likelihoods only if the accumulated scaling factor is non-zero.
     Otherwise, copy the pattern_lh attribute
     @param pattern_lh (OUT) pattern log-likelihoods,
     assuming pattern_lh has the size of the number of patterns
     @param cur_logl current log-likelihood (for sanity check)
     @param pattern_lh_cat (OUT) if not nullptr, store all pattern-likelihood per category
     */
    virtual void computePatternLikelihood(double *pattern_lh = nullptr, double *cur_logl = nullptr,
                                          double *pattern_lh_cat = nullptr, SiteLoglType wsl = WSL_RATECAT) override;

    virtual void initializeAllPartialLh() override;

    virtual void deleteAllPartialLh() override;

    virtual void clearAllPartialLH(bool make_null = false) override;

    /**
            compute pattern posterior probabilities per rate/mixture category
            @param pattern_prob_cat (OUT) all pattern-probabilities per category
            @param wsl either WSL_RATECAT, WSL_MIXTURE or WSL_MIXTURE_RATECAT
     */
    virtual void computePatternProbabilityCategory(double *pattern_prob_cat, SiteLoglType wsl) override;

    /**
     optimize all branch lengths of one tree
     @param my_iterations number of iterations to loop through all branches
     */
    void optimizeAllBranchesOneTree(int whichtree, int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100);
    
    /**
     optimize all branch lengths of all trees
     @param my_iterations number of iterations to loop through all branches
     @return the likelihood of the tree
     */
    virtual double optimizeAllBranches(int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100) override;

    /**
     compute the updated tree weights according to the likelihood values along each site
     prerequisite: computeLikelihood() has been invoked

     */
    double optimizeTreeWeightsByEM(double* pattern_mix_lh, double logl_epsilon, int max_steps, bool& tree_weight_converge);
    double optimizeTreeWeightsByBFGS(double gradient_epsilon);

    double optimizeBranchLensByBFGS(double gradient_epsilon);

    // save branch lengths of all trees
    // node and dad are always nullptr
    void getBranchLengths(vector<DoubleVector> &len, Node *node = nullptr, Node *dad = nullptr);

    // restore branch lengths of all trees
    // node and dad are always nullptr
    void setBranchLengths(vector<DoubleVector> &len, Node *node = nullptr, Node *dad = nullptr);
    
    // virtual void showTree();
    
    /** set the root by name
     @param my_root root node name
     @param multi_taxa TRUE if my_root is a comma-separated list of nodes
     */
    virtual void setRootNode(const char *my_root, bool multi_taxa = false) override;
    
    /**
     @return true if this is a mixture of trees, default: false
     */
    virtual bool isTreeMix() override { return true; }

    /**
     set checkpoint object
     @param checkpoint a checkpoint
     */
    virtual void setCheckpoint(Checkpoint *checkpoint) override;

    virtual void startCheckpoint() override;

    void saveCheckpoint() override;
    
    void restoreCheckpoint() override;
    
    void setMinBranchLen(Params& params);

    /** set pointer of params variable */
    virtual void setParams(Params* params) override;

    /*
     * Generate the branch IDs
     * Branches of different trees with the same partition share the same ID
     */
    void computeBranchID();

    /**
     * Generate the initial tree (usually used for model parameter estimation)
     */
    void computeInitialTree(LikelihoodKernel kernel, istream* in = nullptr) override;

    /**
     * setup all necessary parameters
     */
    virtual void initSettings(Params& params) override;

    /**
     * compute the memory size required for storing partial likelihood vectors
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequired(size_t ncategory = 1, bool full_mem = false) override;

    /**
     * compute the memory size for top partitions required for storing partial likelihood vectors
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequiredThreaded(size_t ncategory = 1, bool full_mem = false) override;

    virtual void setNumThreads(int num_threads) override;

    /**
     test the best number of threads
     */
    virtual int testNumThreads() override;

    /**
     Initialize the tree weights using parsimony scores
     Idea:
     1. Check the parsimony score for each tree along all the sites
     2. Select the sites with different parsimony score between the trees.
     3. For each selected site, we check which parsimony score of the tree is minimum, and assign the site to the tree.
     4. The tree weights are estimated according to the proportion of the sites assigned to each tree.
     */
    void initializeTreeWeights();
    // void initializeTreeWeights2();

    virtual string optimizeModelParameters(bool printInfo, double logl_epsilon) override;

    /**
     print tree to .treefile
     @param suffix suffix of the output file
     */
    virtual void printResultTree(string suffix = "") override;

    /**
     * Return the tree string contining taxon names and branch lengths
     * @return a tree string
     */
    virtual string getTreeString() override;

    /**
     @return the weighted sum of the tree lengths
     @param node the starting node, nullptr to start from the root
     @param dad dad of the node, used to direct the search
     */
    virtual double treeLength(Node *node = nullptr, Node *dad = nullptr) override;

    /**
     @return the weighted sum of the lengths of all internal branches
     @param node the starting node, nullptr to start from the root
     @param dad dad of the node, used to direct the search
     */
    virtual double treeLengthInternal(double epsilon, Node *node = nullptr, Node *dad = nullptr) override;

    virtual int getNParameters();

    virtual void drawTree(ostream &out, int brtype = WT_BR_SCALE + WT_INT_NODE, double zero_epsilon = 2e-6) override;

    /**
     print the tree to the output file in newick format
     @param out the output file.
     @param node the starting node, nullptr to start from the root
     @param dad dad of the node, used to direct the search
     @param brtype type of branch to print
     @return ID of the taxon with smallest ID
     */
    // virtual int printTree(ostream &out, int brtype, Node *node, Node *dad = nullptr);
    /**
            print the tree to the output file in newick format
            @param out the output stream.
            @param brtype type of branch to print
     */
    virtual void printTree(ostream & out, int brtype = WT_BR_LEN) override;

    /**
     *  Return best tree string from the candidate set
     *
     *  @param numTrees
     *      Number of best trees to return
     *  @return
     *      A string vector of trees
     */
    virtual vector<string> getBestTrees(int numTrees = 0) override;
    
    /**
     Read the tree saved with Taxon IDs and branch lengths.
     @param tree_string tree string to read from
     */
    virtual void readTreeString(const string &tree_string) override;

    virtual ModelFactory *getModelFactory() override {
        return at(0)->getModelFactory();
    }

    /**
     get rate heterogeneity
     @return associated rate heterogeneity class
     */
    virtual RateHeterogeneity *getRate() override {
        return at(0)->getRate();
    }

    virtual ModelSubst *getModel() override {
        return at(0)->getModel();
    }

    /**
     get the name of the model
     */
    virtual string getModelName() override;

    /**
     get the name of the model
     */
    virtual string getModelNameParams(bool show_fixed_params = false) override;

    // compute parsimony scores for each tree along the patterns
    // results are stored in the array patn_parsimony
    void computeParsimony();

    /**
     whether pattern is constant
     */
    int* patn_isconst;

    /**
     parsimony scores for each tree along the patterns
     */
    int* patn_parsimony;

    /**
     weights of trees
     */
    vector<double> weights;
    vector<double> tmp_weights; // for optimization
    vector<double> weight_logs; // logarithm of weights

    /**
     ratios of parsimony informative sites with max posterior probability for each tree
     */
    vector<double> max_posterior_ratio;
    
    /**
     ratios of parsimony informative sites with max high-enought likelihood for each tree
     */
    vector<double> max_like_ratio;

    /**
     pattern likelihoods for all trees
     */
    double* ptn_like_cat;

    /**
     scaling factors for all trees
     */
    double* ptn_scale_cat;

    /**
     log-likelihoods of each tree for each pattern
     */
    double* single_ptn_tree_like;

    /**
     log-likelihoods for each pattern
     */
    double* ptn_like;

    /**
     models
     */
    vector<ModelSubst*> models;
    
    /**
     site rates
     */
    vector<RateHeterogeneity*> site_rates;

    /**
     trees assigned for the site rates
     */
    vector<PhyloTree*> site_rate_trees;

    /**
     members of each weight group
     */
    vector<vector<int> > weight_group_member;
    
    /**
     does weight group exist
     */
    bool weightGrpExist;

    /**
     branch ID
     branches of different trees with the same partition are assigned to the same ID
     */
    vector<int> branch_id;
    
    /**
     branch group
     groups of branches of different trees with the same partition
     */
    vector<IntVector> branch_group;
    
    /**
     branch lengths (for optimization)
     */
    vector<DoubleVector> branch_len;
    
    /**
     parsimony scores are computed for each tree along the patterns
     */
    bool parsi_computed;

    /**
     whether the submodels are linked
     */
    bool isLinkModel;
    
    /**
     whether the site rates are linked, if exists
     */
    bool isLinkSiteRate;

    /**
     the estimated average branch length for each tree
     */
    vector<double> estAvgBrlen;

protected:

    // to separate the submodel names and the site rate names from the full model name
    void separateModel(string modelName);
    
    // reset the ptn_freq array to the original frequencies of the patterns
    void resetPtnOrigFreq();

    /**
     get posterior probabilities along each site for each tree
     */
    void getPostProb(double* pattern_mix_lh, bool need_computeLike = true, int updated_which_tree = -1, bool need_multiplyFreq = true);

    /**
     update the ptn_freq array according to the posterior probabilities along each site for each tree

     */
    void computeFreqArray(double* pattern_mix_lh, bool need_computeLike = true, int update_which_tree = -1);
    
    /**
     optimize tree k separately
     */
    void optimizeTreeSeparately(int k, bool printInfo, double logl_epsilon, double gradient_epsilon);
    
    /**
     optimize each tree separately
     */
    void optimizeTreesSeparately(bool printInfo, double logl_epsilon, double gradient_epsilon);

    /**
     If there are multiple tree weights belonging to the same group
     set all the tree weights of the same group to their average
     */
    void checkWeightGrp();

    /**
     If there are multiple branches belonging to the same group
     set all the branches of the same group to their average
     */
    void checkBranchGrp();

    // For the linked RHAS model
    // Store the RHAS variables of tree 0 to the array rhas_var
    void storeTree0RHAS();

    // For the linked RHAS model
    // Replace the RHAS variables of tree t by those of tree 0
    // The array rhas_var should have stored the updated RHAS variables of tree 0
    void copyRHASfrTree0(int t);

    // -------------------------------------
    // for BFGS optimzation on tree weights
    // -------------------------------------

    double targetFunk(double x[]) override;

    // read the tree weights and write into "variables"
    void setVariables(double *variables);

    // read the "variables" and write into tree weights
    void getVariables(double *variables);

    // set the bounds
    void setBounds(double *lower_bound, double *upper_bound, bool* bound_check);
    
    // get the dimension of the variables (for tree weights)
    int getNDim() override;
    
    // optimization on which variable
    // 1 - tree weights
    // 2 - branch lengths
    int optim_type;
    
    /**
     immediate array for pattern likelihoods during computation
     */
    double* _ptn_like_cat;

    /**
     number of optimization steps, default: number of Trees * 2
     */
    int optimize_steps;
    
    /**
     inputted tree-mixture model name
     */
    string treemix_model;
    
    /**
     individual model names (only one item if linked, more than one otherwise)
     */
    vector<string> model_names;
    
    /**
     individual site rate names (only one item if linked, more than one otherwise)
     */
    vector<string> siterate_names;
    
    /**
     whether there is any site rate
     */
    bool anySiteRate;

    /**
     whether it is a edge-len-restricted model in which the edges of different trees having the same partition have the same lengths
     */
    bool isEdgeLenRestrict;

    /**
     number of trees
     */
    size_t ntree;

    /**
     number of tips
     */
    size_t ntip;

    /**
     number of patterns
     */
    size_t nptn;
    
    /**
     number of parsimony informative sites
     */
    size_t ninformsite;
    
    /**
     number of branches of each tree
     */
    size_t nbranch;

    /**
     initial models
     */
    vector<ModelSubst*> initial_models;
    
    /**
     initial site rates
     */
    vector<RateHeterogeneity*> initial_site_rates;
    
    /**
     is the tree weights fixed
     */
    bool isTreeWeightFixed;

    /**
     is nested openmp
     */
    bool isNestedOpenmp;

    /**
     variables for the shared RHAS model
     */
    double *rhas_var;
};

#endif /* iqtreemix_h */
