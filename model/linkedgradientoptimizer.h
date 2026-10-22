/*
 * linkedgradientoptimizer.h
 *
 * The linked step of -S --model-joint with --analytical-gradients (mixture_maker):
 * the exchangeabilities, profiles and weights of a model shared by many partitions
 * are optimised by BFGS on the summed log-likelihood, with the summed analytic
 * gradient instead of finite differences. Branch lengths and rate models stay with
 * their partitions (the alternation in PartitionModel::optimizeParameters).
 *
 * Design: docs/analytical-gradients-design.md, section 17.
 */

#ifndef LINKEDGRADIENTOPTIMIZER_H
#define LINKEDGRADIENTOPTIMIZER_H

#include <memory>
#include <string>
#include <vector>
#include "utils/optimization.h"

class PhyloSuperTree;
class PhyloTree;
class ModelParamMap;

class LinkedGradientOptimizer : public Optimization {
public:
    /** the partitions of `stree` whose model is named `model_name` share that model */
    LinkedGradientOptimizer(PhyloSuperTree *stree, const std::string &model_name);
    virtual ~LinkedGradientOptimizer();

    /**
     * Entry checks, run before every linked step while the copies are unfixed:
     * every partition is supported by the analytic engine and maps to the same
     * parameter layout. Returns "" when the step can run, otherwise a one-line
     * reason (the caller then uses finite differences). Stops with an error if
     * the partitions' copies of the linked model differ.
     */
    std::string prepare();

    /** one linked step from the current model; returns the summed log-likelihood */
    double step(bool write_info, double logl_epsilon);

    /**
     * --ag-gradient-check-only under -S --model-joint: the summed analytic gradient at
     * the current point against central differences of the summed log-likelihood
     * (h and h/2, Richardson); writes <prefix>.gradcheck.tsv and one GRADCHECK line.
     * @return the number of failed parameters
     */
    int gradientCheck();

    /**
     * Writes the shared model as a NEXUS models block reusable with -mdef: for a
     * K-class mixture with linked exchangeabilities `model MM<K>Q` (PAML lower
     * triangle, pinned entry 1, frequencies = weighted mean profile), `frequency
     * MM<K>pi1..K` (normalised, by decreasing weight) and `model MM<K> =
     * MM<K>Q+FMIX{...}`; for one class `model MM1`. 17 significant digits.
     */
    void writeModelFile(const std::string &file) const;

    virtual int getNDim();
    virtual double targetFunk(double x[]);
    virtual double derivativeFunk(double x[], double dfx[]);
    virtual bool stopEarly(int iter, double f_prev, double f_new);
    virtual bool restartParameters(double guess[], int ndim, double lower[], double upper[],
                                   bool bound_check[], int iteration) { return false; }

private:
    /** the linked model's natural parameters (all rate entries, raw frequencies, weights) */
    struct RawState {
        std::vector<std::vector<double>> rates, freqs;
        std::vector<double> weights;
    };
    void snapshot(RawState &st) const;
    /** writes `st` into every partition's copy */
    void restore(const RawState &st);
    /** stops with an error if any partition's copy differs from partition 0's */
    void checkCopiesIdentical() const;
    /** sum over partitions of the current log-likelihood, in partition order */
    double sumLikelihood();

    PhyloSuperTree *stree_;
    std::string name_;
    std::vector<PhyloTree*> trees_;                       // partitions sharing the model, in partition order
    std::vector<int> order_;                              // indices into trees_, largest partition first
    std::vector<std::unique_ptr<ModelParamMap>> maps_;
    std::vector<double> theta_, lh_part_;
    std::vector<std::vector<double>> g_part_;
    double logl_epsilon_ = 0.1, delta_max_ = 0.0;
    int stop_iter_ = 0, steps_ = 0;
    long n_lh_ = 0, n_grad_ = 0, n_fd_ = 0;
    double t_lh_ = 0.0, t_grad_ = 0.0, last_end_ = 0.0;
    bool warned_fd_ = false;
};

#endif /* LINKEDGRADIENTOPTIMIZER_H */
