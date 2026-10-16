/*
 * agtrace.h
 * --ag-trace: write the optimisation history of a tree search to
 * <prefix>.agtrace.jsonl (one JSON object per line) for
 * test_scripts/ag/viz/visualizer.html. Observation only: reads parameters
 * through getters and reuses tree strings IQ-TREE already built.
 */

#ifndef AGTRACE_H
#define AGTRACE_H

#include <string>

class PhyloTree;

namespace AgTrace {
    extern bool on;
    /** open the trace for the main tree search (no-op unless --ag-trace) */
    void open(PhyloTree *tree);
    void close();
    /** start of a phase, e.g. "search" */
    void phase(PhyloTree *tree, const char *name);
    /** one counted tree-search iteration */
    void iteration(PhyloTree *tree, int it, double lnl, const std::string &tree_str);
    /** end of IQTree::optimizeModelParameters */
    void callEnd(PhyloTree *tree, int it, double lnl, const std::string &tree_str);
    /** one EM level or polish round inside GradientOptimizer::optimize */
    void round(PhyloTree *tree, const char *kind, int round, double lnl);
}

#define AG_TRACE(x) do { if (AgTrace::on) AgTrace::x; } while (0)

#endif
