/**********************************************************************
 * Copyright (C) 2017 Christopher Morris <christopher.morris@udo.edu>
 *
 * This file is part of globalwl.
 *
 * globalwl can not be copied and/or distributed without the express
 * permission of Christopher Morris.
 *********************************************************************/

#ifndef WLFAST_WEISFEILERLEHMANTHREELOCAL_H
#define WLFAST_WEISFEILERLEHMANTHREELOCAL_H

#include <cmath>
#include <unordered_map>
#include <queue>

#include "Graph.h"

using Triple = tuple<Node, Node, Node>;

using namespace GraphLibrary;

namespace WeisfeilerLehmanThreeLocal {
    class WeisfeilerLehmanThreeLocal {
    public:
        WeisfeilerLehmanThreeLocal(const GraphDatabase &graph_database);

        // Compute Gram matrix for the 3-LWL.
        GramMatrix compute_gram_matrix(const uint num_iterations, const bool use_sampling, const uint num_samples,
                                       const double eps,
                                       const bool use_labels, const bool use_iso_type);

        ~WeisfeilerLehmanThreeLocal();

    private:

        // Compute lables for each graph in graph database.
        ColorCounter
        compute_colors(const Graph &g, const uint num_iterations, const bool use_labels, const bool use_iso_type);

        // Compute labels for each graph in graph database using sampling.
        ColorCounter compute_colors_sample_adaptive(const Graph &g, const uint num_iterations, const uint num_samples,
                                                    ColorCounter color_counter, const bool use_labels);

        ColorCounter
        compute_colors_sample_adaptive_parallel(const Graph &g, const uint num_iterations, const uint num_samples,
                                                ColorCounter color_counter, const bool use_labels);


        // Get neighborhood of a node.
        inline vector<Triple> explore_neighborhood(const Graph &g, const Triple &triple, const uint num_iterations,
                                                   unordered_map<Triple, uint> &depth,
                                                   unordered_map<Triple, uint> &triple_to_int);

        GraphDatabase m_graph_database;
        ColorCounter m_label_to_index;
        int m_num_labels;
    };
}

#endif //WLFAST_WEISFEILERLEHMANTHREELOCAL_H
