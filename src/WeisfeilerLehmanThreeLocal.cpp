/**********************************************************************
 * Copyright (C) 2017 Christopher Morris <christopher.morris@udo.edu>
 *
 * This file is part of globalwl.
 *
 * globalwl can not be copied and/or distributed without the express
 * permission of Christopher Morris.
 *********************************************************************/

#include "AuxiliaryMethods.h"
#include "WeisfeilerLehmanThreeLocal.h"

namespace WeisfeilerLehmanThreeLocal {
    WeisfeilerLehmanThreeLocal::WeisfeilerLehmanThreeLocal(const GraphDatabase &graph_database) : m_graph_database(
            graph_database), m_label_to_index(), m_num_labels(0) {}

    GramMatrix WeisfeilerLehmanThreeLocal::compute_gram_matrix(const uint num_iterations, const bool use_sampling,
                                                               const uint num_samples, const double eps,
                                                               const bool use_labels,
                                                               bool use_iso_type) {
        vector<ColorCounter> color_counters;
        color_counters.reserve(m_graph_database.size());

        if (!use_sampling) {
            for (auto &graph: m_graph_database) {
                if (!use_sampling) {
                    color_counters.push_back(compute_colors(graph, num_iterations, use_labels, use_iso_type));
                }
            }
        } else {
            size_t num_graphs = m_graph_database.size();
            ColorCounter color_map;
            color_counters.resize(m_graph_database.size(), color_map);

            size_t max_vertices = 0;
            size_t n = 0;

            for (const Graph &g: m_graph_database) {
                n = g.get_num_nodes();
                if (n > max_vertices) {
                    max_vertices = n;
                }
            }

            double sum_samples = num_samples;
            double new_samples = sum_samples;
            double epsilon = 0.0;

            double num_samples_ub = log(max_vertices * num_graphs * (1.0 / 0.01)) / (eps * eps);
            do {
                int c = 0;
                for (Graph &graph: m_graph_database) {
                    color_counters[c] = compute_colors_sample_adaptive(graph, num_iterations, new_samples,
                                                                       color_counters[c], use_labels);
                    c++;
                }

                epsilon = 2.0 * sqrt(sum_samples) * sqrt(2.0 * log(sum_samples)) / sum_samples +
                          sqrt(0.1 / (2.0 * sum_samples));


                new_samples = new_samples * 2.0;
                sum_samples += new_samples;
            } while ((epsilon > eps) and (sum_samples < num_samples_ub));
        }

        ulong num_graphs = m_graph_database.size();
        vector<S> nonzero_compenents;

        for (ulong i = 0; i < num_graphs; ++i) {
            ColorCounter c = color_counters[i];

            for (const auto &j: c) {
                Label key = j.first;
                uint value = j.second;
                uint index = m_label_to_index.find(key)->second;
                nonzero_compenents.push_back(move(S(i, index, value)));
            }
        }

        GramMatrix feature_vectors(num_graphs, m_num_labels);
        feature_vectors.setFromTriplets(nonzero_compenents.begin(), nonzero_compenents.end());

        if (use_sampling) {
            feature_vectors = feature_vectors * (1.0 / m_num_labels);
        }

        GramMatrix gram_matrix(num_graphs, num_graphs);
        gram_matrix = feature_vectors * feature_vectors.transpose();

        return gram_matrix;
    }


    ColorCounter
    WeisfeilerLehmanThreeLocal::compute_colors(const Graph &g, const uint num_iterations, const bool use_labels,
                                               const bool use_iso_type) {
        size_t num_nodes = g.get_num_nodes();

        unordered_map<tuple<uint, uint, uint>, uint> triple_to_int;
        vector<Triple> triples;

        // Generate all three element sets over the nodes of "g".
        size_t num_triples = 0;
        for (Node i = 0; i < num_nodes; ++i) {
            for (Node j = 0; j < num_nodes; ++j) {
                for (Node k = 0; k < num_nodes; ++k) {
                    triples.push_back(make_tuple(i, j, k));
                    triple_to_int.insert({{make_tuple(i, j, k), num_triples}});
                    num_triples++;
                }
            }
        }

        Labels coloring;
        coloring.reserve(num_triples);
        Labels coloring_temp;
        coloring_temp.reserve(num_triples);

        Labels labels;
        if (use_labels) {
            labels = g.get_labels();
        }

        ColorCounter color_map;
        // Assign isomorphism type to each 3-element set.
        for (Triple t: triples) {
            Node i = get<0>(t);
            Node j = get<1>(t);
            Node k = get<2>(t);

            Label new_color;
            if (use_labels) {
                Label c_i = labels[i];
                Label c_j = labels[j];
                Label c_k = labels[k];

                if (use_iso_type) {

                    Labels labels(
                            {{AuxiliaryMethods::pairing(g.get_degree(i) + 1, c_i + 1),
                                     AuxiliaryMethods::pairing(g.get_degree(j) + 1, c_j + 1),
                                     AuxiliaryMethods::pairing(g.get_degree(k) + 1, c_k + 1)
                             }
                            });

                    sort(labels.begin(), labels.end());

                    new_color = g.has_edge(i, j) + g.has_edge(i, k) + g.has_edge(j, k);
                    for (Label d: labels) {
                        new_color = AuxiliaryMethods::pairing(new_color, d);
                    }
                } else {
                    Labels labels({{
                                           AuxiliaryMethods::pairing(AuxiliaryMethods::pairing(c_i, c_j),
                                                                     g.has_edge(i, j)),
                                           AuxiliaryMethods::pairing(AuxiliaryMethods::pairing(c_i, c_k),
                                                                     g.has_edge(i, k)),
                                           AuxiliaryMethods::pairing(AuxiliaryMethods::pairing(c_j, c_k),
                                                                     g.has_edge(j, k))
                                   }
                                  });
                    sort(labels.begin(), labels.end());

                    new_color = 1;
                    for (Label d: labels) {
                        new_color = AuxiliaryMethods::pairing(new_color, d);
                    }
                }
            } else {
                new_color = g.has_edge(i, j) + g.has_edge(i, k) + g.has_edge(j, k);
            }

            coloring[triple_to_int.find(t)->second] = new_color;

            ColorCounter::iterator it(color_map.find(new_color));
            if (it == color_map.end()) {
                color_map.insert({{new_color, 1}});
                m_label_to_index.insert({{new_color, m_num_labels}});
                m_num_labels++;
            } else {
                it->second++;
            }
        }

        uint h = 1;
        while (h <= num_iterations) {
            for (Triple t: triples) {
                Node i = get<0>(t);
                Node j = get<1>(t);
                Node k = get<2>(t);


                Labels colors;
                // Get colors of neighbors.
                // Exchange node 0.
                Nodes neighbors_j = g.get_neighbours(j);
                for (Node &c: neighbors_j) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(
                            make_tuple(c, j, k));


                    // FIX: Do I really need this triple_to_int her???
                    if ((it != triple_to_int.end()) and (c != i)) {
                        colors.push_back(coloring[it->second]);
                    }
                }

                Nodes neighbors_k = g.get_neighbours(k);
                for (Node &c: neighbors_k) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(
                            make_tuple(c, j, k));
                    if (!g.has_edge(c, j)) {
                        if ((it != triple_to_int.end()) and (c != i)) {
                            colors.push_back(coloring[it->second]);
                        }
                    }
                }

                // Exchange node 1.
                Nodes neighbors_i = g.get_neighbours(i);
                for (Node &c: neighbors_i) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(
                            make_tuple(i, c, k));
                    if ((it != triple_to_int.end()) and (c != j)) {
                        colors.push_back(coloring[it->second]);
                    }
                }

                neighbors_k = g.get_neighbours(k);
                for (Node &c: neighbors_k) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(
                            make_tuple(i, c, k));
                    if (!g.has_edge(c, i)) {
                        if ((it != triple_to_int.end()) and (c != j)) {
                            colors.push_back(coloring[it->second]);
                        }
                    }
                }

                // Exchange node 2.
                neighbors_i = g.get_neighbours(i);
                for (Node &c: neighbors_i) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(
                            make_tuple(i, j, c));
                    if ((it != triple_to_int.end()) and (c != k)) {
                        colors.push_back(coloring[it->second]);
                    }
                }

                neighbors_j = g.get_neighbours(j);
                for (Node &c: neighbors_j) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(
                            make_tuple(i, j, c));
                    if (!g.has_edge(c, i)) {
                        if ((it != triple_to_int.end()) and (c != k)) {
                            colors.push_back(coloring[it->second]);
                        }
                    }
                }

                sort(colors.begin(), colors.end());
                Label new_color = coloring[triple_to_int.find(make_tuple(i, j, k))->second];

                for (Label const &c: colors) {
                    new_color = AuxiliaryMethods::pairing(new_color, c);
                }

                coloring_temp[triple_to_int.find(t)->second] = new_color;

                // Keep track how often "new_color" occurs.
                ColorCounter::iterator it(color_map.find(new_color));
                if (it == color_map.end()) {
                    color_map.insert({{new_color, 1}});
                    m_label_to_index.insert({{new_color, m_num_labels}});
                    m_num_labels++;
                } else {
                    it->second++;
                }
            }

            coloring = coloring_temp;
            h++;
        }

        return color_map;
    }


    ColorCounter WeisfeilerLehmanThreeLocal::compute_colors_sample_adaptive(const Graph &g, const uint num_iterations,
                                                                            const uint num_samples,
                                                                            ColorCounter color_map,
                                                                            bool use_labels) {
        size_t num_nodes = g.get_num_nodes();

        random_device rand_dev;
        mt19937 mt(rand_dev());
        uniform_int_distribution<Node> uniform_node_sampler(0, num_nodes - 1);

        Labels coloring;
        for (uint c = 0; c < num_samples; ++c) {
            // Sample a triple.
            Node t_0 = uniform_node_sampler(mt);
            Node t_1 = uniform_node_sampler(mt);
            Node t_2 = uniform_node_sampler(mt);
            Triple triple = make_tuple(t_0, t_1, t_2);

            unordered_map<Triple, uint> depth;
            unordered_map<Triple, uint> triple_to_int;
            vector<Triple> neighborhood = explore_neighborhood(g, triple, num_iterations, depth, triple_to_int);

            coloring.resize(neighborhood.size());

            Labels labels;
            if (use_labels) {
                labels = g.get_labels();
            }

            Label new_color;
            for (Triple t: neighborhood) {
                Node i = get<0>(t);
                Node j = get<1>(t);
                Node k = get<2>(t);

                if (use_labels) {

                    Label c_i = labels[i];
                    Label c_j = labels[j];
                    Label c_k = labels[k];

                    Labels labels({{
                                           AuxiliaryMethods::pairing(AuxiliaryMethods::pairing(c_i, c_j),
                                                                     g.has_edge(i, j)), AuxiliaryMethods::pairing(
                                    AuxiliaryMethods::pairing(c_i, c_k), g.has_edge(i, k)), AuxiliaryMethods::pairing(
                                    AuxiliaryMethods::pairing(c_j, c_k), g.has_edge(j, k))
                                   }
                                  });
                    sort(labels.begin(), labels.end());

                    new_color = 1;
                    for (Label d: labels) {
                        new_color = AuxiliaryMethods::pairing(new_color, d);
                    }
                } else {
                    new_color = g.has_edge(i, j) + g.has_edge(i, k) + g.has_edge(j, k);
                }
                coloring[triple_to_int.find(t)->second] = new_color;

                ColorCounter::iterator it(color_map.find(new_color));
                if (it == color_map.end()) {
                    color_map.insert({{new_color, 1}});
                    m_label_to_index.insert({{new_color, m_num_labels}});
                    m_num_labels++;
                } else {
                    it->second++;
                }
            }

            Labels coloring_temp;
            coloring_temp = coloring;

            uint h = 1;
            while (h <= num_iterations) {
                for (Triple &v: neighborhood) {
                    Node i = get<0>(v);
                    Node j = get<1>(v);
                    Node k = get<2>(v);
                    Labels colors;

                    // Get colors of neighbors.
                    // Exchange node 0.
                    Nodes neighbors_j = g.get_neighbours(j);
                    for (Node &c: neighbors_j) {
                        unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(c, j, k));
                        if ((it != triple_to_int.end()) and (c != i)) {
                            colors.push_back(coloring[it->second]);
                        }
                    }

                    Nodes neighbors_k = g.get_neighbours(k);
                    for (Node &c: neighbors_k) {
                        unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(c, j, k));
                        if (!g.has_edge(c, j)) {
                            if ((it != triple_to_int.end()) and (c != i)) {
                                colors.push_back(coloring[it->second]);
                            }
                        }
                    }

                    // Exchange node 1.
                    Nodes neighbors_i = g.get_neighbours(i);
                    for (Node &c: neighbors_i) {
                        unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(i, c, k));
                        if ((it != triple_to_int.end()) and (c != j)) {
                            colors.push_back(coloring[it->second]);
                        }
                    }

                    neighbors_k = g.get_neighbours(k);
                    for (Node &c: neighbors_k) {
                        unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(i, c, k));
                        if (!g.has_edge(c, i)) {
                            if ((it != triple_to_int.end()) and (c != j)) {
                                colors.push_back(coloring[it->second]);
                            }
                        }
                    }

                    // Exchange node 2.
                    neighbors_i = g.get_neighbours(i);
                    for (Node &c: neighbors_i) {
                        unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(i, j, c));
                        if ((it != triple_to_int.end()) and (c != k)) {
                            colors.push_back(coloring[it->second]);
                        }
                    }

                    neighbors_j = g.get_neighbours(j);
                    for (Node &c: neighbors_j) {
                        unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(i, j, c));
                        if (!g.has_edge(c, i)) {
                            if ((it != triple_to_int.end()) and (c != k)) {
                                colors.push_back(coloring[it->second]);
                            }
                        }
                    }
                    Label new_color = coloring[triple_to_int.find(v)->second];

                    sort(colors.begin(), colors.end());
                    Label new_color_0;
                    if (colors.size() != 0) {
                        new_color_0 = colors.back();
                        colors.pop_back();
                        for (Label const c: colors) {
                            new_color_0 = AuxiliaryMethods::pairing(new_color_0, c);
                        }
                        new_color = AuxiliaryMethods::pairing(new_color, new_color_0);
                    }

                    coloring_temp[triple_to_int.find(v)->second] = new_color;

                    if (h + depth.find(v)->second <= num_iterations) {
                        // Keep track how often "new_label" occurs.
                        ColorCounter::iterator it(color_map.find(new_color));
                        if (it == color_map.end()) {
                            color_map.insert({{new_color, 1}});
                            m_label_to_index.insert({{new_color, m_num_labels}});
                            m_num_labels++;
                        } else {
                            it->second++;
                        }
                    }
                }
                h++;
                coloring = coloring_temp;
            }
        }

        return color_map;
    }

    inline vector<Triple>
    WeisfeilerLehmanThreeLocal::explore_neighborhood(const Graph &g, const Triple &triple, const uint num_iterations,
                                                     unordered_map<Triple, uint> &depth,
                                                     unordered_map<Triple, uint> &triple_to_int) {
        unordered_set<Triple> visited;
        queue<Triple> queue;
        queue.push(triple);
        visited.insert(triple);
        vector<Triple> neighbourhood;
        depth.insert({{triple, 0}});

        uint u = 0;
        while (!queue.empty()) {
            Triple q(queue.front());
            queue.pop();

            Node i = get<0>(q);
            Node j = get<1>(q);
            Node k = get<2>(q);

            neighbourhood.push_back(q);
            triple_to_int.insert({{q, u}});
            u++;
            uint current_depth = depth.find(q)->second;

            if (current_depth < num_iterations) {
                // Exchange node 0.
                Nodes neighbors_j = g.get_neighbours(j);
                vector<Triple> neighbours;
                for (Node &c: neighbors_j) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(c, j, k));
                    if (it == triple_to_int.end() and (c != i)) {
                        neighbours.push_back(make_tuple(c, j, k));
                        triple_to_int.insert({{make_tuple(c, j, k), u}});
                        u++;
                    }
                }

                Nodes neighbors_k = g.get_neighbours(k);
                for (Node &c: neighbors_k) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(c, j, k));
                    if (!g.has_edge(c, j)) {
                        if (it == triple_to_int.end() and (c != i)) {
                            neighbours.push_back(make_tuple(c, j, k));
                            triple_to_int.insert({{make_tuple(c, j, k), u}});
                            u++;
                        }
                    }
                }

                // Exchange node 1.
                Nodes neighbors_i = g.get_neighbours(i);
                for (Node &c: neighbors_i) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(i, c, k));
                    if (it != triple_to_int.end() and (c != j)) {
                        neighbours.push_back(make_tuple(i, c, k));
                        triple_to_int.insert({{make_tuple(i, c, k), u}});
                        u++;
                    }
                }

                neighbors_k = g.get_neighbours(k);
                for (Node &c: neighbors_k) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(i, c, k));
                    if (!g.has_edge(c, i)) {
                        if (it != triple_to_int.end() and (c != j)) {
                            neighbours.push_back(make_tuple(i, c, k));
                            triple_to_int.insert({{make_tuple(i, c, k), u}});
                            u++;
                        }
                    }
                }

                // Exchange node 2.
                neighbors_i = g.get_neighbours(i);
                for (Node &c: neighbors_i) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(i, j, c));
                    if (it != triple_to_int.end() and (c != k)) {
                        neighbours.push_back(make_tuple(i, j, c));
                        triple_to_int.insert({{make_tuple(i, j, c), u}});
                        u++;
                    }
                }

                neighbors_j = g.get_neighbours(j);
                for (Node &c: neighbors_j) {
                    unordered_map<Triple, uint>::const_iterator it = triple_to_int.find(make_tuple(i, j, c));
                    if (!g.has_edge(c, i)) {
                        if (it != triple_to_int.end() and (c != k)) {
                            neighbours.push_back(make_tuple(i, j, c));
                            triple_to_int.insert({{make_tuple(i, j, c), u}});
                            u++;
                        }
                    }
                }

                for (Triple &n: neighbours) {
                    if (visited.find(n) == visited.end()) {
                        depth.insert({{n, current_depth + 1}});
                        queue.push(n);
                        visited.insert(n);
                    }
                }
            }
        }

        return neighbourhood;
    }

    WeisfeilerLehmanThreeLocal::~WeisfeilerLehmanThreeLocal() {}
}
