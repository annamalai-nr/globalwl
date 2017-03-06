/**********************************************************************
 * Copyright (C) 2017 Christopher Morris <christopher.morris@udo.edu>
 *
 * This file is part of globalwl.
 *
 * globalwl can not be copied and/or distributed without the express
 * permission of Christopher Morris.
 *********************************************************************/

#include <chrono>
#include <cstdio>
#include "src/AuxiliaryMethods.h"
#include "src/WeisfeilerLehmanThreeLocal.h"
#include "src/WeisfeilerLehmanThreeGlobal.h"

using namespace chrono;
using namespace std;

int main() {
    string graph_database_name = "ENZYMES";
    printf("%s\n", graph_database_name.c_str());
    GraphDatabase graph_database = AuxiliaryMethods::read_graph_txt_file(graph_database_name);

    // Kernel: 3-LWL
    // Use sampling: False
    // Use labels: True
    // Initial label: Isomorphism type of 3-tuple
    WeisfeilerLehmanThreeLocal::WeisfeilerLehmanThreeLocal crn_0(graph_database);
    GramMatrix gr_0 = crn_0.compute_gram_matrix(3, false, -1, 0.0, true, true);
    AuxiliaryMethods::write_gram_matrix(gr_0, "ENZYMES");

    // Kernel: 3-LWL
    // Use sampling: True (Initial sample size 100, epsilon=0.1)
    // Use labels: True
    // Initial label: Isomorphism type of 3-tuple
    WeisfeilerLehmanThreeLocal::WeisfeilerLehmanThreeLocal crn_1(graph_database);
    GramMatrix gr_1 = crn_1.compute_gram_matrix(3, true, 100, 0.1, true, true);
    AuxiliaryMethods::write_gram_matrix(gr_1, "ENZYMES");

    // Kernel: 3-GWL
    // Use labels: True
    // Initial label: Isomorphism type of 3-tuple
    WeisfeilerLehmanThreeGlobal::WeisfeilerLehmanThreeGlobal crn_2(graph_database);
    GramMatrix gr_2 = crn_2.compute_gram_matrix(3, true, true);
    AuxiliaryMethods::write_gram_matrix(gr_2, "ENZYMES");

    return 0;
}

