#include <iostream>
#include "min-cut.hpp"
#include <cassert>

MinCut::MinCut(std::string path_to_data){
    graph = Graph(path_to_data);
    uint nb_nodes = graph.nodes.size();
    source_node = nb_nodes-2;
    sink_node = nb_nodes-1;
}

double MinCut::evaluate(const Subset& picked) const{
    // Evaluate the value of the cut proposed
    double cut_value = 0;

    // A cut is determined by the sum of all the outgoing arcs from U

    // We are not allowed to pick the sink in U
    assert(picked.find(sink_node)==picked.end());

    // And the source necessarily belongs to U
    Node current_node = graph.nodes[source_node];
    cut_value += cut_from_node(current_node, picked);

    for (const auto& node_idx: picked) {
        current_node = graph.nodes[node_idx];
        cut_value += cut_from_node(current_node, picked);
    }

    return cut_value;
}

double MinCut::cut_from_node(const Node& node, const Subset& picked) const{
    double cut_value = 0;
    for (const auto& edge: node.outarcs) {
        if (picked.find(edge.to)==picked.end()) {
            cut_value += edge.capacity;
        }
    }
    return cut_value;
}

MinCut* create_min_cut_pb(std::string path_to_data){
    MinCut* mc_pb = new MinCut(path_to_data);
    return mc_pb;
}
