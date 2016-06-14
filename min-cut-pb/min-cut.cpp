#include <iostream>
#include <algorithm>
#include <limits>
#include "min-cut.hpp"
#include <vector>
#include <cassert>

MinCut::MinCut(std::string path_to_data){
    graph = Graph(path_to_data);
    nb_nodes = graph.nodes.size();
    source_node = nb_nodes-2;
    sink_node = nb_nodes-1;
}

uint MinCut::dimension() const{
    return nb_nodes-2;
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

    for (const uint& node_idx: picked) {
        current_node = graph.nodes[node_idx];
        cut_value += cut_from_node(current_node, picked);
    }

    return cut_value;
}

double MinCut::cut_from_node(const Node& node, const Subset& picked) const{
    double cut_value = 0;
    for (const Edge& edge: node.outarcs) {
        if (picked.find(edge.to)==picked.end()) {
            cut_value += edge.capacity;
        }
    }
    return cut_value;
}

double MinCut::maxflow(Subset& picked) const {
    // Find the maxflow, using Dinic's algorithm
    double flow_value = 0;

    Graph residual_graph = graph;
    bool finding_st_path = true;

    std::vector<uint> st_path;
    while(true){
        // Finding the shortest st path
        st_path.clear();
        finding_st_path = residual_graph.shortest_path(source_node, sink_node, st_path);
        if(not finding_st_path){
            break;
        }

        // Find the maximum allowable flow
        double max_passable_flow = std::numeric_limits<double>::max();
        uint from_index = 0;
        uint from = st_path[from_index];
        uint to = st_path[from_index+1];
        Node from_node = residual_graph.nodes[from];
        while (from != sink_node) {
            // Find the arc to look at it's capacity
            for (const Edge& edge: from_node.outarcs) {
                if (edge.to==to) {
                    max_passable_flow = std::min(max_passable_flow, edge.capacity);
                    from_index++;
                    from = to;
                    from_node = residual_graph.nodes[from];
                    to = st_path[from_index+1];
                    break;
                }
            }
        }
        // Update the flow value and the residual graph
        residual_graph.pass_flow(st_path, max_passable_flow);

        flow_value += max_passable_flow;
    }

    // Find all the elements that are in S
    residual_graph.reachable(source_node, picked);


    return flow_value;
}


MinCut* create_min_cut_pb(std::string path_to_data){
    MinCut* mc_pb = new MinCut(path_to_data);
    return mc_pb;
}
