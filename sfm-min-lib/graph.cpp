#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <string>
#include <vector>
#include "graph.hpp"

void display_subset(const Subset& picked){
        for (const uint& node: picked) {
                std::cout << node  << ' ';
        }
        std::cout << '\n';
}


Node::Node(uint x, std::vector<Edge> outarcs):id(x), outarcs(outarcs) {}

Edge::Edge(uint from, uint to, double weight): from(from), to(to), capacity(weight) {}

Graph::Graph(uint nb_nodes){
        for (uint i=0; i<nb_nodes; i++) {
                std::vector<Edge> new_outarcs;
                Node new_node = Node(i, new_outarcs);
                nodes.push_back(new_node);
        }
}

void Graph::add_edge(Edge edge){
        for (const Edge& existing_edge: nodes.at(edge.from).outarcs) {
                // We don't want to allow to create edges that already
                // exist
                assert(existing_edge.to != edge.to);// Maybe should replace this with exceptions
        }
        nodes.at(edge.from).outarcs.push_back(edge);
}

bool Graph::exist_edge(uint from, uint to){
        for (const Edge& existing_edge: nodes.at(from).outarcs) {
                if (existing_edge.to == to) {
                        return true;
                }
        }
        return false;
}

Graph::Graph(std::string path){

        std::ifstream file(path.c_str());
        std::string all_edge_weight_str, edge_weight_str;
        uint current_from = 0;
        uint current_to = 0;
        double edge_weight;

        while(getline(file, all_edge_weight_str)){
                std::vector<Edge> outarcs;
                std::stringstream ss(all_edge_weight_str);
                current_to = 0;
                while(getline(ss, edge_weight_str, ' ')){
                        edge_weight = std::stod(edge_weight_str);
                        if(edge_weight>0) {
                                outarcs.push_back(Edge(current_from, current_to, edge_weight));
                        }
                        current_to++;
                }
                Node  new_node = Node(current_from, outarcs);
                nodes.push_back(new_node);
                current_from++;
        }
}

bool Graph::shortest_path(uint source_node, uint target_node, std::vector<uint>& path) const {
        // Perform a breadth first search to find the shortest path
        // from source_node to target_node
        // Returns a boolean indicating whether a path was found
        // Put the index of the Node in the path in the `path` reference
        uint unvisited = std::numeric_limits<uint>::max();
        std::vector<uint> father_node(nodes.size(), unvisited);
        std::queue<uint> to_explore;
        std::vector<Edge> outarcs;
        bool found_target_node = false;
        uint current_node;
        to_explore.push(source_node);


        while(not found_target_node){
                if (to_explore.empty()) {
                        break;
                }
                current_node = to_explore.front();
                to_explore.pop();

                outarcs = nodes[current_node].outarcs;
                for (const Edge& edge: outarcs) {
                        if (father_node[edge.to]==unvisited) {
                                father_node[edge.to] = edge.from;
                                if (edge.to==target_node) {
                                        found_target_node = true;
                                        break;
                                } else {
                                        to_explore.push(edge.to);
                                }
                        }
                }
        }

        if (found_target_node) {
                // There exists a path from the source to the target
                uint node = target_node;
                path.push_back(node);
                while(node!=source_node){
                        node = father_node[node];
                        path.push_back(node);
                }
                std::reverse(path.begin(), path.end());
                return true;
        }
        return false;
}

void Graph::pass_flow(const std::vector<uint>& path, double flow_value){
        // Adapt the residual graph according to the new flow going through path

        uint path_index = 0;
        uint from = path[path_index];
        uint to = path[path_index+1];
        uint end_node_idx = path.back();

        while (from != end_node_idx) {
                // Find the arc to look at it's capacity
                for (std::vector<Edge>::iterator edge = nodes[from].outarcs.begin();
                     edge!=nodes[from].outarcs.end(); ++edge) {
                        if (edge->to==to) {
                                // Found the correct edge
                                // Add the reverse edge
                                // Make sure that we don't duplicate edges
                                bool already_existing_edge = false;
                                for (Edge& rev_edge: nodes[to].outarcs){
                                        if(rev_edge.to == from){
                                                // The edge already exists, just need to augment its value
                                                rev_edge.capacity += flow_value;
                                                already_existing_edge = true;
                                                break;
                                        }
                                }
                                // The reverse edge does not already exist, we can just create it.
                                if (not already_existing_edge) {
                                        Edge new_edge(to, from, flow_value);
                                        nodes[to].outarcs.push_back(new_edge);
                                }


                                if (edge->capacity == flow_value) {
                                        // If the edge is saturated, remove it
                                        nodes[from].outarcs.erase(edge);
                                } else {
                                        // Otherwise, just lower the capacity
                                        assert(flow_value < edge->capacity);
                                        edge->capacity -= flow_value;
                                }

                                // Go handle the next edge
                                path_index++;
                                from = to;
                                to = path[path_index+1];
                                break;
                        }
                }
        }

}

void Graph::reachable(uint source_node, Subset& picked) const{
        // Perform a breadth first search to find all the reachable node from the sources
        // Put the index of all reachable nodes in 'picked'
        std::vector<bool> visited(nodes.size(), false);
        std::queue<uint> to_explore;
        std::vector<Edge> outarcs;
        uint current_node;
        to_explore.push(source_node);


        while(true){
                if (to_explore.empty()) {
                        break;
                }
                current_node = to_explore.front();
                to_explore.pop();

                outarcs = nodes[current_node].outarcs;
                for (const Edge& edge: outarcs) {
                        if (not visited[edge.to]) {
                                visited[edge.to] = true;
                                to_explore.push(edge.to);
                                picked.insert(edge.to);
                        }
                }
        }
}

bool Graph::path_exist(const Subset& from, const Subset& to){
        std::vector<bool> visited(nodes.size(), false);
        std::queue<uint> to_explore;
        std::vector<Edge> outarcs;
        uint current_node;

        for(const uint& src: from){
                visited[src] = true;
                to_explore.push(src);
        }

        while(true){
                if (to_explore.empty()) {
                        break;
                }
                current_node = to_explore.front();
                to_explore.pop();

                outarcs = nodes[current_node].outarcs;
                for (const Edge& edge: outarcs) {
                        if (not visited[edge.to]) {
                                if (to.count(edge.to)>0) { // We can reach an element of the target subset
                                        return true;
                                } else{
                                        visited[edge.to] = true;
                                        to_explore.push(edge.to);
                                }
                        }
                }
        }
        return false;
}


std::vector<uint> Graph::distance_from(const Subset& from){
        uint max_value = std::numeric_limits<uint>::max();
        std::vector<uint> d(nodes.size(), max_value);
        std::queue<uint> to_explore;
        std::vector<Edge> outarcs;
        uint current_node;

        for(const uint& src: from){
                d[src] = 0;
                to_explore.push(src);
        }

        while(true){
                if (to_explore.empty()) {
                        break;
                }
                current_node = to_explore.front();
                to_explore.pop();

                outarcs = nodes[current_node].outarcs;
                for (const Edge& edge: outarcs) {
                        if (d[edge.to] == max_value) {
                                d[edge.to] = d[edge.from] + 1;
                                to_explore.push(edge.to);
                        }
                }
        }
        return d;
}
