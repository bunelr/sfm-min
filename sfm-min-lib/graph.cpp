#include <algorithm>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include "graph.hpp"

Node::Node(uint x, std::vector<Edge> outarcs):id(x), outarcs(outarcs) {}

Edge::Edge(uint from, uint to, double weight): from(from), to(to), capacity(weight) {}


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

bool Graph::shortest_path(uint source_node, uint target_node, std::vector<uint>& path){
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
                for (const auto& edge: outarcs) {
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
