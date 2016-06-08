#include <vector>
#include <fstream>
#include <sstream>
#include <string>
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
