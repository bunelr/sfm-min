#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include "graph.hpp"

Node::Node(uint x, std::vector<Edge> outarcs) {
        id = x;
        outarcs = outarcs;
}


Edge::Edge(uint from, uint to, double weight) {
        from = from;
        to = to;
        capacity = weight;
}


Graph::Graph(std::string path){

        std::ifstream file(path.c_str());
        edge_nb = 0;

        std::string nb_vertices_str;
        getline(file, nb_vertices_str);

        std::string nb_voisins_str;
        unsigned short int current_node =0;
        while(getline(file, nb_voisins_str)){
                int nb_voisins = stoi(nb_voisins_str, nullptr, 10);
                edge_nb += nb_voisins;
                std::vector<Edge> outarcs;

                for(int i=0; i<nb_voisins;++i){
                        std::string neighb_id_str;
                        std::string edge_weight_str;

                        getline(file, neighb_id_str);
                        getline(file, edge_weight_str);

                        uint neighb_id = stoi(neighb_id_str);
                        double edge_weight = stod(edge_weight_str);

                        outarcs.push_back(Edge(current_node, neighb_id, edge_weight));
                        //cout<< "added Edge from "<<current_node<<" to "<<neighb_id<<endl;
                }
                Node  new_node = Node(current_node, outarcs);
                nodes.push_back(new_node);
                //cout<< "Created the node "<< new_node.id<<endl;
                current_node++;
        }
}
