#ifndef GRAPH
#define GRAPH
#include<vector>

typedef unsigned int uint;

class Edge {
public:
        uint from;
        uint to;
        double capacity;

        Edge(uint from, uint to, double weight);
        Edge(){};

};


class Node {

public:
        uint id;
        std::vector<Edge> outarcs;

        Node(uint x, std::vector<Edge> outarcs);
        Node(){};
};


class Graph {

public:
        std::vector<Node> nodes;
        uint edge_nb;

        Graph(std::string path);
        Graph(){};

        bool shortest_path(uint source_node, uint target_node, std::vector<uint>& path) const;
        void pass_flow(const std::vector<uint>& path, double flow_value);
};

#endif
