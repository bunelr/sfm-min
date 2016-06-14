#ifndef GRAPH
#define GRAPH
#include<vector>
#include <unordered_set>

typedef std::unordered_set<uint> Subset;

void display_subset(const Subset&);

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
        // Properties
        std::vector<Node> nodes;
        uint edge_nb;

        // Constructors
        Graph(){};
        Graph(std::string path);
        Graph(uint nb_nodes);

        void add_edge(Edge edge);

        // Simple queries
        bool exist_edge(uint from, uint to);

        // Algorithmic Methods
        bool shortest_path(uint source_node, uint target_node, std::vector<uint>& path) const;
        void reachable(uint source_node, Subset& picked) const;
        void pass_flow(const std::vector<uint>& path, double flow_value);
        bool path_exist(const Subset& from, const Subset& to);
};

#endif
