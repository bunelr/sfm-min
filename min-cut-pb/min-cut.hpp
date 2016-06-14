#include "optim.hpp"


class MinCut: public SF {
    Graph graph;
    uint source_node, sink_node;

    // nb_nodes contains also Source and Sink, that are not really
    // elements of the problem
    uint nb_nodes;
public:
    MinCut(std::string path_to_data);
    uint dimension() const;
    double evaluate(const Subset& picked) const;
    double maxflow(Subset& picked) const;
private:
    double cut_from_node(const Node& node, const Subset& picked) const;
};

MinCut* create_min_cut_pb(std::string path_to_data);
