#include "optim.hpp"


class MinCut: public SF {
    Graph graph;
    uint source_node, sink_node;
public:
    MinCut(std::string path_to_data);
    double evaluate(const Subset& picked) const;
    double maxflow() const;
private:
    double cut_from_node(const Node& node, const Subset& picked) const;
};

MinCut* create_min_cut_pb(std::string path_to_data);
