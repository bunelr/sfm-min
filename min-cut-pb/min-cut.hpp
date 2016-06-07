#include "optim.hpp"


class MinCut: public SF {
    Graph graph;
public:
    MinCut(std::string path_to_data);
    double evaluate(Subset picked);
};

MinCut* create_min_cut_pb(std::string path_to_data);
