#include "min-cut.hpp"

MinCut::MinCut(std::string path_to_data){
    graph = Graph(path_to_data);
}

double MinCut::evaluate(Subset picked){
    return 2;
}

MinCut* create_min_cut_pb(std::string path_to_data){
    MinCut* mc_pb = new MinCut(path_to_data);
    return mc_pb;
}
