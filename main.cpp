#include <iostream>
#include "min-cut.hpp"


int main(int argc, char *argv[])
{
    MinCut* mc_pb = create_min_cut_pb(argv[1]);
    Subset picked;
    double max_value = mc_pb->maxflow(picked);
    std::cout << "Value computed by Dinic Algorithm: " << max_value << '\n';
    std::cout << "Evaluation by the oracle: " << mc_pb->evaluate(picked) << '\n';
    display_subset(picked);
    return 0;
}
