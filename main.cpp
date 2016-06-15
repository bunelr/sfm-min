#include <iostream>
#include "min-cut.hpp"


int main(int argc, char *argv[])
{
    // Min-cut, knowing that it is a min-cut,
    // we solve it by max-flow
    MinCut* mc_pb = create_min_cut_pb(argv[1]);
    Subset picked;
    double max_value = mc_pb->maxflow(picked);
    // std::cout << "Value computed by Dinic Algorithm: " << max_value << '\n';
    // std::cout << "Evaluation by the oracle: " << mc_pb->evaluate(picked) << '\n';
    // display_subset(picked);

    // Min-cut, as a submodular function that we don't know
    // Using Schrijver's algorithm
    double min_value = mc_pb->minimize();
    // std::cout << "Value computed by Schrijver's algorithm: " << min_value << '\n';

    return 0;
}
