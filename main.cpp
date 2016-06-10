#include <iostream>
#include "min-cut.hpp"
#include "optim.hpp"


int main(int argc, char *argv[])
{
    MinCut* mc_pb = create_min_cut_pb(argv[1]);
    Subset picked = {0};
    double max_value = mc_pb->maxflow();
    std::cout << max_value << '\n';
    return 0;
}
