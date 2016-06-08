#include <iostream>
#include "min-cut.hpp"
#include "optim.hpp"


int main(int argc, char *argv[])
{
    MinCut* mc_pb = create_min_cut_pb(argv[1]);
    Subset picked = {0};
    double cut_value = mc_pb->evaluate(picked);
    std::cout << cut_value << '\n';
    return 0;
}
