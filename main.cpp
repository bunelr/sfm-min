#include <iostream>
#include "min-cut.hpp"
#include "optim.hpp"


int main(int argc, char *argv[])
{
    MinCut* mc_pb = create_min_cut_pb(argv[1]);
    std::cout << mc_pb->evaluate(Subset()) << '\n';
    return 0;
}
