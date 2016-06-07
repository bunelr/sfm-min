#ifndef OPTIM
#define OPTIM
#include <unordered_set>
#include "graph.hpp"

typedef std::unordered_set<uint> Subset;


class SF{

public:
    virtual double evaluate(Subset picked) = 0;
};

#endif
