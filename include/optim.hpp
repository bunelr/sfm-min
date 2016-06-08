#ifndef OPTIM
#define OPTIM
#include <unordered_set>
#include "graph.hpp"

typedef std::unordered_set<uint> Subset;


class SF{

public:
    virtual double evaluate(const Subset& picked) const= 0;
};

#endif
