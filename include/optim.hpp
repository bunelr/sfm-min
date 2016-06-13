#ifndef OPTIM
#define OPTIM
#include "graph.hpp"


class SF{

public:
    virtual double evaluate(const Subset& picked) const= 0;
};

#endif
