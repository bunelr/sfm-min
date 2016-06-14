#ifndef OPTIM
#define OPTIM
#include "graph.hpp"
#include <vector>


typedef std::vector<uint> ordering;
typedef std::vector<double> vec;

class SF{

public:
    virtual double evaluate(const Subset& picked) const= 0;
    virtual uint dimension() const =0;
    double minimize();
};

class Order{
    ordering _ordered_elts;
    vec greedy_vec;

public:
    Order(ordering ordered_elts, const SF *const  problem);
};


#endif
