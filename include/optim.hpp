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
public:
    vec greedy_vec;

    Order(ordering ordered_elts, const SF *const  problem);
    uint nb_intermediary(uint from, uint to) const;
    uint at(uint pos) const;
    std::vector<ordering> generate_new_orderings(uint from, uint to) const;
    vec arrange_vector_according_to_ordering(const vec& node_idx_sorted) const;
    uint pos_in_ordering(uint node_idx) const;
};


#endif
