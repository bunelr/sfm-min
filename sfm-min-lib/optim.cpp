#include "optim.hpp"

Order::Order(ordering ordered_elts, const SF *const problem): _ordered_elts(ordered_elts){
    // We also want to compute the associated greedy_vec
    greedy_vec = vec(problem->dimension());
    Subset picked;
    double prev, current;

    prev = problem->evaluate(picked);
    for(auto& elt : ordered_elts){
        picked.insert(elt);
        current = problem->evaluate(picked);
        greedy_vec.at(elt) = current - prev;
        prev = current;
    }
}

double SF::minimize(){

    // What we are going to maintain:
    std::vector<Order> all_orders;
    std::vector<double> order_weights;

    // Useful constants
    uint nb_elements = dimension();

    // Initialisation - Create a first dumb ordering
    ordering initial_ordering(nb_elements);
    for (uint i=0; i < nb_elements; i++) {
        initial_ordering.at(i) = i;
    }
    Order first_order(initial_ordering, this);
    all_orders.push_back(first_order);
    order_weights.push_back(1);


    return 2;
}
