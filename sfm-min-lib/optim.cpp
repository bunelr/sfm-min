#include "optim.hpp"

Order::Order(ordering ordered_elts, const SF *const problem): _ordered_elts(ordered_elts){
    // We also want to compute the associated greedy_vec
    greedy_vec = vec(problem->dimension());
    Subset picked;
    double prev, current;

    prev = problem->evaluate(picked);
    for(const uint& elt : ordered_elts){
        picked.insert(elt);
        current = problem->evaluate(picked);
        greedy_vec.at(elt) = current - prev;
        prev = current;
    }
}

uint Order::at(uint pos) const{
    return _ordered_elts.at(pos);
}

double SF::minimize(){

    // What we are going to maintain:
    std::vector<Order> all_orders;
    std::vector<double> order_weights;
    vec x;

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
    x = first_order.greedy_vec;


    // Start iterating
    while (true) {
        // Build directed graph based on the ordering
        Graph ordering_graph(nb_elements);
        uint node_from, node_to;
        for(const Order& order: all_orders){ // For each of the ordering
            for (uint from=0;from<nb_elements; ++from) {  // Add an edge from each node to all nodes after him
                for (uint to=0; to<nb_elements; ++to) {
                    node_from = order.at(from);
                    node_to = order.at(to);
                    if (not ordering_graph.exist_edge(node_from, node_to)) {
                        Edge new_edge(node_from, node_to, 1);
                        ordering_graph.add_edge(new_edge);
                    }
                }
            }
        }

        // Identify subsets P {u s.t. x(u) > 0}  and N {u s.t. x(u) < 0}
        Subset P,N;
        for (uint i=0; i<nb_elements; ++i) {
            if (x.at(i) > 0) {
                P.insert(i);
            } else if (x.at(i) < 0) {
                N.insert(i);
            }
        }


        // Find a path from P to N
        if (not ordering_graph.path_exist(P, N)) {
            // There is no path from P to N
            // We have found the minimiser.

            // Find all vertices that have a path to N.
            // This is easy to find because this is all of N, and some
            // of the elements that are neither in P nor N.
            // This corresponds to the solution.
            break;
        } else {
            // There is a path from P to N
            // We will do updates to the ordering

            // Compute the distance d(u) from P to each element
            // For each ordering, find t -> latest in ordering that has max d(u)
            //                    find s -> (s,t) is in the graph and d(s) = d(t)-1, maximally big in the ordering
            // Take the ordering with maximal difference between s and t
            // Generate new orderings, where in each, one element from between s and t goes just in front of s (no other changes)

            // Find the convex combination of greedy-vector of new_ordering that gives \mu ( \delta_t - \delta_s)

            // Find y

            // Find the new value of x -> x' by taking the point closest to y that has x'(t) =< 0

            // Adjust the convex combination to have less than nb_elements Order

        }
        break; // Safety break during development
    }
    return 2;
}
