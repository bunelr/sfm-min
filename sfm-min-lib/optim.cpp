#include "optim.hpp"
#include <limits>
#include <cassert>

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

uint Order::nb_intermediary(uint from, uint to) const{
    uint s,t;
    uint i=0;
    for (const uint& elt: _ordered_elts) {
        if (elt == from) {
            s = i;
        } else if (elt==to) {
            t = i;
        }
        ++i;
    }
    return t - s;
}

std::vector<ordering> generate_new_orderings(uint from, uint to) const{
    // Find the position of s and t in the ordering
    uint s_pos, t_pos;
    uint i=0;
    for (const uint& elt: _ordered_elts) {
        if (elt == from) {
            s_pos = i;
        } else if (elt==to) {
            t_pos = i;
        }
        ++i;
    }

    // Now that we have both, create our (empty) vector of ordering
    std::vector<ordering> new_orderings;
    for (uint pos_of_moved_element=s_pos+1;
         pos_of_moved_element < t_pos+1;
         pos_of_moved_element++) {
        // Do some copy
    }


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
                std::vector<uint> d = ordering_graph.distance_from(P);
                // For each ordering, find t -> biggest that has max d(u)
                //                    find s -> (s,t) is in the graph and d(s) = d(t)-1, biggest
                uint to_ignore = std::numeric_limits<uint>::max();
                uint max_d=0;
                uint val;
                uint s, t;
                for (uint i=0; i < nb_elements; i++) {
                    val = d[i];
                    if (val!=to_ignore) {
                        if (val>max_d) {
                            max_d = val;
                            s = t; // d[s] might not be d[t]-1 here but
                            // will be eventually because s now
                            // currently holds the second biggest
                            // ever seen and we know that there
                            // exists some with value one-less than
                            // d[t].
                            t = i;
                        } else if (val == max_d) {
                            t = i;
                        } else if (val == (max_d - 1)) {
                            s = i;
                        }
                    }
                }
                assert(d[t]==(d[s]+1));
                // Take the ordering with maximal difference between s and t
                uint max_nb_intermediary=0;
                uint argmax_nb_intermediary=0;
                uint nb_intermediary;
                for (uint i = 0; i < all_orders.size(); i++) {
                    nb_intermediary = all_orders[i].nb_intermediary(s, t);
                    if (nb_intermediary > max_nb_intermediary) {
                        max_nb_intermediary = nb_intermediary;
                        argmax_nb_intermediary = i;
                    }
                }
                assert(max_nb_intermediary>0);

                // Generate new orderings, where in each, one element from
                // between s and t goes just in front of s (no other
                // changes)
                std::vector<ordering> new_orderings = all_orders[argmax_nb_intermediary].generate_new_orderings(s, t);
                // Convert them to Order to get their associated greedy vectors


                // Find the convex combination of greedy-vector of new_ordering that gives \mu ( \delta_t - \delta_s)

                // Find y

                // Find the new value of x -> x' by taking the point closest to y that has x'(t) =< 0

                // Adjust the convex combination to have less than nb_elements Order

            }
            break; // Safety break during development
        }
        return 2;
    }
