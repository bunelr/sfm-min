#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cassert>
#include "optim.hpp"
#include "linalg.hpp"

ordering generate_random_ordering(uint size){
    ordering initial_ordering(size);
    for (uint i=0; i < size; i++) {
        initial_ordering.at(i) = i;
    }
    std::random_shuffle(initial_ordering.begin(), initial_ordering.end());
    return initial_ordering;
}

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

void Order::show_ordering() const{
    for (const uint& node: _ordered_elts ) {
        std::cout << node << ' ';
    }
    std::cout << '\n';

}

uint Order::at(uint pos) const{
    return _ordered_elts.at(pos);
}

uint Order::pos_in_ordering(uint node_idx) const {
    uint i=0;
    for (const uint& elt: _ordered_elts) {
        if (elt == node_idx) {
            return i;
        }
        ++i;
    }
}

uint Order::nb_intermediary(uint from, uint to) const{
    uint s_pos = pos_in_ordering(from);
    uint t_pos = pos_in_ordering(to);
    if (t_pos < s_pos) {
        return 0;
    }
    return t_pos - s_pos;
}

vec Order::arrange_vector_according_to_ordering(const vec& node_idx_sorted) const{
    vec ordering_sorted(node_idx_sorted.size());
    for (uint i=0; i < node_idx_sorted.size(); i++) {
        ordering_sorted[i] = node_idx_sorted[_ordered_elts[i]];
    }
    return ordering_sorted;
}

std::vector<ordering> Order::generate_new_orderings(uint from, uint to) const{
    // Find the position of s and t in the ordering
    uint s_pos = pos_in_ordering(from);
    uint t_pos = pos_in_ordering(to);

    // Now that we have both, create our (empty) vector of ordering
    std::vector<ordering> new_orderings;
    for (uint pos_of_moved_element=s_pos+1;
         pos_of_moved_element < t_pos+1;
         pos_of_moved_element++) {
        ordering new_ord(_ordered_elts.size(), 0);
        std::copy(_ordered_elts.begin(), _ordered_elts.begin()+s_pos,
                  new_ord.begin());
        new_ord[s_pos] = _ordered_elts[pos_of_moved_element];
        std::copy(_ordered_elts.begin()+s_pos, _ordered_elts.begin()+pos_of_moved_element,
                  new_ord.begin()+s_pos+1);
        std::copy(_ordered_elts.begin()+pos_of_moved_element+1, _ordered_elts.end(),
                  new_ord.begin()+pos_of_moved_element+1);
        assert(new_ord.size()==_ordered_elts.size());
        new_orderings.push_back(new_ord);
        // for (const uint& node: new_ord) {
        //     std::cout << node  << ' ';
        // }
        // std::cout << '\n' << '\n';

    }
    return new_orderings;

}

double SF::minimize(Subset& picked){

    // Useful constants
    uint nb_elements = dimension();

    // What we are going to maintain:
    std::vector<Order> all_orders;
    std::vector<double> order_weights;
    vec x, old_x, new_x;

    // Useful for the invariant
    bool invariant_pass;
    std::vector<uint> old_d, d(nb_elements, 0);
    Graph old_ordering_graph, ordering_graph;
    Subset P,N, neutral, old_P, old_N, old_neutral;
    uint old_dt, dt=std::numeric_limits<uint>::max();
    uint old_t, t=std::numeric_limits<uint>::max();
    uint old_s, s=std::numeric_limits<uint>::max();
    uint old_alpha, alpha=std::numeric_limits<uint>::max();
    uint old_beta, beta=std::numeric_limits<uint>::max();


    // Initialisation - Create a first dumb ordering
    ordering initial_ordering = generate_random_ordering(nb_elements);
    Order first_order(initial_ordering, this);
    all_orders.push_back(first_order);
    order_weights.push_back(1);
    old_x = first_order.greedy_vec;
    x = first_order.greedy_vec;


    // Start iterating
    while (true) {
        // Move the current value to the old_* container.
        old_x = x;
        old_d = d;
        old_ordering_graph = ordering_graph;
        old_dt = dt;
        old_t = t;
        old_s = s;
        old_alpha = alpha;
        old_beta = beta;
        old_P = P;
        old_N = N;
        old_neutral = neutral;



        // Build directed graph based on the ordering
        ordering_graph = Graph(nb_elements);
        uint node_from, node_to;
        for(const Order& order: all_orders){ // For each of the ordering
            for (uint from=0;from<nb_elements; ++from) {  // Add an edge from each node to all nodes after him
                for (uint to=from+1; to<nb_elements; ++to) {
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
        P.clear();N.clear();neutral.clear();
        for (uint i=0; i<nb_elements; ++i) {
            if (x.at(i) > 0) {
                P.insert(i);
            } else if (x.at(i) < 0) {
                N.insert(i);
            } else {
                neutral.insert(i);
            }
        }


        // Find a path from P to N
        if (not ordering_graph.path_exist(P, N)) {
            // There is no path from P to N
            // We have found the minimiser.


            // Find all vertices that have a path to N.
            // This is easy to find because this is all of N, and some
            // of the elements that are neither in P nor N.

            // For every element of Neutral, check whether it can reach N
            Subset dummy_subset;
            for (const uint& src: neutral) {
                dummy_subset = {src};
                assert(dummy_subset.size()==1);
                if (ordering_graph.path_exist(dummy_subset, N)) {
                    N.insert(src);
                }
            }

            // Compute the value by summming all x(N) now
            picked = N;
            return evaluate(picked);
            // TODO
            break;
        } else {
            // There is a path from P to N
            // We will do updates to the ordering

            // Compute the distance d(u) from P to each element
            d = ordering_graph.distance_from(P);
            // For each ordering, find t -> biggest that has max d(u)
            uint to_ignore = std::numeric_limits<uint>::max();
            uint max_d=0;
            uint val;
            for (uint i=0; i < nb_elements; i++) {
                if (N.find(i)==N.end()) { // Only interested in finding it in N
                    continue;
                }
                val = d[i];
                if (val!=to_ignore) {
                    if (val>max_d) {
                        max_d = val;
                        t = i;
                    } else if (val == max_d) {
                        t = i;
                    }
                }
            }
            dt = d[t];
            //                    find s -> (s,t) is in the graph and d(s) = d(t)-1, biggest
            for (uint i=0; i < nb_elements; i++) { // PERF: look backward instead of forward -> earlier stopping.
                val = d[i];
                if (val == (max_d-1)) {
                    bool st_in_graph;
                    if (ordering_graph.exist_edge(i, t)) {
                        s = i;
                    }
                }
            }
            // std::cout << "s is: " << s << ", Distance is " << d[s]<< '\n';
            // std::cout << "t is: " << t << ", Distance is " << d[t]<< '\n';

            assert(d[t]==(d[s]+1));
            // Take the ordering with maximal difference between s and t
            uint max_nb_intermediary=0;
            uint nb_max_intermediary=0;
            uint argmax_nb_intermediary=0;
            uint nb_intermediary;
            //std::cout << "All the intermediaries:" << '\n';
            for (uint i = 0; i < all_orders.size(); i++) {
                nb_intermediary = all_orders[i].nb_intermediary(s, t);
                //std::cout << nb_intermediary << ' ';
                assert(nb_intermediary < nb_elements);
                // I assumed that if s > t, the number of intermediary
                // is considered to be zero.

                if (nb_intermediary > max_nb_intermediary) {
                    max_nb_intermediary = nb_intermediary;
                    argmax_nb_intermediary = i;
                    nb_max_intermediary = 1;
                } else if (nb_intermediary == max_nb_intermediary ) {
                    nb_max_intermediary++;
                }
            }
            //std::cout << '\n';
            alpha = max_nb_intermediary;
            beta = nb_max_intermediary;



            // Generate new orderings, where in each, one element from
            // between s and t goes just in front of s (no other
            // changes)
            std::vector<ordering> new_orderings = all_orders[argmax_nb_intermediary].generate_new_orderings(s, t);
            // Convert them to Order to get their associated greedy vectors
            std::vector<Order> new_orders;
            for (const ordering& ord: new_orderings) {
                Order new_order(ord, this);
                new_orders.push_back(new_order);
            }

            double delta;
            std::vector<double> new_ord_weights(new_orders.size());
            // Find the convex combination of greedy-vector of new_ordering that gives
            //                greedy_of_ordering + \delta ( I_t - I_s)
            // Generate the matrix, the rows correspond to the orders, indexed by u in the appropriate order.
            // The colons contains the greedy vec, in the appropriate order
            Order& base_order = all_orders[argmax_nb_intermediary];
            vec ref = base_order.arrange_vector_according_to_ordering(base_order.greedy_vec);
            uint s_pos = base_order.pos_in_ordering(s);
            uint t_pos = base_order.pos_in_ordering(t);
            uint nb_elem_considered = (t_pos+1 - s_pos);
            std::vector<double> working_matrix(new_orders.size() * nb_elem_considered);
            bool zero_valid_delta = false;
            uint row_zero_valid_delta;
            // std::cout << "s_pos " << s_pos << " \t t_pos " << t_pos << '\n';
            // std::cout << "Nb orders " << new_orders.size() << "\t nb_elem: " << nb_elem_considered << '\n';
            for (uint i=0; i < new_orders.size(); i++) {
                vec order_greedy = base_order.arrange_vector_according_to_ordering(new_orders[i].greedy_vec);
                // Check for degenerate case
                uint diag_term = (i+1) + s_pos;
                if (order_greedy[diag_term]==ref[diag_term]) {
                    zero_valid_delta = true;
                    row_zero_valid_delta = i;
                    break;
                }
                for (uint j=s_pos; j < t_pos+1; j++) {
                    uint adj_index = j - s_pos;
                    working_matrix[i*nb_elem_considered + adj_index] = order_greedy[j] - ref[j];
                }
                // for (uint j=s_pos; j < t_pos+1; j++){
                //     uint adj_index = j - s_pos;
                //     std::cout << working_matrix[i*nb_elem_considered +adj_index] << '\t';
                // }
                // std::cout << '\n';
            }
            if (zero_valid_delta) {
                new_ord_weights[row_zero_valid_delta] = 1;
                delta = 0;
                // Verify that this means that they had the same greedy vector
                // std::cout << "Zero is a valid delta" << '\n';
                vec& ref = base_order.greedy_vec;
                vec& rep = new_orders[row_zero_valid_delta].greedy_vec;
                for (uint i=0; i<ref.size(); i++) {
                    assert(abs(ref[i]-rep[i]) < std::numeric_limits<double>::epsilon());
                }
            } else{
                std::vector<double> current(nb_elem_considered, 0);
                for (int i=new_orders.size()-1; i >= 0; i--) {
                    double target;
                    if (i==new_orders.size()-1) { // Refactor this
                        target = 1;
                    } else{
                        target = 0;
                    }
                    double diag_term = working_matrix[i*nb_elem_considered + (i+1)];

                    new_ord_weights[i] = (target - current[(i+1)]) / diag_term;
                    for (uint elem=0; elem < nb_elem_considered; elem++) {
                        double comp = working_matrix[i*nb_elem_considered + elem];
                        if (elem > (i+1)) {
                            assert(comp==0);
                        }
                        current[elem] += comp * new_ord_weights[i];
                    }
                }
                // for (const double& node: current) {
                //     std::cout << node  << ' ';
                // }
                // std::cout << '\n';
                assert(abs(current[0]- (-1)) < std::numeric_limits<double>::epsilon());
                assert(abs(current[current.size()-1] -1) < std::numeric_limits<double>::epsilon());
                for (uint i=1; i < current.size()-1; i++) {
                    // if (current[i] != 0) {
                    //     for (const double& node: current) {
                    //         std::cout << node  << ' ';
                    //     }
                    //     std::cout << '\n';
                    // }
                    assert(abs(current[i] - 0)< std::numeric_limits<double>::epsilon());
                }
                delta = 0;
                for (const double& w: new_ord_weights) {
                    assert(w>=0);
                    delta += w;
                }
                delta = 1/delta; // We want a convex combinations so all things must sum to one
                for (double& w: new_ord_weights) {
                    w *= delta;
                }
            }
            // std::cout << "Delta is: " << delta << '\n';
            // std::cout << "Weights on the new orderings are: ";
            // for (const double& node: new_ord_weights) {
            //     std::cout << node  << ' ';
            // }
            // std::cout << '\n';


            // Find y
            vec y = x;
            double lambda_1 = order_weights[argmax_nb_intermediary];
            y[s] -= lambda_1 * delta;
            y[t] += lambda_1 * delta;

            // Find the new value of x -> x' by taking the point closest to y that has x'(t) =< 0
            double partial_multiplier = 1;
            if (y[t] <= 0) {
                // std::cout << "Gone full y" << '\n';
                new_x = y;
                order_weights[argmax_nb_intermediary] = 0;
            } else {
                new_x = x;
                new_x[s] += new_x[t];
                partial_multiplier = -new_x[t] / (delta*lambda_1);
                new_x[t] = 0;
                // Because we are only doing a partial application of our update vector,
                // we need to keep the original ordering, to preserve convexity of the combination.
                // (As opposed to the other branch where we can get rid of it)
                order_weights[argmax_nb_intermediary] = lambda_1 * (1-partial_multiplier);
            }


            // Add the new ordering and the related weights to the maintained list
            for (Order new_ord: new_orders) {
                all_orders.push_back(new_ord);
            }
            for (double new_weight: new_ord_weights) {
                order_weights.push_back(new_weight * lambda_1 * partial_multiplier);
            }

            // Check that we haven't made any mistake and that we still have x as a convex combination.
            double sum_order_weights = 0;
            for (double& w: order_weights) {
                sum_order_weights += w;
            }
            // std::cout << "Error on the sum of the weights: " << sum_order_weights-1 << '\n';
            // std::cout << "Smallest admissible error: " << std::numeric_limits<double>::epsilon() << '\n';
            assert(abs(sum_order_weights-1)< std::numeric_limits<double>::epsilon());
            vec reconstructed_x(nb_elements, 0);
            for (uint i=0; i < all_orders.size(); i++) {
                double weight = order_weights[i];
                vec& greedy_vec = all_orders[i].greedy_vec;
                for (uint j = 0; j < nb_elements; j++) {
                    reconstructed_x[j] += weight * greedy_vec[j];
                }
            }

            for (uint i=0; i < nb_elements; i++) {
                if (abs(reconstructed_x[i] - new_x[i]) > std::numeric_limits<double>::epsilon()) {
                    std::cout << "Error at position " << i <<
                        ": "<<  reconstructed_x[i] <<
                        " versus " << new_x[i] << '\n';
                }
                assert(abs(reconstructed_x[i] - new_x[i]) < std::numeric_limits<double>::epsilon());
            }

            // We now have a valid combination of greedy vectors associated with our orderings.
            // Let's make it diet a bit

            // All the orderings that have zero weight in the current linear combination should be removed
            uint ordering_index = 0;
            while (ordering_index < order_weights.size()) {
                if (order_weights[ordering_index]==0) {
                    all_orders.erase(all_orders.begin() + ordering_index);
                    order_weights.erase(order_weights.begin() + ordering_index);
                } else {
                    ordering_index++;
                }
            }

            while (all_orders.size() > nb_elements) {
                // This is horribly inefficient
                // We reduce the dimension by 1 every time, while we could do much more.
                // We probably would want to reduce it as much as possible.
                std::vector<vec> all_greedys;
                for (Order ord: all_orders) {
                    all_greedys.push_back(ord.greedy_vec);
                }
                reduce_form(all_greedys, order_weights);

                // All the orderings that have been set to zero
                uint ordering_index = 0;
                bool removed_one = false;
                while (ordering_index < order_weights.size()) {
                    if (order_weights[ordering_index]==0) {
                        all_orders.erase(all_orders.begin() + ordering_index);
                        order_weights.erase(order_weights.begin() + ordering_index);
                        removed_one = true;
                        break;
                    } else {
                        ordering_index++;
                    }
                }
                assert(removed_one);
            }


            // Verification code that nothing problematic happened here
            // while we were removing redundant ordering
            std::fill(reconstructed_x.begin(), reconstructed_x.end(), 0);
            for (uint i=0; i < all_orders.size(); i++) {
                double weight = order_weights[i];
                vec& greedy_vec = all_orders[i].greedy_vec;
                for (uint j = 0; j < nb_elements; j++) {
                    reconstructed_x[j] += weight * greedy_vec[j];
                }
            }
            for (uint i=0; i < nb_elements; i++) {
                if (abs(reconstructed_x[i] - new_x[i]) > std::numeric_limits<double>::epsilon()) {
                    std::cout << "Error at position " << i <<
                        ": "<<  reconstructed_x[i] <<
                        " versus " << new_x[i] << '\n';
                }
                assert(abs(reconstructed_x[i] - new_x[i]) < std::numeric_limits<double>::epsilon());
            }




        }

        // Invariant checking
        // std::cout << "(dt,t, s, a, b)" << '\n';
        std::cout << "(" <<dt << ", " <<
            t << ", " <<
            s << ", " <<
            alpha << ", "<<
            beta << ")\n";
        // std::cout << "x: [";
        // for (double x_val: old_x) {
        //     std::cout << x_val << ',';
        // }; std::cout << "]\n";

        bool distance_improvement = false;
        // No distance regression
        for (uint elt= 0; elt < nb_elements; elt++) {
            assert(d[elt] >= old_d[elt]);
            if (d[elt] > old_d[elt]) {
                distance_improvement = true;
                invariant_pass = true;
            }
        }
        if (not distance_improvement) {
            if(invariant_pass){
                // I'm not entirely sure of the offset
                invariant_pass=false;
            } else {
                assert(x[t]<0);
                assert(old_x[t]<0 or t==old_s);
                assert(old_d[s]<old_d[t]);
                assert(old_d[t] <= old_d[old_t]);
                if (old_d[t]==old_d[old_t]) {
                    assert(t <= old_t);
                }
                if (old_d[t]==old_d[old_t] and t == old_t) {
                    // Assert s, old_t in A'
                    assert(ordering_graph.exist_edge(s, old_t));
                    // Assert s, old_t was in A
                    assert(old_ordering_graph.exist_edge(s, old_t));
                    // From it, we conclude that
                    assert(s <= old_s);
                }


                // Lexicographical ordering should hold
                assert(dt <= old_dt);
                if (dt == old_dt) {
                    assert(t <= old_t);
                    if (t == old_t) {
                        assert(s <= old_s);
                        if (s == old_s) {
                            assert(alpha <= old_alpha);
                            if (alpha == old_alpha) {
                                assert(beta < old_beta);
                            }
                        }
                    }
                }
            }

        } else {
            std::cout << "Distance improvement" << '\n';
        }
        x = new_x;
    }
}
