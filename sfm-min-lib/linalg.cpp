#include "linalg.hpp"
#include <armadillo>
#include <limits>
#include <iostream>
#include <vector>
#include <cassert>

void reduce_form(std::vector<vec>& vectors, std::vector<double>& coefficient){
    uint nb_vectors = vectors.size();
    uint dimension  = vectors[0].size();

    // std::cout << "Nb vectors: "<< nb_vectors << '\n';
    // std::cout << "Dimension: " << dimension  << '\n';

    arma::mat data = arma::ones(dimension+1, nb_vectors);
    for (int i=0; i < nb_vectors; i++) {
        for (int j=0; j < dimension; j++) {
            data(j,i) = vectors[i][j]; // should use batch-insertion constructor
        }
    }

    arma::mat null_base = arma::null(data);

    // Each row correspond to one vector
    // Each col correspond to one dimension of the null_space
    // std::cout << null_base.n_rows << '\n';
    // std::cout << null_base.n_cols << '\n';

    // Use the first dimension of the null space to bring one coefficient to zero
    double smallest_ratio = std::numeric_limits<double>::max();
    bool positive_coeff = true;
    for (uint vec_idx = 0; vec_idx < nb_vectors; vec_idx++) {
        if (null_base(vec_idx, 0)!=0) {
            double ratio = fabs(coefficient[vec_idx] / null_base(vec_idx, 0));
            if (ratio < smallest_ratio) {
                smallest_ratio = ratio;
                positive_coeff = (null_base(vec_idx, 0) > 0);
            }
        }
    }
    // std::cout << "The smallest ratio is: " << smallest_ratio << '\n';
    if (positive_coeff) {
        smallest_ratio *= -1;
    }
    for (uint vec_idx = 0; vec_idx < nb_vectors; vec_idx++) {
        coefficient[vec_idx] += smallest_ratio * null_base(vec_idx, 0);
        if (coefficient[vec_idx]<=std::numeric_limits<double>::epsilon()) {
            // std::cout << "Bringing coeff " << vec_idx << "to zero."<< '\n';
            coefficient[vec_idx] = 0;
        }
    }
}
