/*
 * R package sparvaride by
 *     Darjus Hosszejni and Sylvia Frühwirth-Schnatter Copyright (C) 2023-
 *
 *  This file is part of the R package sparvaride: Variance Identification
 *  in Sparse Factor Analysis
 *
 *  The R package sparvaride is free software: you can redistribute it
 *  and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 2 or
 *  any later version of the License.
 *
 *  The R package sparvaride is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with the R package sparvaride. If that is not the case, please
 *  refer to <http://www.gnu.org/licenses/>.
 */


#include <RcppArmadillo.h>
#include "edmonds_karp.h"

//' Verify that the counting rule CR(r,1) holds
//'
//' This is an implementation of the algorithm described in Section 3 of
//' Hosszejni and Fruehwirth-Schnatter (2022). The algorithm is used to verify
//' that the counting rule CR(r,1) holds for the sparsity pattern of the transpose
//' of a factor loading matrix. As detailed in Section 2 of the same paper, if
//' CR(r,1) holds, then the idiosyncratic variances are generically identified.
//' If CR(r,1) does not hold, then we do not know whether the idiosyncratic
//' variances are identified or not.
//'
//' @param delta an `m` x `r` matrix of `0`s and `1`s, where `delta(i,j) == 1` if and only if
//'    the i-th observation loads on the j-th factor
//' @returns `TRUE` if CR(`r`,`1`) holds, `FALSE` otherwise
//' @keywords models multivariate
//' @concept factor analysis variance identification
//' @seealso [stats::factanal()]
//' @references Hosszejni and Fruehwirth-Schnatter (2022). "Cover It Up! Bipartite
//'    Graphs Uncover Identifiability in Sparse Factor Analysis". arXiv:2211.00671.
//'    <doi:10.48550/arXiv.2211.00671>
//' @example inst/examples/counting_rule_holds.R
//' @export
// [[Rcpp::export(name="counting_rule_holds", rng=false)]]
bool counting_rule_holds(
    const arma::mat& delta);

static
Graph delta_to_graph(
    const arma::umat& delta) {
  const unsigned int n_factors = delta.n_rows,
                     n_observations = delta.n_cols;
  // allocate memory
  Graph graph(2 + n_factors + n_observations);
  graph[0].reserve(n_factors);
  for (unsigned int i = 0; i < n_factors; i++) {
    graph[1 + i].reserve(1 + arma::accu(delta.row(i)));
  }
  for (unsigned int i = 0; i < n_observations; i++) {
    graph[1 + n_factors + i].reserve(arma::accu(delta.col(i)) + 1);
  }
  graph[1 + n_factors + n_observations].reserve(n_observations);
  // source
  EdgeList& edges_source = graph[0];
  for (unsigned int i = 0; i < n_factors; i++) {
    EdgeList& edges_neighbor = graph[i + 1];
    edges_source.push_back({0, i + 1, (int)(2 * n_factors + 1), 0});
    edges_neighbor.push_back({i + 1, 0, 0, 0});

    edges_source.back().reverse = &edges_neighbor.back();
    edges_neighbor.back().reverse = &edges_source.back();
  }
  // factor nodes
  for (unsigned int i = 0; i < n_factors; i++) {
    const unsigned int index_factor = i + 1;
    EdgeList& edges_factor = graph[index_factor];  // simplify notation
    const arma::uvec neighbors = arma::find(delta.row(i)) + n_factors + 1;
    for (unsigned int j = 0; j < neighbors.n_elem; j++) {
      EdgeList& edges_neighbor = graph[neighbors(j)];  // simplify notation
      edges_factor.push_back({index_factor, neighbors(j), std::numeric_limits<int>::max(), 0});
      edges_neighbor.push_back({neighbors(j), index_factor, 0, 0});

      edges_factor.back().reverse = &edges_neighbor.back();
      edges_neighbor.back().reverse = &edges_factor.back();
    }
  }
  // observation nodes
  EdgeList& edges_target = graph.back();
  for (unsigned int i = 0; i < n_observations; i++) {
    const unsigned int index_observation = n_factors + i + 1;
    const unsigned int index_target = n_factors + n_observations + 1;
    EdgeList& edges_observation = graph[index_observation];
    edges_observation.push_back({index_observation, index_target, (int)n_factors, 0});
    edges_target.push_back({index_target, index_observation, 0, 0});

    edges_observation.back().reverse = &edges_target.back();
    edges_target.back().reverse = &edges_observation.back();
  }
  return graph;
}

static
arma::umat graph_to_delta(
    const Graph& graph) {
  const unsigned int n_factors = graph[0].size();
  const unsigned int n_observations = graph.size() - 2 - n_factors;
  arma::umat delta(n_factors, n_observations, arma::fill::zeros);
  // validate factor capacities
  for (unsigned int i = 0; i < n_factors; i++) {
    if (graph[0][i].capacity != 2 * (int)n_factors + 1) {
      Rf_error("capacities to factors are incorrect");
    }
  }
  // validate and store delta edges
  for (unsigned int i = 0; i < n_factors; i++) {
    const unsigned int index_factor = i + 1;
    for (unsigned int j = 1; j < graph[index_factor].size(); j++) {
      delta(i, graph[index_factor][j].end - n_factors - 1) = 1u;
      if (graph[index_factor][j].capacity < std::numeric_limits<int>::max()) {
        Rf_error("capacities between factors and observations incorrect");
      }
    }
  }
  // validate edges to the sink/target node
  for (unsigned int i = 0; i < n_observations; i++) {
    const unsigned int index_observation = i + n_factors + 1;
    if (graph[index_observation].back().capacity != (int)n_factors) {
      Rf_error("incorrect capacities from observations");
    }
    if (graph[index_observation].back().end != graph.size() - 1) {
      Rf_error("edge end incorrect (observations)");
    }
  }
  // validate reverse edges
  for (unsigned int i = 0; i < graph.size(); i++) {
    for (unsigned int j = 0; j < graph[i].size(); j++) {
      if (graph[i][j].reverse -> start != graph[i][j].end or
          graph[i][j].reverse -> end != graph[i][j].start) {
        Rcpp::Rcout << i << ' ' << j << std::endl;
        Rcpp::Rcout << graph[i][j].start << ' ' << graph[i][j].end << std::endl;
        Rcpp::Rcout << graph[i][j].reverse -> start << ' ' << graph[i][j].reverse -> end << std::endl;
        Rf_error("incorrect reverse edge pointer");
      }
    }
  }

  return delta;
}

bool counting_rule_holds(
    const arma::mat& delta) {  // input is not transposed
  const arma::umat udelta = arma::trans(delta > 0);
  const unsigned int r = udelta.n_rows;
  // quick check
  if (udelta.n_cols <= 2 * r) {
    return false;
  }
  if (not arma::all(arma::sum(udelta, 0) > 0) or
      not arma::all(arma::sum(udelta, 1) > 0)) {
    Rf_error("delta has zero rows or zero columns");
  }

  // resort to the vertex covering method
  Graph graph = delta_to_graph(udelta);
#ifndef NDEBUG
  const arma::umat udelta_recovered = graph_to_delta(graph);
  if (arma::any(arma::any(udelta != udelta_recovered))) {
    Rf_error("delta and the recovered delta are not equal");
  }
#endif
  const unsigned int max_flow = edmonds_karp(graph);
  return max_flow == r * (2 * r + 1);
}

