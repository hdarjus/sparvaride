#include <RcppArmadillo.h>
#include "edmonds_karp.h"

int max_flow_sufficient(
    const arma::umat& delta);

// [[Rcpp::export("max_flow_sufficient")]]
int max_flow_sufficient_r(
    const Rcpp::IntegerMatrix& delta_in) {
  arma::umat delta (delta_in.nrow(), delta_in.ncol());
  std::copy(delta_in.cbegin(), delta_in.cend(), delta.begin());
  return max_flow_sufficient(delta);
}

bool counting_rule_holds(
    const arma::umat& delta);

// [[Rcpp::export("variance_is_identified_indicator")]]
bool counting_rule_holds_r(
    const Rcpp::IntegerMatrix& delta_in) {
  arma::umat delta (delta_in.nrow(), delta_in.ncol());
  std::copy(delta_in.cbegin(), delta_in.cend(), delta.begin());  // why copy??
  return counting_rule_holds(delta);
}

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
    const arma::umat& delta) {
  const unsigned int r = delta.n_rows;
  // quick check
  if (delta.n_cols <= 2 * r) {
    return false;
  }
  if (not arma::all(arma::sum(delta, 0) > 0) or
      not arma::all(arma::sum(delta, 1) > 0)) {
    Rf_error("delta has zero rows or zero columns");
  }

  // resort to the vertex covering method
  Graph graph = delta_to_graph(delta);
#ifndef NDEBUG
  const arma::umat delta_recovered = graph_to_delta(graph);
  if (arma::any(arma::any(delta != delta_recovered))) {
    Rf_error("delta and the recovered delta are not equal");
  }
#endif
  const unsigned int max_flow = edmonds_karp(graph);
  return max_flow == r * (2 * r + 1);
}

int max_flow_sufficient(
    const arma::umat& delta) {
  const unsigned int r = delta.n_rows;
  // quick check
  if (delta.n_cols <= 2 * r) {
    return false;
  }
  if (not arma::all(arma::sum(delta, 0) > 0) or
      not arma::all(arma::sum(delta, 1) > 0)) {
    Rf_error("delta has zero rows or zero columns");
  }

  // resort to the vertex covering method
  Graph graph = delta_to_graph(delta);
#ifndef NDEBUG
  const arma::umat delta_recovered = graph_to_delta(graph);
  if (arma::any(arma::any(delta != delta_recovered))) {
    Rf_error("delta and the recovered delta are not equal");
  }
#endif
  return edmonds_karp(graph);
}

