#include <RcppArmadillo.h>
#include "edmonds_karp.h"

static const int infinity = std::numeric_limits<int>::max();

int max_flow_necessary(
    const arma::umat& delta);

// [[Rcpp::export("max_flow_necessary")]]
int max_flow_necessary_r(
    const Rcpp::IntegerMatrix& delta_in) {
  arma::umat delta (delta_in.nrow(), delta_in.ncol());
  std::copy(delta_in.cbegin(), delta_in.cend(), delta.begin());
  return max_flow_necessary(delta);
}

bool necessary_condition_holds(
    const arma::umat& delta);

// [[Rcpp::export("variance_may_be_identified_indicator")]]
bool necessary_condition_holds_r(
    const Rcpp::IntegerMatrix& delta_in) {
  arma::umat delta (delta_in.nrow(), delta_in.ncol());
  std::copy(delta_in.cbegin(), delta_in.cend(), delta.begin());  // why copy??
  return necessary_condition_holds(delta);
}

// we always create the representation of a dense delta and do not optimize for graph size
// @param delta binary matrix with dimensions factors x observations
// @return network graph; variables are in column-major ordering (of transposed delta)
//   and equations are in lexicographical ordering of column indices (of transposed delta)
static
Graph delta_to_graph(
    const arma::umat& delta) {
  const unsigned int n_variables = delta.n_rows * delta.n_cols,
                     n_equations = (delta.n_cols * (delta.n_cols - 1u)) / 2u;
  // allocate memory
  Graph graph(2 + n_variables + n_equations);
  graph[0].reserve(n_variables);
  for (unsigned int i = 0; i < n_variables; i++) {
    const arma::uvec index = arma::ind2sub(arma::size(delta), i);
    const unsigned int row = index(0),
                       col = index(1);
    graph[1 + i].reserve(delta(row, col) ? arma::accu(delta.row(row)) : 1u);
  }
  // (one_col, other_col) is the lexicographical ordering
  for (unsigned int one_col = 0; one_col < delta.n_cols; one_col++) {
    for (unsigned int other_col = one_col + 1u; other_col < delta.n_cols; other_col++) {
      const unsigned int i = (one_col * delta.n_cols + other_col) - (one_col + 1u + ((one_col + 1u) * one_col) / 2u);
      graph[1 + n_variables + i].reserve(1u + 2u * arma::accu(delta.col(one_col) % delta.col(other_col)));
    }
  }
  graph[1 + n_variables + n_equations].reserve(n_equations);
  // source - variables
  EdgeList& edges_source = graph[0];
  for (unsigned int i = 0; i < n_variables; i++) {
    EdgeList& edges_neighbor = graph[i + 1];
    edges_source.push_back({0, i + 1, 1, 0});
    edges_neighbor.push_back({i + 1, 0, 0, 0});

    edges_source.back().reverse = &edges_neighbor.back();
    edges_neighbor.back().reverse = &edges_source.back();
  }
  // variables - equations
  for (unsigned int one_col = 0; one_col < delta.n_cols; one_col++) {
    for (unsigned int other_col = one_col + 1u; other_col < delta.n_cols; other_col++) {
      const unsigned int index_equation = (1u + n_variables + one_col * delta.n_cols + other_col) - (one_col + 1u + ((one_col + 1u) * one_col) / 2u);
      EdgeList& edges_equation = graph[index_equation];
      for (unsigned int row = 0; row < delta.n_rows; row++) {
        if (delta(row, one_col) and delta(row, other_col)) {
          const unsigned int index_one_variable = 1u + one_col * delta.n_rows + row,
                             index_other_variable = 1u + other_col * delta.n_rows + row;
          EdgeList& edges_one_variable = graph[index_one_variable],  // simplify notation
                  & edges_other_variable = graph[index_other_variable];  // simplify notation

          edges_one_variable.push_back({index_one_variable, index_equation, infinity, 0});
          edges_equation.push_back({index_equation, index_one_variable, 0, 0});

          edges_one_variable.back().reverse = &edges_equation.back();
          edges_equation.back().reverse = &edges_one_variable.back();

          edges_other_variable.push_back({index_other_variable, index_equation, infinity, 0});
          edges_equation.push_back({index_equation, index_other_variable, 0, 0});

          edges_other_variable.back().reverse = &edges_equation.back();
          edges_equation.back().reverse = &edges_other_variable.back();
        }
      }
    }
  }
  // equations - target
  EdgeList& edges_target = graph.back();
  for (unsigned int i = 0; i < n_equations; i++) {
    const unsigned int index_equation = n_variables + i + 1;
    const unsigned int index_target = n_variables + n_equations + 1;
    EdgeList& edges_equation = graph[index_equation];
    edges_equation.push_back({index_equation, index_target, 1, 0});
    edges_target.push_back({index_target, index_equation, 0, 0});

    edges_equation.back().reverse = &edges_target.back();
    edges_target.back().reverse = &edges_equation.back();
  }
  return graph;
}

static
arma::umat graph_to_delta(
    const Graph& graph) {
  const unsigned int n_variables = graph[0].size(),
                     n_equations = graph.size() - 2 - n_variables,
                     m = std::round(0.5 * (1 + std::sqrt(1 + 8. * n_equations))),
                     r = n_variables / m;
  arma::umat delta(r, m, arma::fill::zeros);
  // validate variable capacities
  for (unsigned int i = 0; i < n_variables; i++) {
    if (graph[0][i].capacity != 1) {
      Rf_error("capacities to variables are incorrect");
    }
  }
  // validate and store delta edges
  for (unsigned int i = 0; i < n_variables; i++) {
    const unsigned int index_variable = i + 1u;
    const arma::uvec index_in_delta = arma::ind2sub(arma::size(delta), index_variable - 1u);
    delta(index_in_delta(0), index_in_delta(1)) = graph[index_variable].size() > 1;
    for (const Edge edge : graph[index_variable]) {
      if (edge.capacity < infinity and edge.end > 0) {
        Rcpp::Rcout << edge.start << ' ' <<
          edge.end << "; " <<
          edge.capacity << ' ' <<
          edge.flow <<
          std::endl;
        Rf_error("capacities between variables and equations incorrect");
      }
    }
  }
  // validate edges to the sink/target node
  for (unsigned int i = 0; i < n_equations; i++) {
    const unsigned int index_equation = i + n_variables + 1;
    if (graph[index_equation].back().capacity != 1) {
      Rf_error("incorrect capacities from equations");
    }
    if (graph[index_equation].back().end != graph.size() - 1) {
      Rf_error("edge end incorrect (equations)");
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

bool necessary_condition_holds(
    const arma::umat& delta) {
  // quick check
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
  return max_flow == arma::accu(delta);
}

int max_flow_necessary(
    const arma::umat& delta) {
  // quick check
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

