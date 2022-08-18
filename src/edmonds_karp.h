#ifndef EDMONDS_KARP_EFI_H
#define EDMONDS_KARP_EFI_H

#include <vector>

using Vertex = unsigned int;

struct Edge {
  Vertex start,
         end;
  int capacity,
      flow;
  Edge* reverse = nullptr;
};

using EdgeList = std::vector<Edge>;
using Graph = std::vector<EdgeList>;

unsigned int edmonds_karp (Graph& graph);

#endif //EDMONDS_KARP_EFI_H
