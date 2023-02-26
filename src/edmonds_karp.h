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
