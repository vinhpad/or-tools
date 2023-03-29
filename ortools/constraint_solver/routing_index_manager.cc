// Copyright 2010-2022 Google LLC
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "ortools/constraint_solver/routing_index_manager.h"

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

// không gian sbsl cung cấp flat_hash_set
// nghi ngờ đây là một cây BST
// https://github.com/abseil/abseil-cpp/blob/master/absl/container/flat_hash_set.h
#include "absl/container/flat_hash_set.h"
#include "ortools/base/logging.h"
#include "ortools/base/map_util.h"

namespace operations_research {

const int64_t RoutingIndexManager::kUnassigned = -1;

//TH3: Xác định số lượng tập điểm di chuyển(num_nodes)
//     số lượng xe cần tìm hành trình di chuyển(num_vehicles)
//     với tất cả các xe thứ i ta sẽ hiểu cần tìm hành trình đi từ 
//     start[i] = depot -> end[i] = depot
//     -> gọi lại hàm hàm khởi tạo 
//RoutingIndexManager(num_nodes,num_vehicles,[
//  {start[1]=depot           ,end[1]=depot},
//  {start[2]=depot           ,end[2]=depot},
//  ...
//  {start[num_vehicles]=depot,end[num_vehicles]=depot}
//])
RoutingIndexManager::RoutingIndexManager(int num_nodes, int num_vehicles,
                                         NodeIndex depot)
    : RoutingIndexManager(num_nodes, num_vehicles,
                          std::vector<std::pair<NodeIndex, NodeIndex>>(
                              num_vehicles, {depot, depot})) {}


//TH1: Xác định số lượng tập điểm di chuyển(num_nodes),
//     số lượng xe cần tìm hành trình di chuyển(num_vehicles)
//     tương ứng với xe thứ i ta có start[i] và end[i] ta cần tìm hành trình
//     thoả mãn cho xe thứ i đi từ điểm start[i] tới end[i]
RoutingIndexManager::RoutingIndexManager(int num_nodes, int num_vehicles,
                                         const std::vector<NodeIndex>& starts,
                                         const std::vector<NodeIndex>& ends) {
  CHECK_EQ(starts.size(), num_vehicles);
  CHECK_EQ(ends.size(), num_vehicles);
  std::vector<std::pair<NodeIndex, NodeIndex>> starts_ends(num_vehicles);
  for (int v = 0; v < num_vehicles; ++v) {
    starts_ends[v] = {starts[v], ends[v]};
  }
  Initialize(num_nodes, num_vehicles, starts_ends);
}

// tương tự trường hợp 1  nhưng đưa tập diểm {start,end}
//     là pair  
RoutingIndexManager::RoutingIndexManager(
    int num_nodes, int num_vehicles,
    const std::vector<std::pair<NodeIndex, NodeIndex>>& starts_ends) {
  Initialize(num_nodes, num_vehicles, starts_ends);
}

// Hàm khởi tạo chung 
void RoutingIndexManager::Initialize(
    int num_nodes, int num_vehicles,
    const std::vector<std::pair<NodeIndex, NodeIndex>>& starts_ends) {
  static const NodeIndex kZeroNode(0);

  num_nodes_ = num_nodes;
  num_vehicles_ = num_vehicles;
  CHECK_EQ(num_vehicles_, starts_ends.size());
  absl::flat_hash_set<NodeIndex> starts;
  absl::flat_hash_set<NodeIndex> ends;
  absl::flat_hash_set<NodeIndex> unique_depots;
  
  for (const std::pair<NodeIndex, NodeIndex>& start_end : starts_ends) {
    const NodeIndex start = start_end.first;
    const NodeIndex end = start_end.second;

    // ? có thể là kiểm tra start < 0 || start > num_nodes || end < 0 || end < num_nodes
    CHECK_GE(start, 0);
    CHECK_GE(end, 0);
    CHECK_LE(start, num_nodes_);
    CHECK_LE(end, num_nodes_);

    // Thêm vào tập quản lý các điểm bắt đầu và kết thúc giống với Set của C++ thông thường   
    starts.insert(start);
    ends.insert(end);

    // Kiểm tra xem có phải unique depot không khi unique_depot = 1
    unique_depots.insert(start);
    unique_depots.insert(end);
  }

  num_unique_depots_ = unique_depots.size();

  const int size = num_nodes_ + num_vehicles_ - num_unique_depots_;


  index_to_node_.resize(size + num_vehicles_);
  
  node_to_index_.resize(num_nodes_, kUnassigned);
  
  vehicle_to_start_.resize(num_vehicles_);
  
  vehicle_to_end_.resize(num_vehicles_);
  
  int64_t index = 0;
  

  // mapping struct {int} to int
  for (NodeIndex i = kZeroNode; i < num_nodes_; ++i) {
    if (starts.contains(i) || !ends.contains(i)) {
      index_to_node_[index] = i;
      node_to_index_[i] = index;
      ++index;
    }
  }

  // Nếu có start trùng nhau thì tạo ra nút start ảo và không ánh xạ lại
  // -> nếu là node thật thì index_to_node[] và node_to_index[] cùng tồn tại
  // -> nếu là node out thì chỉ tồn tại index_to_node[]
  // tương tự như vậy với end
  // ind
  
  absl::flat_hash_set<NodeIndex> seen_starts;
  for (int i = 0; i < num_vehicles_; ++i) {
    const NodeIndex start = starts_ends[i].first;
    if (!seen_starts.contains(start)) {
      seen_starts.insert(start);
      const int64_t start_index = node_to_index_[start];
      vehicle_to_start_[i] = start_index;
      CHECK_NE(kUnassigned, start_index);
    } else {
      vehicle_to_start_[i] = index;
      index_to_node_[index] = start;
      ++index;
    }
  }
  for (int i = 0; i < num_vehicles_; ++i) {
    NodeIndex end = starts_ends[i].second;
    index_to_node_[index] = end;
    vehicle_to_end_[i] = index;
    CHECK_LE(size, index);
    ++index;
  }

  // Logging model information.
  VLOG(1) << "Number of nodes: " << num_nodes_;
  VLOG(1) << "Number of vehicles: " << num_vehicles_;
  for (int64_t index = 0; index < index_to_node_.size(); ++index) {
    VLOG(2) << "Variable index " << index << " -> Node index "
            << index_to_node_[index];
  }
  for (NodeIndex node = kZeroNode; node < node_to_index_.size(); ++node) {
    VLOG(2) << "Node index " << node << " -> Variable index "
            << node_to_index_[node];
  }
}

std::vector<int64_t> RoutingIndexManager::NodesToIndices(
    const std::vector<NodeIndex>& nodes) const {
  std::vector<int64_t> indices;
  indices.reserve(nodes.size());
  for (const NodeIndex node : nodes) {
    const int64_t index = NodeToIndex(node);
    CHECK_NE(kUnassigned, index);
    indices.push_back(index);
  }
  return indices;
}

std::vector<RoutingIndexManager::NodeIndex> RoutingIndexManager::IndicesToNodes(
    const std::vector<int64_t>& indices) const {
  std::vector<NodeIndex> nodes;
  nodes.reserve(indices.size());
  for (const int64_t index : indices) {
    nodes.push_back(IndexToNode(index));
  }
  return nodes;
}

}  // namespace operations_research
