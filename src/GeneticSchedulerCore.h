#ifndef GENETICS_CHEDULER_CORE_H
#define GENETICS_CHEDULER_CORE_H

#include <algorithm>
#include <limits>
#include <set>
#include <unordered_set>
#include <chrono>

#include "common/make_unique.h"
#include "common/common.h"

#include "node.h"

namespace Scheduler
{

constexpr int16_t inf = -1;

struct CellWithScore
{
    std::unique_ptr<std::vector<int16_t>> cell;
    int16_t score;
};

class GeneticSchedulerCore
{
  public:
    GeneticSchedulerCore(std::vector<Node> nodes, int16_t all_core_num, uint32_t max_cell_num, uint16_t death_rate, float mutation_rate, uint16_t max_loop);
    int16_t getBest(){return best_result_;};

  private:
    std::unique_ptr<std::vector<int16_t>> generate_new_cell();
    int16_t evaluate(std::unique_ptr<std::vector<int16_t>>& cell);
    std::unique_ptr<std::vector<int16_t>>  create_next_generation(std::unique_ptr<std::vector<int16_t>>& cell);

  private: 
    std::vector<Node> nodes_;
    int16_t all_core_num_;
    int16_t best_result_;
    uint32_t max_cell_num_;
    uint16_t death_rate_;
    float mutation_rate_;
    uint16_t max_loop_;
    std::vector<uint32_t> seeds_;
};

}

#endif
