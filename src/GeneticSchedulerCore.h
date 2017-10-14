#ifndef GENETICS_CHEDULER_CORE_H
#define GENETICS_CHEDULER_CORE_H

#include <algorithm>
#include <limits>
#include <set>
#include <unordered_set>
#include <chrono>

#include "common/make_unique.h"
#include "common/common.h"
#include "config.h"
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
    GeneticSchedulerCore(std::vector<Node> nodes, Scheduler::Config* config_ptr);
    int16_t getBest(){return best_result_;};

  private:
    std::unique_ptr<std::vector<int16_t>> generate_new_cell();
    int16_t evaluate(std::unique_ptr<std::vector<int16_t>>& cell);
    std::unique_ptr<std::vector<int16_t>>  create_next_generation(std::unique_ptr<std::vector<int16_t>>& cell);
    std::vector<double> create_roulette(std::vector<CellWithScore>& cellWithScores);  
    int32_t spin_roulette(std::vector<double>& roulette, uint32_t& seed);
      


  private: 
    std::vector<Node> nodes_;
    Scheduler::Config* config_ptr_;
    int16_t best_result_;
    std::vector<uint32_t> seeds_;
};

}

#endif
