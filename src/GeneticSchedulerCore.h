#ifndef GENETICS_CHEDULER_CORE_H
#define GENETICS_CHEDULER_CORE_H

#include <algorithm>
#include <limits>
#include <chrono>

#include "common/make_unique.h"
#include "common/common.h"
#include "config.h"
#include "node.h"
#include "roulette.h"

namespace Scheduler
{

constexpr int16_t inf = -1;

class GeneticSchedulerCore
{
public:
  GeneticSchedulerCore(std::vector<Node> nodes, Scheduler::Config *config_ptr);
  int16_t getBest() { return best_result_; };

private:
  std::unique_ptr<std::vector<int16_t>> generate_new_chromosome();
  int16_t evaluate(std::unique_ptr<std::vector<int16_t>> &chromosome);
  void mutation(std::unique_ptr<std::vector<int16_t>> &chromosome);
  std::unique_ptr<std::vector<int16_t>> crossover(std::unique_ptr<std::vector<int16_t>> &chromosome_a,
                                                  std::unique_ptr<std::vector<int16_t>> &chromosome_b,
                                                  uint32_t &seed);

private:
  std::vector<Node> nodes_;
  Scheduler::Config *config_ptr_;
  int16_t best_result_;
  std::vector<uint32_t> seeds_;
};
}

#endif
