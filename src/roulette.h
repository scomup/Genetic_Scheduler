#ifndef ROULETTE_H
#define ROULETTE_H

#include <vector>

#include "node.h"

namespace Scheduler
{

class Roulette
{
  public:
    Roulette(std::vector<chromosomeWithScore> &chromosomeWithScores, double aphla);
    int32_t spin_roulette(uint32_t &seed);

  private:
    std::vector<double> roulette_;

};

}

#endif