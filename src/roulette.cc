#include "roulette.h"

#include <algorithm>
#include <assert.h>
#include <omp.h>

namespace Scheduler
{
Roulette::Roulette(std::vector<chromosomeWithScore> &chromosomeWithScores, double aphla)
{
    roulette_.resize(chromosomeWithScores.size());
    int16_t best_result = chromosomeWithScores[0].score;
    double sum = 0;
    #pragma omp parallel for reduction(+ : sum)
    for (size_t i = 0; i < roulette_.size(); i++)
    {
        int16_t current_result = chromosomeWithScores[i].score;
        int16_t diff = current_result - best_result;
        assert(diff >= 0);
        roulette_[i] = exp(-(aphla * diff));
        sum += roulette_[i];
    }
    assert(sum != 0);
    roulette_[0] = roulette_[0] / sum;
    //#pragma omp parallel for
    for (size_t i = 1; i < roulette_.size(); i++)
    {

        roulette_[i] = roulette_[i - 1] + roulette_[i] / sum;
        //std::cout<<roulette[i]<<std::endl;
    }
    //for(size_t i = 1; i < roulette.size(); i++){
    //
    //    std::cout<<roulette[i]<<std::endl;
    //}
}

int32_t Roulette::spin_roulette(uint32_t &seed)
{
    double r = static_cast<double>(rand_r(&seed)) / static_cast<double>(RAND_MAX);
    int32_t index = roulette_.size();
    int32_t step = index;
    while (true)
    {
        if (step > 1)
            step /= 2;

        if (index == 0 || (r >= roulette_[index - 1] && r < roulette_[index]))
        {
            break;
        }
        if (r >= roulette_[index - 1])
        {
            index += step;
        }
        else
        {
            index -= step;
        }
    }
    return index;
}

}
