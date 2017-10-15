#ifndef NODE_H
#define NODE_H

#include <memory>
#include <stdint.h>

namespace Scheduler
{

struct Node
{
    int16_t iD;
    int16_t core_num;
    int16_t time;
    std::vector<int16_t> sub_nodes;
};

struct chromosomeWithScore
{
    std::unique_ptr<std::vector<int16_t>> chromosome;
    int16_t score;
};

}

#endif