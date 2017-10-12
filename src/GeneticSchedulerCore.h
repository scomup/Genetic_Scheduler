#ifndef SCHEDULE_STATE
#define SCHEDULE_STATE

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
    GeneticSchedulerCore(std::vector<Node> nodes, int16_t all_core_num, uint32_t max_cell_num, uint16_t death_rate, float mutation_rate);
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
};

GeneticSchedulerCore::GeneticSchedulerCore(std::vector<Node> nodes, int16_t all_core_num, uint32_t max_cell_num = 10000, uint16_t death_rate = 10, float mutation_rate = 0.1) 
:nodes_(nodes)
,all_core_num_(all_core_num)
,best_result_(std::numeric_limits<int16_t>::max())
,max_cell_num_(max_cell_num)
,death_rate_(death_rate)
,mutation_rate_(mutation_rate)
{
    std::vector<CellWithScore> cellWithScores;

    auto p0 = std::chrono::system_clock::now();
    for(uint32_t i = 0; i < max_cell_num_; i++){
        auto cell = generate_new_cell();
        int16_t score = evaluate(cell);
        CellWithScore cell_with_score{std::move(cell), score};
        cellWithScores.emplace_back(std::move(cell_with_score));
    }
    auto p1= std::chrono::system_clock::now();
    auto diff1 = p1 - p0;
    std::cout << "generate_new_cell:" << std::chrono::duration_cast<std::chrono::seconds>(diff1).count() <<" sec." <<std::endl;

    while (true)
    {
        std::sort(cellWithScores.begin(), cellWithScores.end(), [](const CellWithScore &a, const CellWithScore &b) { return a.score < b.score; });
        
        std::cout << "best:" << cellWithScores[0].score <<std::endl;

        cellWithScores.erase(cellWithScores.begin() + (max_cell_num_ / death_rate), cellWithScores.end());
        std::vector<CellWithScore> new_cellWithScores;
        for (size_t i = 0; i < cellWithScores.size(); i++)
        {
            for (size_t j = 0; j < death_rate; j++)
            {
                auto cell = create_next_generation(cellWithScores[0].cell);
                int16_t score = evaluate(cell);
                CellWithScore cell_with_score{std::move(cell), score};
                new_cellWithScores.emplace_back(std::move(cell_with_score));
            }
        }
        new_cellWithScores.swap(cellWithScores);
    }

    /*
    auto p1= std::chrono::system_clock::now();
    auto diff1 = p1 - p0;
    std::cout << "generate_new_cell:" << std::chrono::duration_cast<std::chrono::seconds>(diff1).count() <<" sec." <<std::endl;
  
    for(uint32_t i = 0; i < max_cell_num_; i++){
        int16_t evaluate_result = evaluate(cells[i]);
        if (evaluate_result < best_result_)
        {
            best_result_ = evaluate_result;
        }
    }
    auto p2= std::chrono::system_clock::now();
    auto diff2 = p2 - p1;
    std::cout << "evaluate_results:" << std::chrono::duration_cast<std::chrono::seconds>(diff2).count() <<" sec." <<std::endl;
    std::cout << "best_results:" << best_result_ <<std::endl;
    */
    /*
    uint64_t c = 0;
    while (true)
    {
        auto cell = generate_new_cell();

        int16_t evaluate_result = evaluate(std::move(cell));
        if (evaluate_result < best_result_)
        {
            best_result_ = evaluate_result;
            std::cout << "new:" << best_result_ << std::endl;
        }
        if (c % 100000 == 0)
        {
            std::cout << c << ":" << best_result_ << std::endl;
        }
        c++;
    }*/
};

std::unique_ptr<std::vector<int16_t>> GeneticSchedulerCore::generate_new_cell()
{
    int16_t node_num = nodes_.size() - 1;
    auto cell = common::make_unique<std::vector<int16_t>>(node_num);
    (*cell)[0] = 0;
    for (int16_t i = 1; i < node_num; i++)
    {

        int16_t min_index = 0;
        for (int16_t j : nodes_[i].sub_nodes)
        {
            min_index = std::max((*cell)[j], min_index);
        }
        min_index++;
        int16_t index = Scheduler::common::randi<int16_t>(min_index, i);
        (*cell)[i] = index;
        for (int16_t j = 1; j < i; j++)
        {
            if ((*cell)[j] >= index)
            {
                (*cell)[j]++;
            }
        }
    }
    return cell;
}

int16_t GeneticSchedulerCore::evaluate(std::unique_ptr<std::vector<int16_t>>& cell)
{
    std::vector<int16_t> idx = Scheduler::common::sort_indexes<int16_t>(*cell);
    std::vector<int16_t> nodes_finish_time(nodes_.size());
    std::vector<int16_t> cores_ocuppied_time(all_core_num_);

    std::fill(nodes_finish_time.begin(), nodes_finish_time.end(), inf);
    std::fill(cores_ocuppied_time.begin(), cores_ocuppied_time.end(), 0);

    for (int16_t i : idx)
    {
        int16_t schedulabe_time = 0;
        for (size_t j = 0; j < nodes_[i].sub_nodes.size(); j++)
        {
            int16_t sub_node_finished_time = nodes_finish_time[nodes_[i].sub_nodes[j]];
            assert(sub_node_finished_time != inf);
            schedulabe_time = std::max(sub_node_finished_time, schedulabe_time);
        }

        int16_t core_num = nodes_[i].core_num;
        std::sort(cores_ocuppied_time.begin(), cores_ocuppied_time.end());
        schedulabe_time = std::max(schedulabe_time, cores_ocuppied_time[core_num - 1]);
        nodes_finish_time[i] = schedulabe_time + nodes_[i].time;
        for (int16_t j = 0; j < nodes_[i].core_num; j++)
        {
            cores_ocuppied_time[j] = nodes_finish_time[i];
        }
    }
    std::sort(cores_ocuppied_time.begin(), cores_ocuppied_time.end());
    return cores_ocuppied_time.back();
}

std::unique_ptr<std::vector<int16_t>> GeneticSchedulerCore::create_next_generation(std::unique_ptr<std::vector<int16_t>>& cell)
{
    auto new_cell = common::make_unique<std::vector<int16_t>>(*cell);
    int16_t node_num = nodes_.size() - 1;

    for (int16_t i = 1; i < node_num; i++)
    {
        float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (r < mutation_rate_)
        {
            int16_t min_index = 0;
            for (int16_t j : nodes_[i].sub_nodes)
            {
                min_index = std::max((*new_cell)[j], min_index);
            }
            min_index++;
            int16_t current_index = (*new_cell)[i];
            int16_t index = Scheduler::common::randi<int16_t>(min_index, current_index);
            for (int16_t j = 1; j < node_num; j++)
            {
                if ((*new_cell)[j] >= index && (*new_cell)[j] < current_index)
                {
                    (*new_cell)[j]++;
                }
            }
            (*new_cell)[i] = index;
        }
    }
    return new_cell;
}
}

#endif
