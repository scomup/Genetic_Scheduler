#include <iostream>
#include <assert.h>
#include <omp.h>

#include "GeneticSchedulerCore.h"

namespace Scheduler
{

GeneticSchedulerCore::GeneticSchedulerCore(std::vector<Node> nodes,
                                           int16_t all_core_num, uint32_t max_cell_num = 10000,
                                           uint16_t death_rate = 10,
                                           float mutation_rate = 0.1,
                                           uint16_t max_loop = 100)
    : nodes_(nodes)
    , all_core_num_(all_core_num)
    , best_result_(std::numeric_limits<int16_t>::max())
    , max_cell_num_(max_cell_num)
    , death_rate_(death_rate)
    , mutation_rate_(mutation_rate)
    , max_loop_(max_loop)
{

    auto p0 = std::chrono::system_clock::now();
    std::vector<CellWithScore> cellWithScores;

    for (int i = 0, N = omp_get_max_threads(); i < N; ++i)
    {
        seeds_.emplace_back(time(NULL));
    }

    //#pragma omp parallel for
    for (uint32_t i = 0; i < max_cell_num_; i++)
    {
        auto cell = generate_new_cell();
        int16_t score = evaluate(cell);
        CellWithScore cell_with_score{std::move(cell), score};
        //#pragma omp critical
        {
            cellWithScores.emplace_back(std::move(cell_with_score));
        }
    }

    std::sort(cellWithScores.begin(), cellWithScores.end(), [](const CellWithScore &a, const CellWithScore &b) { return a.score < b.score; });

    auto p1 = std::chrono::system_clock::now();
    auto diff1 = p1 - p0;

    for (uint32_t loop = 0; loop < max_loop_; loop++)
    {
        auto p2 = std::chrono::system_clock::now();

        cellWithScores.erase(cellWithScores.begin() + (max_cell_num_ / death_rate), cellWithScores.end());
        std::vector<CellWithScore> new_cellWithScores;

        #pragma omp parallel for
        for (size_t i = 0; i < cellWithScores.size(); i++)
        {
            for (size_t j = 0; j < death_rate; j++)
            {
                auto cell = create_next_generation(cellWithScores[i].cell);
                int16_t score = evaluate(cell);
                CellWithScore cell_with_score{std::move(cell), score};
                #pragma omp critical
                {
                    new_cellWithScores.emplace_back(std::move(cell_with_score));
                }
            }
        }

        new_cellWithScores.swap(cellWithScores);
        std::sort(cellWithScores.begin(), cellWithScores.end(), [](const CellWithScore &a, const CellWithScore &b) { return a.score < b.score; });
        if(best_result_ > cellWithScores[0].score){
            best_result_ = cellWithScores[0].score;
        }

        auto p3 = std::chrono::system_clock::now();
        auto diff2 = p3 - p2;
        printf("time:%5ld     best:%4d   worst%4d   survivor:%4d~%4d\n",
               std::chrono::duration_cast<std::chrono::milliseconds>(diff2).count(),
               cellWithScores[0].score,
               cellWithScores.back().score,
               cellWithScores[0].score,
               cellWithScores[max_cell_num_ / death_rate].score);
    }


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
        int16_t index = Scheduler::common::randi<int16_t>(min_index, i, seeds_[omp_get_thread_num()]);
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
        float r = static_cast <float> (rand_r(&seeds_[omp_get_thread_num()])) / static_cast <float> (RAND_MAX);
        
        if (r < mutation_rate_ )
        {
            int16_t min_index = 0;
            for (int16_t j : nodes_[i].sub_nodes)
            {
                min_index = std::max((*new_cell)[j], min_index);
            }
            min_index++;
            int16_t current_index = (*new_cell)[i];
            assert(min_index <= current_index);
            int16_t index = Scheduler::common::randi<int16_t>(min_index, current_index, seeds_[omp_get_thread_num()]);
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
