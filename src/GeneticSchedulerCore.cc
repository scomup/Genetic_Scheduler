#include <iostream>
#include <assert.h>
#include <omp.h>
#include <fstream>

#include "GeneticSchedulerCore.h"

namespace Scheduler
{

GeneticSchedulerCore::GeneticSchedulerCore(std::vector<Node> nodes, Scheduler::Config* config_ptr)
    : nodes_(nodes)
    , config_ptr_(config_ptr)
    , best_result_(std::numeric_limits<int16_t>::max())
{

    auto p0 = std::chrono::system_clock::now();
    std::vector<CellWithScore> cellWithScores;

    for (int i = 0, N = omp_get_max_threads(); i < N; ++i)
    {
        seeds_.emplace_back(time(NULL)^(i+1));
    }

    #pragma omp parallel for
    for (uint32_t i = 0; i < config_ptr_->max_cell_num; i++)
    {
        auto cell = generate_new_cell();
        int16_t score = evaluate(cell);
        CellWithScore cell_with_score{std::move(cell), score};
        #pragma omp critical
        {
            cellWithScores.emplace_back(std::move(cell_with_score));
        }
    }

    std::sort(cellWithScores.begin(), cellWithScores.end(), [](const CellWithScore &a, const CellWithScore &b) { return a.score < b.score; });

    auto p1 = std::chrono::system_clock::now();
    auto diff1 = p1 - p0;

    for (uint32_t loop = 0; loop < config_ptr_->max_loop; loop++)
    {
        auto p2 = std::chrono::system_clock::now();

        std::vector<CellWithScore> new_cellWithScores(config_ptr_->max_cell_num);

        auto roulette = create_roulette(cellWithScores);

        //char buffer [50];
        //sprintf (buffer, "%d.txt", loop);
        //std::ofstream ofs (buffer, std::ofstream::out);
        //std::vector<int16_t> nnn(config_ptr_->max_cell_num);
        //std::fill(nnn.begin(),nnn.end(),0);
        #pragma omp parallel for
        for (size_t i = 0; i < config_ptr_->max_cell_num; i++)
        {
            uint32_t select = spin_roulette(roulette,seeds_[omp_get_thread_num()] );

            //nnn[select] ++;
            
            //std::cout<<index<<std::endl;
            //double r = std::abs(Scheduler::common::rand_normal<double>(0, 1./death_rate, seeds_[omp_get_thread_num()]));
            //if(r > 1)
            //r = fmod (r,1.);
            //uint32_t select = r*config_ptr_->max_cell_num;
            auto cell = create_next_generation(cellWithScores[select].cell);
            int16_t score = evaluate(cell);
            CellWithScore cell_with_score{std::move(cell), score};
            new_cellWithScores[i] = std::move(cell_with_score);
        }

        /*
        for (size_t i = 0; i < config_ptr_->max_cell_num; i++)
        {
            if(i>=1){
                ofs<<roulette[i]-roulette[i-1]<<",";
            }
            else{
                ofs<<roulette[i]<<",";
            }
            ofs <<cellWithScores[i].score<<","<< nnn[i] <<std::endl;
        }
        ofs.close();
        */

        new_cellWithScores.swap(cellWithScores);
        std::sort(cellWithScores.begin(), cellWithScores.end(), [](const CellWithScore &a, const CellWithScore &b) { return a.score < b.score; });
        if (best_result_ > cellWithScores[0].score)
        {
            best_result_ = cellWithScores[0].score;
        }


    auto p3 = std::chrono::system_clock::now();
    auto diff2 = p3 - p2;
    //printf("time:%5ld\n",
    //       std::chrono::duration_cast<std::chrono::milliseconds>(diff2).count());

    printf("time:%5ld     best:%4d   worst%4d\n",
           std::chrono::duration_cast<std::chrono::milliseconds>(diff2).count(),
           cellWithScores[0].score,
           cellWithScores.back().score);
    }


};

int32_t GeneticSchedulerCore::spin_roulette(std::vector<double> &roulette, uint32_t &seed)
{
    double r = static_cast<double>(rand_r(&seed)) / static_cast<double>(RAND_MAX);
    int32_t index = config_ptr_->max_cell_num / 2;
    int32_t step = index;
    while (true)
    {
        if (step > 1)
            step /= 2;

        if (index == 0 || (r >= roulette[index - 1] && r < roulette[index]))
        {
            break;
        }
        if (r >= roulette[index - 1])
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

std::vector<double> GeneticSchedulerCore::create_roulette(std::vector<CellWithScore> &cellWithScores)
{
    std::vector<double> roulette = std::vector<double>(cellWithScores.size());
    int16_t best_result = cellWithScores[0].score;
    double sum = 0;
    #pragma omp parallel for reduction(+ : sum)
    for (size_t i = 0; i < roulette.size(); i++)
    {
        int16_t current_result = cellWithScores[i].score;
        int16_t diff = current_result - best_result;
        assert(diff >= 0);
        roulette[i] = exp(-(config_ptr_->aphla * diff));
        sum += roulette[i];
    }
    assert(sum != 0);
    roulette[0] = roulette[0] / sum;
    //#pragma omp parallel for
    for (size_t i = 1; i < roulette.size(); i++)
    {

        roulette[i] = roulette[i - 1] + roulette[i] / sum;
        //std::cout<<roulette[i]<<std::endl;
    }
    //for(size_t i = 1; i < roulette.size(); i++){
    //
    //    std::cout<<roulette[i]<<std::endl;
    //}
    return roulette;
}

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

int16_t GeneticSchedulerCore::evaluate(std::unique_ptr<std::vector<int16_t>> &cell)
{
    std::vector<int16_t> idx = Scheduler::common::sort_indexes<int16_t>(*cell);
    std::vector<int16_t> nodes_finish_time(nodes_.size());
    std::vector<int16_t> cores_ocuppied_time(config_ptr_->all_core_num);

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
        
        if (r < config_ptr_->mutation_rate )
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
