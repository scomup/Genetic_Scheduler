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
    //Initial seed(for rand) for each openmp thread.
    //this step is very essential because multiple threads access one seed may leads to cache issues.
    auto p0 = std::chrono::system_clock::now();
    std::vector<chromosomeWithScore> chromosomeWithScores;
    for (int i = 0, N = omp_get_max_threads(); i < N; ++i)
    {
        seeds_.emplace_back(time(NULL)^(i+1));
    }

    //Generate the initial chromosomes and evaluate them.
    #pragma omp parallel for
    for (uint32_t i = 0; i < config_ptr_->max_chromosome_num; i++)
    {
        auto chromosome = generate_new_chromosome();
        int16_t score = evaluate(chromosome);
        chromosomeWithScore chromosome_with_score{std::move(chromosome), score};
        #pragma omp critical
        {
            chromosomeWithScores.emplace_back(std::move(chromosome_with_score));
        }
    }

    //To easy find better chromosomes, we just sort them by score.
    //If use roulette may not necessary...  
    std::sort(chromosomeWithScores.begin(), chromosomeWithScores.end(), [](const chromosomeWithScore &a, const chromosomeWithScore &b) { return a.score < b.score; });

    auto p1 = std::chrono::system_clock::now();
    auto diff1 = p1 - p0;

    for (uint32_t loop = 0; loop < config_ptr_->max_loop; loop++)
    {
        auto p2 = std::chrono::system_clock::now();

        std::vector<chromosomeWithScore> new_chromosomeWithScores(config_ptr_->max_chromosome_num);

        //Create a roulette by the score of chromosomes 
        Roulette roulette = Roulette(chromosomeWithScores, config_ptr_->aphla);

        //char buffer [50];
        //sprintf (buffer, "%d.txt", loop);
        //std::ofstream ofs (buffer, std::ofstream::out);
        //std::vector<int16_t> nnn(config_ptr_->max_chromosome_num);
        //std::fill(nnn.begin(),nnn.end(),0);
        #pragma omp parallel for
        for (size_t i = 0; i < config_ptr_->max_chromosome_num; i++)
        {
            
            //takes two chromosomes (parents) by roulette and swaps part of 
            //their  genetic  information  to  create  new chromosome
            uint32_t select_a;
            uint32_t select_b;
            do{
                select_a = roulette.spin_roulette(seeds_[omp_get_thread_num()] );
                select_b = roulette.spin_roulette(seeds_[omp_get_thread_num()] );
            }
            while(select_a == select_b);
            auto chromosome = crossover(chromosomeWithScores[select_a].chromosome,
                                        chromosomeWithScores[select_b].chromosome,
                                        seeds_[omp_get_thread_num()]);
            //change the values of one or more genes in the chromosome by probability
            mutation(chromosome);

            //nnn[select] ++;
            //Evaluate the chromosome.
            int16_t score = evaluate(chromosome);
            chromosomeWithScore chromosome_with_score{std::move(chromosome), score};
            new_chromosomeWithScores[i] = std::move(chromosome_with_score);
        }

        /*
        for (size_t i = 0; i < config_ptr_->max_chromosome_num; i++)
        {
            if(i>=1){
                ofs<<roulette[i]-roulette[i-1]<<",";
            }
            else{
                ofs<<roulette[i]<<",";
            }
            ofs <<chromosomeWithScores[i].score<<","<< nnn[i] <<std::endl;
        }
        ofs.close();
        */

        new_chromosomeWithScores.swap(chromosomeWithScores);
        std::sort(chromosomeWithScores.begin(), chromosomeWithScores.end(), [](const chromosomeWithScore &a, const chromosomeWithScore &b) { return a.score < b.score; });
        if (best_result_ > chromosomeWithScores[0].score)
        {
            best_result_ = chromosomeWithScores[0].score;
        }


    auto p3 = std::chrono::system_clock::now();
    auto diff2 = p3 - p2;

    
    //printf("time:%5ld\n",
    //       std::chrono::duration_cast<std::chrono::milliseconds>(diff2).count());
    //return;

    printf("time:%5ld     best:%4d   worst%4d\n",
           std::chrono::duration_cast<std::chrono::milliseconds>(diff2).count(),
           chromosomeWithScores[0].score,
           chromosomeWithScores.back().score);
    //return;
    }


};

//-----------------------------------------------------------------------------
//This function used to generate new chromosomes.
//A valid chromosome is an int array of length N (N = number of nodes ).
//denote as C[N]
//C[i] indicate the scheduling order of node_i. 
//if node_i is a parent node of node_j, C[i] < C[j]
//-----------------------------------------------------------------------------
std::unique_ptr<std::vector<int16_t>> GeneticSchedulerCore::generate_new_chromosome()
{
    int16_t node_num = nodes_.size() - 1;
    auto chromosome = common::make_unique<std::vector<int16_t>>(node_num);
    (*chromosome)[0] = 0;
    for (int16_t i = 1; i < node_num; i++)
    {
        int16_t min_index = 0;
        for (int16_t j : nodes_[i].sub_nodes)
        {
            min_index = std::max((*chromosome)[j], min_index);
        }
        min_index++;
        int16_t index = Scheduler::common::randi<int16_t>(min_index, i, seeds_[omp_get_thread_num()]);
        (*chromosome)[i] = index;
        for (int16_t j = 1; j < i; j++)
        {
            if ((*chromosome)[j] >= index)
            {
                (*chromosome)[j]++;
            }
        }
    }
    return chromosome;
}

//-----------------------------------------------------------------------------
//This funciton try to assembly of the scheduling result according to the genetic information.

//We choose node by scheduling order carried by chromosome in each step.
//We calcuate the finished time for this node_i.
//  .time_a = Max(the finished time of all parent task of node_i)
//  .time_b = The most recent time when system can provide enough core resources for this node.
//  .The finished time of this node = Max(time_a, time_a) + the runing time of this node.
//Update the occupied time of selected cores.
//Finished if all node be processed.
//Return the longest occupied time of cores.
//-----------------------------------------------------------------------------

int16_t GeneticSchedulerCore::evaluate(std::unique_ptr<std::vector<int16_t>> &chromosome)
{
    std::vector<int16_t> idx = Scheduler::common::sort_indexes<int16_t>(*chromosome);
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

//-----------------------------------------------------------------------------
//This funciton try to create new generation with mutation.
//-----------------------------------------------------------------------------

void GeneticSchedulerCore::mutation(std::unique_ptr<std::vector<int16_t>>& chromosome)
{
    auto new_chromosome = common::make_unique<std::vector<int16_t>>(*chromosome);
    int16_t node_num = nodes_.size() - 1;

    for (int16_t i = 1; i < node_num; i++)
    {
        float r = static_cast <float> (rand_r(&seeds_[omp_get_thread_num()])) / static_cast <float> (RAND_MAX);
        
        if (r < config_ptr_->mutation_rate )
        {
            int16_t min_index = 0;
            for (int16_t j : nodes_[i].sub_nodes)
            {
                min_index = std::max((*chromosome)[j], min_index);
            }
            min_index++;
            int16_t current_index = (*chromosome)[i];
            assert(min_index <= current_index);
            int16_t index = Scheduler::common::randi<int16_t>(min_index, current_index, seeds_[omp_get_thread_num()]);
            for (int16_t j = 1; j < node_num; j++)
            {
                if ((*chromosome)[j] >= index && (*chromosome)[j] < current_index)
                {
                    (*chromosome)[j]++;
                }
            }
            (*chromosome)[i] = index;
        }
    }
}

std::unique_ptr<std::vector<int16_t>> GeneticSchedulerCore::crossover(std::unique_ptr<std::vector<int16_t>> &chromosome_a,
                                                                      std::unique_ptr<std::vector<int16_t>> &chromosome_b,
                                                                      uint32_t &seed)
{
    std::vector<int16_t> idx_a = Scheduler::common::sort_indexes<int16_t>(*chromosome_a);
    std::vector<int16_t> idx_b = Scheduler::common::sort_indexes<int16_t>(*chromosome_b);
    int16_t node_num = (*chromosome_a).size();
    std::vector<bool> mark(node_num);
    //std::fill(mark.begin(),mark.end(),false);
    int16_t cutting_point = Scheduler::common::randi<int16_t>(1, node_num-1, seed);

    for(int16_t i = 0; i<cutting_point; i++)
    {
        mark[idx_a[i]] = true;
    }
    for(int16_t i = 0; i<node_num; i++)
    {
        if(mark[idx_b[i]]){
            continue;
        }
        idx_a[cutting_point++] = idx_b[i];
    }

    std::vector<int16_t> chromosome = Scheduler::common::sort_indexes<int16_t>(idx_a);
    auto chromosome_ptr = common::make_unique<std::vector<int16_t>>(chromosome);
    //evaluate(chromosome_ptr);
    return chromosome_ptr;
}
}
