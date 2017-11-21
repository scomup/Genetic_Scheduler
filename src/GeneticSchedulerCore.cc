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
    std::vector<chromosomeWithScore> chromosomeWithScores;
    for (int i = 0, N = omp_get_max_threads(); i < N; ++i)
    {
        seeds_.emplace_back(time(NULL)^(i+1));
    }

};

void GeneticSchedulerCore::run()
{
    auto p0 = std::chrono::system_clock::now();
    std::vector<chromosomeWithScore> chromosomeWithScores;
    for (int i = 0, N = omp_get_max_threads(); i < N; ++i)
    {
        seeds_.emplace_back(time(NULL) ^ (i + 1));
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

        //Create a roulette by the score of all chromosomes
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
            do
            {
                select_a = roulette.spin_roulette(seeds_[omp_get_thread_num()]);
                select_b = roulette.spin_roulette(seeds_[omp_get_thread_num()]);
            } while (select_a == select_b);
            auto chromosome_ptr = common::make_unique<std::vector<int16_t>>(*chromosomeWithScores[select_a].chromosome);
            if (config_ptr_->use_crossover)
            {
                chromosome_ptr = crossover(chromosomeWithScores[select_a].chromosome,
                                           chromosomeWithScores[select_b].chromosome,
                                           seeds_[omp_get_thread_num()]);
            }

            //change the values of one or more genes in the chromosome by probability
            mutation(chromosome_ptr);

            //nnn[select] ++;
            //Evaluate the chromosome.
            int16_t score = evaluate(chromosome_ptr);
            chromosomeWithScore chromosome_with_score{std::move(chromosome_ptr), score};
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
        if (config_ptr_->show_step_info)
        {
            printf("time:%5ld     best:%4d   worst%4d\n",
                   std::chrono::duration_cast<std::chrono::milliseconds>(diff2).count(),
                   chromosomeWithScores[0].score,
                   chromosomeWithScores.back().score);
        }
    }
}
//-----------------------------------------------------------------------------
//This function used to generate new chromosomes.
//A valid chromosome is an int array of length N (N = number of nodes ).
//denote as C[N]
//C indicate the order of tasks are scheduled.
//ex. task C[i] is the i_th one be scheduled.
//-----------------------------------------------------------------------------
std::unique_ptr<std::vector<int16_t>> GeneticSchedulerCore::generate_new_chromosome()
{
    int16_t node_num = nodes_.size() - 1;
    std::vector<int16_t> node_order(node_num);

    node_order[0] = 0;
    for (int16_t i = 1; i < node_num; i++)
    {
        int16_t min_index = 0;
        for (int16_t j : nodes_[i].sub_nodes)
        {
            min_index = std::max(node_order[j], min_index);
        }
        min_index++;
        int16_t index = Scheduler::common::randi<int16_t>(min_index, i, seeds_[omp_get_thread_num()]);
        node_order[i] = index;
        for (int16_t j = 1; j < i; j++)
        {
            if (node_order[j] >= index)
            {
                node_order[j]++;
            }
        }
    }
    auto chromosome_ptr = common::make_unique<std::vector<int16_t>>(Scheduler::common::sort_indexes<int16_t>(node_order));
    return chromosome_ptr;
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

int16_t GeneticSchedulerCore::evaluate(std::unique_ptr<std::vector<int16_t>> &chromosome_ptr)
{
    std::vector<int16_t> nodes_finish_time(nodes_.size());
    std::vector<int16_t> cores_ocuppied_time(config_ptr_->all_core_num);

    std::fill(nodes_finish_time.begin(), nodes_finish_time.end(), inf);
    std::fill(cores_ocuppied_time.begin(), cores_ocuppied_time.end(), 0);

    for (int16_t i : *chromosome_ptr)
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

//----------------------------------------------------------------------------
//New generation created by crossover.
//The new chromosome will mutated by a certain probability.
//----------------------------------------------------------------------------
void GeneticSchedulerCore::mutation(std::unique_ptr<std::vector<int16_t>>& chromosome_ptr)
{
    std::vector<int16_t> node_order = Scheduler::common::sort_indexes<int16_t>(*chromosome_ptr);
    int16_t node_num = nodes_.size() - 1;
    for (int16_t i = 1; i < node_num; i++)
    {
        float r = static_cast <float> (rand_r(&seeds_[omp_get_thread_num()])) / static_cast <float> (RAND_MAX);
        
        if (r < config_ptr_->mutation_rate )
        {
            int16_t min_order = 0;
            for (int16_t j : nodes_[i].sub_nodes)
            {
                min_order = std::max(node_order[j], min_order);
            }
            min_order++;
            int16_t current_order = node_order[i];
            assert(min_order <= current_order);
            int16_t new_order = Scheduler::common::randi<int16_t>(min_order, current_order, seeds_[omp_get_thread_num()]);
            if(new_order == current_order){
                continue;
            }
            for (int16_t j = current_order - 1; j >= new_order; j--)
            {
                int16_t node_to_be_moved = (*chromosome_ptr)[j];
                (*chromosome_ptr)[j+1] = node_to_be_moved;
                node_order[node_to_be_moved]++;
            }
            (*chromosome_ptr)[new_order] = i;
            node_order[i] = new_order;
        }
    }
}
//----------------------------------------------------------------------------
//Randomly select two chromosome as parent.
//Exchange part of the gene of the selected chromosome, and generate a new one.
//Due to the special definition of our chromosome, we designed a unique crossover mechanism...
//----------------------------------------------------------------------------
std::unique_ptr<std::vector<int16_t>> GeneticSchedulerCore::crossover(std::unique_ptr<std::vector<int16_t>> &chromosome_a_ptr,
                                                                      std::unique_ptr<std::vector<int16_t>> &chromosome_b_ptr,
                                                                      uint32_t &seed)
{
    auto chromosome_ptr = common::make_unique<std::vector<int16_t>>(*chromosome_a_ptr); 
    int16_t node_num = (*chromosome_a_ptr).size();
    std::vector<bool> mark(node_num);
    int16_t cutting_point = Scheduler::common::randi<int16_t>(1, node_num, seed);

    for(int16_t i = 0; i<cutting_point; i++)
    {
        mark[(*chromosome_a_ptr)[i]] = true;
    }
    for(int16_t i = 0; i<node_num; i++)
    {
        if(mark[(*chromosome_b_ptr)[i]]){
            continue;
        }
        (*chromosome_ptr)[cutting_point++] = (*chromosome_b_ptr)[i];
    }

    return chromosome_ptr;
}
}
