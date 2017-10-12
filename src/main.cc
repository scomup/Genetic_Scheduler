#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include "FileReader.h"
#include <assert.h>
#include <random>

int16_t inf = -1;

std::vector<int16_t> sort_indexes(const std::vector<int16_t> &v)
{
    std::vector<int16_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&v](int16_t i1, int16_t i2) { return v[i1] < v[i2]; });
    return idx;
}

int randi(int lo, int hi)
{
    std::random_device rnd;                          // 非決定的な乱数生成器を生成
    std::mt19937 mt(rnd());                          //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    std::uniform_int_distribution<> rand100(lo, hi); // [0, 99] 範囲の一様乱数
    return rand100(mt);
    /*
    int n = hi - lo + 1;
    int i = rand() % n;
    if (i < 0) i = -i;
    return lo + i;
    */
}

void generate_DNA(std::vector<Scheduler::Node> &nodes, std::vector<int16_t> &dna)
{
    int16_t node_num = nodes.size() - 1;
    dna.resize(node_num);
    for (int16_t i = 0; i < node_num; i++)
    {
        if (i == 0)
        {
            dna[i] = 0;
            continue;
        }
        int16_t min_index = 0;
        for (int16_t j : nodes[i].sub_nodes)
        {
            min_index = std::max(dna[j], min_index);
        }
        min_index++;
        int16_t index = randi(min_index, i);
        dna[i] = index;
        for (int16_t j = 1; j < i; j++)
        {
            if (dna[j] >= index)
            {
                dna[j]++;
            }
        }
    }
}

int16_t evaluate(std::vector<Scheduler::Node> &nodes, std::vector<int16_t> &dna, int16_t all_core_num)
{
    std::vector<int16_t> idx = sort_indexes(dna);
    std::vector<int16_t> nodes_finish_time(nodes.size());
    std::vector<int16_t> cores_ocuppied_time(all_core_num);

    std::fill(nodes_finish_time.begin(), nodes_finish_time.end(), inf);
    std::fill(cores_ocuppied_time.begin(), cores_ocuppied_time.end(), 0);

    for (int16_t i : idx)
    {
        int16_t schedulabe_time = 0;
        for (size_t j = 0; j < nodes[i].sub_nodes.size(); j++)
        {
            int16_t sub_node_finished_time = nodes_finish_time[nodes[i].sub_nodes[j]];
            assert(sub_node_finished_time != inf);
            schedulabe_time = std::max(sub_node_finished_time, schedulabe_time);
        }

        int16_t core_num = nodes[i].core_num;
        std::sort(cores_ocuppied_time.begin(), cores_ocuppied_time.end());
        schedulabe_time = std::max(schedulabe_time, cores_ocuppied_time[core_num - 1]);
        nodes_finish_time[i] = schedulabe_time + nodes[i].time;
        for (int16_t j = 0; j < nodes[i].core_num; j++)
        {
            cores_ocuppied_time[j] = nodes_finish_time[i];
        }
    }
    std::sort(cores_ocuppied_time.begin(), cores_ocuppied_time.end());
    return cores_ocuppied_time.back();
}

int main(int argc, char *argv[])
{
    auto reader = Scheduler::FileReader(std::string("/home/liu/workspace/Scheduler_GA/STG/50/rand0001.stg"));

    std::vector<Scheduler::Node> &nodes = reader.getNodes();

    int16_t best_result = 10000;
    uint64_t c = 0;
    while (true)
    {
        std::vector<int16_t> dna;
        generate_DNA(nodes, dna);

        /*
        for (int16_t m : dna)
            std::cout << m << " ";
        std::cout << std::endl;
        std::vector<int16_t> idx = sort_indexes(dna);
        for (int16_t m : idx)
            std::cout << m << " ";
        std::cout << std::endl;
        */

        int16_t evaluate_result = evaluate(nodes, dna, 4);
        if (evaluate_result < best_result)
        {
            best_result = evaluate_result;
            std::cout << "new:"<< best_result << std::endl;
        }
        if (c % 100000 == 0)
        {
            std::cout << c <<":"<< best_result<< std::endl;
        }
        c++;
    }

    //auto schedule_ptr = new Scheduler::Schedule(reader.getNodes(), 4);
    //schedule_ptr->solve();
    //std::cout<<"OK!"<<std::endl;
    //delete schedule_ptr;
    //auto schedule_ptr2 = new Scheduler::Schedule(reader.getNodes(), 4);
    //schedule_ptr2->solve();
    std::cout << "OK!" << std::endl;
    //std::this_thread::sleep_for(std::chrono::minutes(100));

    //auto all_unscheduled = state.get_unscheduled();
}
