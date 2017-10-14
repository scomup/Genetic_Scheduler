#include "FileReader.h"

namespace Scheduler
{

FileReader::FileReader(std::string filename, int16_t all_core_num)
    : filename_(filename), all_core_num_(all_core_num)
{
    read_node_file();
}

void FileReader::read_node_file()
{
    std::ifstream fin;
    int16_t node_num;
    fin.open(filename_);
    if (fin.fail())
    {
        std::cout << "Open file failed!\n";
        exit(1);
    }
    fin >> node_num;

    for (int16_t k = 0; k < node_num + 2; k++)
    {
        //Node node;
        int16_t id;
        int16_t runing_time;
        int16_t sub_num;
        int16_t core;
        std::vector<int16_t> sub_nodes;
        fin >> id >> runing_time >> sub_num;
        core = id % all_core_num_ + 1;
        sub_nodes.resize(sub_num);
        for (int16_t i = 0; i < sub_num; i++)
        {
            fin >> sub_nodes[i];
        }
        Node node{id, core, runing_time, sub_nodes};
        Nodes_.push_back(node);
    }
}
}
