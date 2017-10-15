#ifndef FILEREADER_H
#define FILEREADER_H

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

#include "node.h"

namespace Scheduler
{

class FileReader
{
  public:
    FileReader(std::string, int16_t all_core_num_);
    inline std::vector<Node>& getNodes(){return Nodes_;};

  private:
    void read_node_file();

    std::vector<Node> Nodes_;
    std::string filename_;
    int16_t all_core_num_;
};

}

#endif