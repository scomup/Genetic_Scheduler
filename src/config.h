#ifndef CONFIG_H
#define CONFIG_H


#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

namespace Scheduler
{

class Config
{
  public:
    Config(std::string file_name);

    uint32_t all_core_num;
    uint32_t max_chromosome_num;
    float aphla;
    float mutation_rate;
    uint32_t max_loop;
    bool use_crossover;
    bool show_step_info;  

  private:
    void parse(std::istream &file);
};
}

#endif
