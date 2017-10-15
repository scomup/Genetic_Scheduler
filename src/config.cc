#include "config.h"


namespace Scheduler
{
    Config::Config(std::string file_name)
    {
        std::ifstream config_file;
        config_file.open(file_name);
        if (!config_file.is_open())
        {
            std::cout << "Cannot open config file!" << std::endl;
            exit(1);
        }
        parse(config_file);
    }

    void Config::parse(std::istream &file)
    {
        {
            for (std::string line; std::getline(file, line);)
            {
                std::istringstream iss(line);
                std::string id, eq, val;

                bool error = false;

                if (!(iss >> id))
                {
                    error = true;
                }
                else if (id[0] == '#')
                {
                    continue;
                }
                else if (!(iss >> eq >> val >> std::ws) || eq != "=" || iss.get() != EOF)
                {
                    error = true;
                }

                if (error)
                {
                    std::cout << "parsing error!" << std::endl;
                    exit(1);
                }
                else
                {
                    //options_[id] = val;
                    if (id == "all_core_num")
                    {
                        all_core_num = atoi(val.c_str());
                    }
                    else if (id == "max_chromosome_num")
                    {
                        max_chromosome_num = atoi(val.c_str()); 
                    }
                    else if (id == "aphla")
                    {
                        aphla = atof(val.c_str()); 
                    }
                    else if (id == "mutation_rate")
                    {
                        mutation_rate = atof(val.c_str());  
                    }
                    else if (id == "max_loop")
                    {
                        max_loop = atoi(val.c_str()); 
                    }
                    else if (id == "use_crossover")
                    {
                        use_crossover = atoi(val.c_str()); 
                    }
                    else if (id == "show_step_info")
                    {
                        show_step_info = atoi(val.c_str()); 
                    }
                }
            }
        };
        //std::map<std::string, std::string> options_;
    }



}