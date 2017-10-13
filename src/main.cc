#include "FileReader.h"
#include "GeneticSchedulerCore.h"

int main(int argc, char *argv[])
{
    if (argc == 1)
    {
        auto reader = Scheduler::FileReader(std::string("/home/liu/workspace/Genetic_Scheduler/STG/50/rand0000.stg"));
        auto p0 = std::chrono::system_clock::now();
        Scheduler::GeneticSchedulerCore scheduler(reader.getNodes(), 4, pow(2, 14), 8, 0.05, 50);
        auto p1 = std::chrono::system_clock::now();
        auto diff1 = p1 - p0;
        std::cout << "time:" << std::chrono::duration_cast<std::chrono::milliseconds>(diff1).count() << " millisec.  "
                  << "best: " << scheduler.getBest() << std::endl;
    }
    else if(argc == 2)
    {
        auto reader = Scheduler::FileReader(argv[1]);        
        auto p0 = std::chrono::system_clock::now();
        Scheduler::GeneticSchedulerCore scheduler(reader.getNodes(), 4, pow(2, 14), 8, 0.05, 50);
        auto p1 = std::chrono::system_clock::now();
        auto diff1 = p1 - p0;
        std::cout << "time:" << std::chrono::duration_cast<std::chrono::milliseconds>(diff1).count() << " millisec.  "
                  << "best: " << scheduler.getBest() << std::endl;
    }
}
