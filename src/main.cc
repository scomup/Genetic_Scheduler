#include "common/common.h"
#include "FileReader.h"
#include "GeneticSchedulerCore.h"

#include <time.h>
#include <random>
#include<stdio.h>
#include<omp.h>


int main(int argc, char *argv[])
{
    auto reader = Scheduler::FileReader(std::string("/home/liu/workspace/Genetic_Scheduler/STG/50/rand0007.stg"));
    auto p0 = std::chrono::system_clock::now();
    Scheduler::GeneticSchedulerCore scheduler(reader.getNodes(), 4,  pow(2,14),  8,  0.05, 50);
    auto p1 = std::chrono::system_clock::now();
    auto diff1 = p1 - p0;
    std::cout << "time:" << std::chrono::duration_cast<std::chrono::milliseconds>(diff1).count() << " millisec.  " << "best: "<<scheduler.getBest()<< std::endl;
}
