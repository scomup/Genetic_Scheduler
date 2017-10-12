#include <assert.h>

#include "common/common.h"
#include "FileReader.h"
#include "GeneticSchedulerCore.h"

int16_t inf = -1;


int main(int argc, char *argv[])
{
    srand(time(NULL));
    auto reader = Scheduler::FileReader(std::string("/home/liu/workspace/Genetic_Scheduler/STG/50/rand0007.stg"));

    Scheduler::GeneticSchedulerCore scheduler(reader.getNodes(), 4,  100000,  100,  0.2);
    
    std::cout << "OK!" << std::endl;

}
