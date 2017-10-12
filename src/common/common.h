#ifndef COMMON_COMMON_H_
#define COMMON_COMMON_H_

#include <algorithm>
#include <functional>
#include <random>

namespace Scheduler {
namespace common {

  template <typename T>
  std::vector<T> sort_indexes(const std::vector<T> &v) {
  
    // initialize original index locations
    std::vector<T> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
  
    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
         [&v](T i1, T i2) {return v[i1] < v[i2];});
  
    return idx;
  }

  template <typename T>
  T randi(T lo, T hi)
  {
    #if 0
      std::random_device rnd;  
      std::mt19937 mt(rnd()); 
      std::uniform_int_distribution<> rand(lo, hi);
      return rand(mt);
      #else
      //srand (time(NULL));
      int n = hi - lo + 1;
      int i = rand() % n;
      if (i < 0) i = -i;
      return lo + i;
      #endif
      
  }
}


}

#endif
