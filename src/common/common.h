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
  T randi(T lo, T hi, uint32_t seed = 0)
  {
    #if 0
      std::random_device rnd;  
      std::mt19937 mt(rnd()); 
      std::uniform_int_distribution<> rand(lo, hi);
      return rand(mt);
      #else
      //srand (time(NULL));
      int n = hi - lo + 1;
      int i = rand_r(&seed) % n;
      if (i < 0) i = -i;
      return lo + i;
      #endif
      
  }
}

/*
double rand_gaussian(double sigma, uint32_t seed = 0)
{
  double x1, x2, w, r;

  do
  {
    do { r = drand48_r(&seed); } while (r==0.0);
    x1 = 2.0 * r - 1.0;
    do { r = drand48_r(&seed); } while (r==0.0);
    x2 = 2.0 * r - 1.0;
    w = x1*x1 + x2*x2;
  } while(w > 1.0 || w==0.0);

  return(sigma * x2 * sqrt(-2.0*log(w)/w));
}
*/



}

#endif
