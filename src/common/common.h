#ifndef COMMON_COMMON_H_
#define COMMON_COMMON_H_

#include <algorithm>
#include <functional>
#include <random>

namespace Scheduler
{
namespace common
{


template <typename T>
std::vector<T> sort_indexes(const std::vector<T> &v)
{

  // initialize original index locations
  std::vector<T> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](T i1, T i2) { return v[i1] < v[i2]; });

  return idx;
}

template <typename T>
T randi(T lo, T hi, uint32_t& seed)
{
  int n = hi - lo + 1;
  int i = rand_r(&seed) % n;
  if (i < 0)
    i = -i;
  return lo + i;
}


template <typename T>
T rand_normal(T mean, T stddev, uint32_t& seed)
{ //Box muller method
  static T n2 = 0.0;
  static int n2_cached = 0;
  if (!n2_cached)
  {
    T x, y, r;
    do
    {
      x = 2.0 * static_cast <float> (rand_r(&seed)) / static_cast <float> (RAND_MAX) - 1;
      y = 2.0 * static_cast <float> (rand_r(&seed)) / static_cast <float> (RAND_MAX) - 1;
      r = x * x + y * y;
    } while (r == 0.0 || r > 1.0);
    {
      T d = sqrt(-2.0 * log(r) / r);
      T n1 = x * d;
      n2 = y * d;
      T result = n1 * stddev + mean;
      n2_cached = 1;
      return result;
    }
  }
  else
  {
    n2_cached = 0;
    return n2 * stddev + mean;
  }
}


}
}

#endif
