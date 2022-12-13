/*****************************************************************************************
 *              MIT License                                                              *
 *                                                                                       *
 * Copyright (c) 2022 G. Cherchi, M. Livesu, R. Scateni, M. Attene and F. Pellacini      *
 *                                                                                       *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this  *
 * software and associated documentation files (the "Software"), to deal in the Software *
 * without restriction, including without limitation the rights to use, copy, modify,    *
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to    *
 * permit persons to whom the Software is furnished to do so, subject to the following   *
 * conditions:                                                                           *
 *                                                                                       *
 * The above copyright notice and this permission notice shall be included in all copies *
 * or substantial portions of the Software.                                              *
 *                                                                                       *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,   *
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT    *
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION     *
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE        *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                *
 *                                                                                       *
 * Authors:                                                                              *
 *      Gianmarco Cherchi (g.cherchi@unica.it)                                           *
 *      https://www.gianmarcocherchi.com                                                 *
 *                                                                                       *
 *      Marco Livesu (marco.livesu@ge.imati.cnr.it)                                      *
 *      http://pers.ge.imati.cnr.it/livesu/                                              *
 *                                                                                       *
 *      Riccardo Scateni (riccardo@unica.it)                                             *
 *      https://people.unica.it/riccardoscateni/                                         *
 *                                                                                       *
 *      Marco Attene (marco.attene@ge.imati.cnr.it)                                      *
 *      https://www.cnr.it/en/people/marco.attene/                                       *
 *                                                                                       *
 *      Fabio Pellacini (fabio.pellacini@uniroma1.it)                                    *
 *      https://pellacini.di.uniroma1.it                                                 *
 *                                                                                       *
 * ***************************************************************************************/

#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <deque>
#include <algorithm>

#include <absl/container/flat_hash_map.h>

template<typename T>
inline void remove_duplicates(std::vector<T>& values) {
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
}

template<typename T, size_t N>
inline size_t remove_duplicates(std::array<T, N>& values) {
  std::sort(values.begin(), values.end());
  return (size_t)(std::unique(values.begin(), values.end()) - values.begin());
}

template<typename T>
inline void parallel_remove_duplicates(std::vector<T>& values) {
  tbb::parallel_sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
}

template<typename T>
inline T my_unique(std::vector<T>& data) {
  if(true) {
    return std::unique(data.begin(), data.end());
  } else {
    auto to_keep = std::vector<bool>(data.size());
    for(auto idx = 0; idx < data.size(); idx++) {
      to_keep[idx] = idx == 0 || data[idx] != data[idx-1];
    }
    auto count = 0;
    for(auto idx = 0; idx < data.size(); idx++) {
      if(to_keep[idx]) data[count++] = data[idx];
    }
    return data.begin() + count;
  }
}

template<typename T>
inline void fast_remove_duplicates(std::vector<T>& values) {
  tbb::parallel_sort(values.begin(), values.end());
  values.erase(my_unique(values), values.end());
}

template<typename T>
inline bool contains(const std::vector<T>& values, T value) {
  for(auto& value_ : values) if(value == value_) return true; 
  return false;
}

namespace std {
template<typename T, size_t N>
struct hash<std::array<T, N>>
{
    std::hash<T> hasher;
    inline size_t operator()(const std::array<T, N> &p) const
    {
        size_t seed = 0;
        for(auto i = 0; i < N; i ++) 
          seed ^= hasher(p[i]) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
    }
};

template<typename T1, typename T2>
struct hash<std::pair<T1, T2>>
{
    std::hash<T1> hasher1;
    std::hash<T2> hasher;
    inline size_t operator()(const std::pair<T1, T2> &p) const
    {
        size_t seed = 0;
        seed ^= hasher(p.first) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= hasher(p.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
    }
};

template<typename T>
struct hash<std::vector<T>>
{
    std::hash<T> hasher;
    inline size_t operator()(const std::vector<T> &p) const
    {
        size_t seed = 0;
        for(auto i = 0; i < p.size(); i ++) 
        seed ^= hasher(p[i]) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
    }
};
}

#include <absl/container/inlined_vector.h>

template<typename T, size_t N>
inline void remove_duplicates(absl::InlinedVector<T, N>& values) {
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
}

template<typename T, size_t N>
inline bool contains(const absl::InlinedVector<T, N>& values, T value) {
  for(auto& value_ : values) if(value == value_) return true; 
  return false;
}

#if 1
/*
* 这是一个二维数组，行是固定大小的，大小由N来确定，填满一行就新增一行，后面的数据就添加到新行里面。
*/
template<typename T, size_t N>
struct bucket_arena {
  std::vector<std::vector<T>> buckets;

  bucket_arena() {
    buckets.reserve(16);
  }

  template<typename ... Args>
  T& emplace_back(Args&& ... args) {
    if(buckets.empty() || buckets.back().capacity() == buckets.back().size()) {
      auto& bucket = buckets.emplace_back();
      bucket.reserve(N);
      return bucket.emplace_back(std::forward<Args>(args)...);
    } else {
      auto& bucket = buckets.back();
      return bucket.emplace_back(std::forward<Args>(args)...);
    }
  }

  void pop_back() {
    buckets.back().pop_back();
    if(buckets.back().empty()) buckets.pop_back();
  }
};

/*
* 能装下各种类型的顶点的容器
* init: 存储去除了重复顶点之后的顶点数据，其他地方会索引这些顶点数据。
*/
struct point_arena {
  std::vector<explicitPoint3D> init;//保存去除了重复顶点的顶点，其他顶点指针指向这里
  bucket_arena<implicitPoint3D_LPI, 1024 * 1024> edges;//由直线与平面相交确定的隐式点
  bucket_arena<explicitPoint3D, 1024> jolly;
  bucket_arena<implicitPoint3D_TPI, 1024 * 1024> tpi;//由三个平面相交确定的隐式点
};

#else

struct point_arena {
  std::vector<explicitPoint3D> init;
  std::deque<implicitPoint3D_LPI> edges;
  std::deque<explicitPoint3D> jolly;
  std::deque<implicitPoint3D_TPI> tpi;
};

#endif

#endif
