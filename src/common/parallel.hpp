#pragma once

#include <algorithm>
#include <vector>
#include <cstddef>

namespace ntic::lbm::common
{

template<typename T>
void parallelSort(
    std::vector<T>& data)
{
    std::sort(
        data.begin(),
        data.end());
}


template<typename T, typename Compare>
void parallelSort(
    std::vector<T>& data,
    Compare compare)
{
    std::sort(
        data.begin(),
        data.end(),
        compare);
}

template<typename InputIt,
         typename OutputIt,
         typename T>
void parallelScan(
    InputIt first,
    InputIt last,
    OutputIt result,
    T init)
{
    T sum = init;

    for (; first != last; ++first, ++result)
    {
        *result = sum;
        sum += *first;
    }
}

template<typename T>
std::size_t parallelScan(
    const std::vector<T>& input,
    std::vector<std::size_t>& offsets)
{
    offsets.resize(input.size());

    std::size_t sum = 0;

    for (std::size_t i = 0;
         i < input.size();
         ++i)
    {
        offsets[i] = sum;
        sum += static_cast<std::size_t>(input[i]);
    }

    return sum;
}

}