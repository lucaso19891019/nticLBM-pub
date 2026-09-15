#pragma once

#include <algorithm>
#include <vector>


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

}