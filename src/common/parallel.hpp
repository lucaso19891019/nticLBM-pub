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

}