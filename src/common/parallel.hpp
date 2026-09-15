#pragma once

#include <algorithm>
#include <vector>
#include <cstddef>
#include <omp.h>


namespace ntic::lbm::common
{

//=============================================================================
// Sort a vector in parallel using a bitonic sorting network.
//
// The input is temporarily padded to the next power-of-two size required by
// the bitonic sorting network. Compare-exchange operations at each network
// stage are independent and are executed in parallel with OpenMP.
//
// The original vector size is restored after sorting.
//=============================================================================

template<typename T, typename Compare>
void parallelSort(
    std::vector<T>& data,
    Compare compare)
{
    const std::size_t originalSize =
        data.size();


    if(originalSize <= 1)
    {
        return;
    }


    //-------------------------------------------------------------------------
    // Find next power of two.
    //-------------------------------------------------------------------------

    std::size_t size = 1;

    while(size < originalSize)
    {
        size *= 2;
    }


    T maxValue =
        data[0];

    for(const auto& value : data)
    {
        if(compare(maxValue,value))
        {
            maxValue = value;
        }
    }


    data.resize(
        size,
        maxValue);



    //-------------------------------------------------------------------------
    // Bitonic sorting network.
    //-------------------------------------------------------------------------

    for(std::size_t k = 2;
        k <= size;
        k *= 2)
    {
        for(std::size_t j = k / 2;
            j > 0;
            j /= 2)
        {

#pragma omp parallel for
            for(std::size_t i = 0;
                i < size;
                ++i)
            {
                const std::size_t ixj =
                    i ^ j;


                if(ixj > i)
                {
                    const bool ascending =
                        ((i & k) == 0);


                    if(ascending)
                    {
                        if(compare(
                               data[ixj],
                               data[i]))
                        {
                            std::swap(
                                data[i],
                                data[ixj]);
                        }
                    }
                    else
                    {
                        if(compare(
                               data[i],
                               data[ixj]))
                        {
                            std::swap(
                                data[i],
                                data[ixj]);
                        }
                    }
                }
            }
        }
    }


    data.resize(
        originalSize);
}



//=============================================================================
// Parallel exclusive scan
//
// Example:
// input:
//     [3,5,2,7]
//
// output:
//     [0,3,8,10]
//
// return:
//     17
//
// Used by STL geometry construction:
// - thread local vector merging
// - representative ID generation
//=============================================================================

template<typename T>
std::size_t parallelScan(
    const std::vector<T>& input,
    std::vector<std::size_t>& offsets)
{
    const std::size_t size =
        input.size();


    offsets.resize(
        size);


    if(size == 0)
    {
        return 0;
    }


    const int threadCount =
        omp_get_max_threads();


    std::vector<std::size_t> partial(
        static_cast<std::size_t>(threadCount),
        0);



    //-------------------------------------------------------------------------
    // Phase 1:
    // Each thread computes local prefix sum.
    //-------------------------------------------------------------------------

    #pragma omp parallel
    {
        const int tid =
            omp_get_thread_num();


        const std::size_t begin =
            static_cast<std::size_t>(tid) *
            size /
            static_cast<std::size_t>(threadCount);


        const std::size_t end =
            static_cast<std::size_t>(tid + 1) *
            size /
            static_cast<std::size_t>(threadCount);


        std::size_t localSum = 0;


        for(std::size_t i = begin;
            i < end;
            ++i)
        {
            offsets[i] =
                localSum;

            localSum +=
                static_cast<std::size_t>(
                    input[i]);
        }


        partial[
            static_cast<std::size_t>(tid)]
            =
            localSum;
    }



    //-------------------------------------------------------------------------
    // Phase 2:
    // Compute offset of every thread block.
    //-------------------------------------------------------------------------

    std::size_t total = 0;


    for(int t = 0;
        t < threadCount;
        ++t)
    {
        const std::size_t tmp =
            partial[
                static_cast<std::size_t>(t)];

        partial[
            static_cast<std::size_t>(t)]
            =
            total;

        total += tmp;
    }



    //-------------------------------------------------------------------------
    // Phase 3:
    // Add block offsets.
    //-------------------------------------------------------------------------

    #pragma omp parallel
    {
        const int tid =
            omp_get_thread_num();


        const std::size_t add =
            partial[
                static_cast<std::size_t>(tid)];


        const std::size_t begin =
            static_cast<std::size_t>(tid) *
            size /
            static_cast<std::size_t>(threadCount);


        const std::size_t end =
            static_cast<std::size_t>(tid + 1) *
            size /
            static_cast<std::size_t>(threadCount);



        for(std::size_t i = begin;
            i < end;
            ++i)
        {
            offsets[i] += add;
        }
    }


    return total;
}


}