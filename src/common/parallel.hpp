#pragma once

#include <algorithm>
#include <vector>
#include <cstddef>
#include <omp.h>


namespace ntic::lbm::common
{


//=============================================================================
// Parallel sort
//
// Strategy:
// 1. Split input vector into thread-local chunks.
// 2. Sort every chunk independently in parallel.
// 3. Merge sorted chunks serially.
//
// The interface remains identical to std::sort style.
//=============================================================================

template<typename T, typename Compare>
void parallelSort(
    std::vector<T>& data,
    Compare compare)
{
    const std::size_t size =
        data.size();


    if(size < 1024)
    {
        std::sort(
            data.begin(),
            data.end(),
            compare);

        return;
    }


    const int threadCount =
        omp_get_max_threads();


    const std::size_t chunkSize =
        (size +
         static_cast<std::size_t>(threadCount) -
         1)
        /
        static_cast<std::size_t>(threadCount);



    //-------------------------------------------------------------------------
    // Sort independent chunks.
    //-------------------------------------------------------------------------

    #pragma omp parallel for schedule(static)
    for(int t = 0;
        t < threadCount;
        ++t)
    {
        const std::size_t begin =
            static_cast<std::size_t>(t) *
            chunkSize;


        const std::size_t end =
            std::min(
                begin + chunkSize,
                size);


        if(begin < end)
        {
            std::sort(
                data.begin() + begin,
                data.begin() + end,
                compare);
        }
    }



    //-------------------------------------------------------------------------
    // Merge sorted chunks.
    //
    // This stage is currently serial.
    // It keeps the implementation simple and deterministic.
    // It can later be replaced by parallel merge tree.
    //-------------------------------------------------------------------------

    std::vector<T> buffer(
        size);


    for(std::size_t width = chunkSize;
        width < size;
        width *= 2)
    {
        for(std::size_t begin = 0;
            begin < size;
            begin += 2 * width)
        {
            const std::size_t middle =
                std::min(
                    begin + width,
                    size);


            const std::size_t end =
                std::min(
                    begin + 2 * width,
                    size);


            if(middle >= end)
            {
                continue;
            }


            std::merge(
                data.begin() + begin,
                data.begin() + middle,
                data.begin() + middle,
                data.begin() + end,
                buffer.begin() + begin,
                compare);
        }


        data.swap(buffer);
    }
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