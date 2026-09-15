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
// The sorting network operates on the next power-of-two number of positions.
// Positions beyond the original vector size are treated as virtual padding
// elements that compare greater than every real element and therefore migrate
// to the end of the network.
//
// Compare-exchange operations within each network stage are independent and
// are executed in parallel with OpenMP.
//=============================================================================

template<typename T, typename Compare>
void parallelSort(
    std::vector<T>& data,
    Compare compare)
{
    const std::size_t originalSize =
        data.size();


    if (originalSize <= 1)
    {
        return;
    }


    //-------------------------------------------------------------------------
    // Determine the next power-of-two network size.
    //-------------------------------------------------------------------------

    std::size_t networkSize = 1;


    while (networkSize < originalSize)
    {
        networkSize <<= 1;
    }


    //-------------------------------------------------------------------------
    // Allocate the sorting network.
    //
    // A separate validity array marks real and padding positions. This avoids
    // requiring a sentinel value or any additional property of T.
    //-------------------------------------------------------------------------

    std::vector<T> network(
        networkSize);

    std::vector<unsigned char> valid(
        networkSize,
        0);


#pragma omp parallel for
    for (std::size_t i = 0;
         i < originalSize;
         ++i)
    {
        network[i] =
            data[i];

        valid[i] =
            1;
    }


    //-------------------------------------------------------------------------
    // Execute the bitonic sorting network.
    //-------------------------------------------------------------------------

    for (std::size_t sequenceSize = 2;
         sequenceSize <= networkSize;)
    {
        for (std::size_t stride = sequenceSize >> 1;
             stride > 0;
             stride >>= 1)
        {
#pragma omp parallel for
            for (std::size_t i = 0;
                 i < networkSize;
                 ++i)
            {
                const std::size_t partner =
                    i ^ stride;


                if (partner <= i)
                {
                    continue;
                }


                const bool ascending =
                    ((i & sequenceSize) == 0);


                bool swapRequired = false;


                if (valid[i] != valid[partner])
                {
                    if (ascending)
                    {
                        swapRequired =
                            (valid[i] == 0);
                    }
                    else
                    {
                        swapRequired =
                            (valid[partner] == 0);
                    }
                }
                else if (valid[i] != 0)
                {
                    if (ascending)
                    {
                        swapRequired =
                            compare(
                                network[partner],
                                network[i]);
                    }
                    else
                    {
                        swapRequired =
                            compare(
                                network[i],
                                network[partner]);
                    }
                }


                if (swapRequired)
                {
                    std::swap(
                        network[i],
                        network[partner]);

                    std::swap(
                        valid[i],
                        valid[partner]);
                }
            }
        }


        if (sequenceSize == networkSize)
        {
            break;
        }


        sequenceSize <<= 1;
    }


    //-------------------------------------------------------------------------
    // Copy the sorted real elements back to the original vector.
    //-------------------------------------------------------------------------

#pragma omp parallel for
    for (std::size_t i = 0;
         i < originalSize;
         ++i)
    {
        data[i] =
            std::move(
                network[i]);
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