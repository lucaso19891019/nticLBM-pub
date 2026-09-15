#pragma once

#include <algorithm>
#include <vector>
#include <cstddef>
#include <omp.h>


namespace ntic::lbm::common
{

//=============================================================================
// Sort a vector using an OpenMP task-based parallel merge sort.
//
// The algorithm recursively divides the input range into smaller sorted
// segments. Independent segments are processed concurrently using OpenMP
// tasks. The sorted segments are combined using a parallel merge-path based
// merge operation.
//
// Unlike a traditional merge sort, the merge stage is also parallelized.
// Each worker is assigned an independent output range, eliminating the serial
// merge bottleneck.
//
// The implementation supports arbitrary data types and user-defined comparison
// operators.
//=============================================================================

template<typename T, typename Compare>
void parallelSort(
    std::vector<T>& data,
    Compare compare)
{
    const std::size_t size =
        data.size();


    if (size <= 1)
    {
        return;
    }


    std::vector<T> buffer(
        size);



    //-------------------------------------------------------------------------
    // Merge-path partition.
    //
    // Given two sorted ranges:
    //
    //     A = [aBegin, aEnd)
    //     B = [bBegin, bEnd)
    //
    // Find the partition point on diagonal k:
    //
    //     i + j = k
    //
    // such that:
    //
    //     A[i-1] <= B[j]
    //     B[j-1] <  A[i]
    //
    //-------------------------------------------------------------------------
    
    auto mergePath =
        [&](std::size_t k,
            std::size_t aBegin,
            std::size_t aEnd,
            std::size_t bBegin,
            std::size_t bEnd)
        {
            std::size_t low =
                (k > (bEnd - bBegin))
                ?
                k - (bEnd - bBegin)
                :
                0;


            std::size_t high =
                std::min(
                    k,
                    aEnd - aBegin);



            while (low < high)
            {
                const std::size_t i =
                    (low + high) / 2;


                const std::size_t j =
                    k - i;


                if (i < aEnd - aBegin &&
                    j > 0 &&
                    compare(
                        data[aBegin + i],
                        data[bBegin + j - 1]))
                {
                    low =
                        i + 1;
                }
                else
                {
                    high =
                        i;
                }
            }


            return low;
        };



    //-------------------------------------------------------------------------
    // Parallel merge.
    //
    // The output range is divided into independent diagonal segments.
    // Each thread performs a local sequential merge inside its segment.
    //-------------------------------------------------------------------------

    auto parallelMerge =
        [&](std::size_t begin,
            std::size_t middle,
            std::size_t end)
        {
            const std::size_t leftSize =
                middle - begin;


            const std::size_t rightSize =
                end - middle;


            const std::size_t total =
                leftSize + rightSize;


            const int threads =
                omp_get_max_threads();


            const std::size_t chunk =
                (total + threads - 1)
                /
                threads;



#pragma omp parallel
            {
                const int tid =
                    omp_get_thread_num();


                const std::size_t outBegin =
                    std::min(
                        static_cast<std::size_t>(tid) * chunk,
                        total);


                const std::size_t outEnd =
                    std::min(
                        outBegin + chunk,
                        total);


                if (outBegin < outEnd)
                {
                    const std::size_t a0 =
                        mergePath(
                            outBegin,
                            begin,
                            middle,
                            middle,
                            end);


                    const std::size_t b0 =
                        outBegin - a0;

                    std::size_t i =
                        begin + a0;


                    std::size_t j =
                        middle + b0;


                    std::size_t k =
                        begin + outBegin;


                    while (k < begin + outEnd)
                    {
                        if (i < middle &&
                            (j >= end ||
                             compare(
                                 data[i],
                                 data[j])))
                        {
                            buffer[k++] =
                                std::move(
                                    data[i++]);
                        }
                        else
                        {
                            buffer[k++] =
                                std::move(
                                    data[j++]);
                        }
                    }
                }
            }


#pragma omp parallel for
            for (std::size_t i = begin;
                 i < end;
                 ++i)
            {
                data[i] =
                    std::move(
                        buffer[i]);
            }
        };



    //-------------------------------------------------------------------------
    // Recursive parallel merge sort.
    //-------------------------------------------------------------------------

    auto mergeSort =
        [&](auto&& self,
            std::size_t begin,
            std::size_t end,
            int depth) -> void
        {
            const std::size_t length =
                end - begin;


            constexpr std::size_t threshold =
                2048;


            if (length <= threshold)
            {
                std::sort(
                    data.begin() + begin,
                    data.begin() + end,
                    compare);

                return;
            }


            const std::size_t middle =
                begin + length / 2;


#pragma omp task if(depth < 4)
            {
                self(
                    self,
                    begin,
                    middle,
                    depth + 1);
            }


#pragma omp task if(depth < 4)
            {
                self(
                    self,
                    middle,
                    end,
                    depth + 1);
            }


#pragma omp taskwait


            parallelMerge(
                begin,
                middle,
                end);
        };



#pragma omp parallel
    {
#pragma omp single
        {
            mergeSort(
                mergeSort,
                0,
                size,
                0);
        }
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