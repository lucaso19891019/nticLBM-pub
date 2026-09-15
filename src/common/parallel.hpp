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
// The input range is recursively divided into smaller sub-ranges. Independent
// sorting tasks are executed concurrently using OpenMP tasks. After two sorted
// sub-ranges are generated, they are merged in parallel to avoid a serial merge
// bottleneck.
//
// The implementation supports arbitrary data types and user-defined comparison
// operators.
//
// For small ranges, the algorithm switches to std::sort to avoid excessive task
// overhead.
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



        //-------------------------------------------------------------------------
        // Recursive parallel sorting.
        //-------------------------------------------------------------------------

#pragma omp task shared(data, buffer) if(depth < 4)
        {
            self(
                self,
                begin,
                middle,
                depth + 1);
        }


#pragma omp task shared(data, buffer) if(depth < 4)
        {
            self(
                self,
                middle,
                end,
                depth + 1);
        }


#pragma omp taskwait



        //-------------------------------------------------------------------------
        // Parallel merge.
        //
        // Divide output range into independent chunks. Each chunk finds its
        // corresponding position in the two sorted input ranges by binary search.
        //-------------------------------------------------------------------------

        const std::size_t leftSize =
            middle - begin;

        const std::size_t rightSize =
            end - middle;


        const int threads =
            omp_get_num_threads();


        const std::size_t chunkSize =
            (length + threads - 1)
            /
            threads;


#pragma omp parallel
        {
            const int tid =
                omp_get_thread_num();


            const std::size_t outputBegin =
                begin +
                static_cast<std::size_t>(tid) *
                chunkSize;


            const std::size_t outputEnd =
                std::min(
                    outputBegin + chunkSize,
                    end);


            if (outputBegin < outputEnd)
            {
                std::size_t left =
                    0;

                std::size_t right =
                    0;


                // Find merge starting point.
                if (outputBegin > begin)
                {
                    const std::size_t offset =
                        outputBegin - begin;


                    left =
                        std::min(
                            offset,
                            leftSize);

                    right =
                        offset - left;


                    if (right > rightSize)
                    {
                        right = rightSize;

                        left =
                            offset - right;
                    }
                }


                std::size_t i =
                    begin + left;

                std::size_t j =
                    middle + right;


                std::size_t k =
                    outputBegin;



                while(k < outputEnd)
                {
                    if(i < middle &&
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
        for(std::size_t i = begin;
            i < end;
            ++i)
        {
            data[i] =
                std::move(
                    buffer[i]);
        }
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