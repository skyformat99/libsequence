#ifndef SEQUENCE_PARALLEL_TRANSFORM_HPP__
#define SEQUENCE_PARALLEL_TRANSFORM_HPP__

#include <algorithm>
#include <thread>
#include <future>
#include <type_traits>
#include "thread_joiner.hpp"

namespace Sequence
{
    namespace parallel
    {
        template <typename InputIterator, typename OutputIterator,
                  typename Function>
        OutputIterator
        transform(InputIterator first1, InputIterator last1,
                  OutputIterator first2, Function f, const int maxthreads = -1)
        {
            using value = typename InputIterator::value_type;
            using return_type = typename std::result_of<Function(value)>::type;
            if (maxthreads == 1)
                return std::transform(first1, last1, first2, f);
            auto dist = std::distance(first1, last1);
            if (!dist)
                return first2;

            unsigned nthreads = maxthreads;
            if (maxthreads < 1)
                {
                    int min_job_load = 25;
                    int max_threads_from_load
                        = (dist + min_job_load - 1) / min_job_load;
                    int available_threads
                        = std::thread::hardware_concurrency();
                    nthreads
                        = std::min(max_threads_from_load, available_threads);
                }
            auto stepsize = dist / nthreads;

            std::vector<std::thread> threads(nthreads - 1);
            thread_joiner joiner(threads);
            std::vector<std::future<void>> futures(nthreads - 1);
            auto start_in = first1;
            auto start_out = first2;
            for (int i = 0; i < nthreads - 1; ++i)
                {
                    auto in_chunk_end = start_in;
                    std::advance(in_chunk_end, stepsize);
                    std::packaged_task<void(void)> pt(
                        [start_in, in_chunk_end, start_out, f]() -> void {
                            std::transform(start_in, in_chunk_end, start_out,
                                           f);
                        });
                    futures[i] = pt.get_future();
                    threads[i] = std::thread(std::move(pt));
                    start_in = in_chunk_end;
                    std::advance(start_out, stepsize);
                }
            auto rv = std::transform(start_in, last1, start_out, f);
            for (auto& f : futures)
                {
                    f.get();
                }
            return rv;
        }
    } // namespace parallel
} // namespace Sequence

#endif
