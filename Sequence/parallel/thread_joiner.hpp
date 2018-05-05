#ifndef SEQUENCE_PARALLEL_THREAD_JOINER_HPP__
#define SEQUENCE_PARALLEL_THREAD_JOINER_HPP__

#include <vector>
#include <thread>

namespace Sequence
{
    namespace parallel
    {
        class thread_joiner
        {
          private:
            std::vector<std::thread>& threads;

          public:
            thread_joiner(std::vector<std::thread>& threads_)
                : threads(threads_)
            {
            }
            ~thread_joiner()
            {
                for (auto& t : threads)
                    {
                        t.join();
                    }
            }
        };
    }
}

#endif
