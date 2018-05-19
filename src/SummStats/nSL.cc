#include <Sequence/SimData.hpp>
#include <Sequence/SummStats/nSL.hpp>
#include <algorithm>
#include <numeric>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>
// For parallelizing nSL:
#include <functional>
#include <Sequence/parallel/transform.hpp>

using namespace std;

namespace
{
    inline double
    maxabs(double score, double mean, double sd, double rv)
    {
        if (isfinite(score))
            {
                double zscore = (score - mean) / sd;
                if (!isfinite(rv) || fabs(zscore) > fabs(rv))
                    return zscore;
            }
        return rv;
    }

    double
    update_s2(std::string::const_iterator start,
              std::string::const_iterator left,
              std::string::const_iterator right, const Sequence::SimData &d,
              const std::unordered_map<double, double> &gmap)
    {
        auto p1 = d.position(
            std::vector<double>::size_type(distance(start, right)) - 1);
        auto p2 = d.position(
            std::vector<double>::size_type(distance(start, left)));
        if (gmap.empty())
            {
                // return phyisical distance
                return fabs(p1 - p2);
            }
        // return distance along genetic map,
        // in whatever units those are.
        auto fp1 = gmap.find(p1);
        auto fp2 = gmap.find(p2);
        if (fp1 == gmap.end() || fp2 == gmap.end())
            {
                throw std::runtime_error(
                    "position could not be found in genetic map, "
                    + std::string(__FILE__) + " line "
                    + std::to_string(__LINE__));
            }
        return fabs(fp1->second - fp2->second);
    }
    /*
      Mechanics of the nSL statistic

      RV = nSL,iHS, as defined in doi:10.1093/molbev/msu077
    */
    std::array<double, 4>
    __nSLdetails(const std::size_t &core, const Sequence::SimData &d,
                 // const vector<size_t> &coretype,
                 const std::unordered_map<double, double> &gmap)
    {
        auto csize = d.size();
        // This tracks s,s2 for ancestral and derived
        // mutation, resp:
        std::array<double, 4> rv = { 0., 0., 0., 0. };
        // number of comparisons for ancestral and
        // derived, resp:
        std::array<unsigned, 2> nc = { 0u, 0u };
        for (size_t i = 0; i < csize; ++i)
            {
                auto start = d[i].cbegin();
                auto bi
                    = start
                      + static_cast<decltype(start)::difference_type>(core);
                size_t iIsDer = (*bi == '1');
                for (size_t j = i + 1; j < csize; ++j)
                    {
                        auto bj
                            = d[j].cbegin()
                              + static_cast<decltype(start)::difference_type>(
                                    core);
                        size_t jIsDer = (*bj == '1');
                        if (iIsDer == jIsDer)
                            {
                                auto eri = d[i].crend();
                                auto ei = d[i].cend();
                                auto right = mismatch(bi, ei, bj);
                                string::const_reverse_iterator ri1(bi),
                                    ri2(bj);
                                auto left = mismatch(ri1, eri, ri2);
                                if (left.first != eri && right.first != ei)
                                    {
                                        rv[2 * iIsDer] += static_cast<double>(
                                            distance(left.first.base(),
                                                     right.first)
                                            + 1);
                                        rv[2 * iIsDer + 1] += update_s2(
                                            start, left.first.base(),
                                            right.first, d, gmap);
                                        nc[iIsDer]++;
                                    }
                            }
                    }
            }
        rv[0] /= static_cast<double>(nc[0]);
        rv[1] /= static_cast<double>(nc[0]);
        rv[2] /= static_cast<double>(nc[1]);
        rv[3] /= static_cast<double>(nc[1]);
        return rv;
    }
} // namespace

namespace Sequence
{
    /*
      The nSL statistic of doi:10.1093/molbev/msu077
    */
    pair<double, double>
    nSL(const std::size_t &core, const SimData &d,
        const std::unordered_map<double, double> &gmap)
    {
        auto nsl = __nSLdetails(core, d, gmap);
        return make_pair(log(nsl[0]) - log(nsl[2]), log(nsl[1]) - log(nsl[3]));
    }

    vector<tuple<double, double, uint32_t>>
    nSL_t(const SimData &d, const std::unordered_map<double, double> &gmap)
    {
        vector<size_t> core_snps(d.numsites());
        for (size_t core = 0; core < core_snps.size(); ++core)
            core_snps[core] = core;
        return nSL_t(d, core_snps, gmap);
    }

    vector<tuple<double, double, uint32_t>>
    nSL_t(const SimData &d, const std::vector<size_t> &core_snps,
          const std::unordered_map<double, double> &gmap,
          const int nthreads = -1)
    {
        using offset_type = SimData::const_site_iterator::difference_type;
        std::vector<std::tuple<double, double, uint32_t>> rv(core_snps.size());
        parallel::transform(
            core_snps.begin(), core_snps.end(), rv.begin(),
            [&d, &gmap](const std::size_t i) {
                auto temp = __nSLdetails(i, d, gmap);
                std::uint32_t dcount = static_cast<std::uint32_t>(
                    std::count((d.sbegin() + i)->second.begin(),
                               (d.sbegin() + i)->second.end(), '1'));
                return std::make_tuple(std::log(temp[0]) - std::log(temp[2]),
                                       std::log(temp[1]) - std::log(temp[3]),
                                       dcount);
            },
            nthreads);
        return rv;
    }
} // namespace Sequence
