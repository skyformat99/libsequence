#include <cmath>
#include <limits>
#include <Sequence/summstats/algorithm.hpp>
#include <Sequence/StateCounts.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/auxillary.hpp>

namespace Sequence
{
    double
    tajd(const VariantMatrix& m)
    {
        unsigned S = 0;
        double pi = 0.0;
        const auto agg = [&S, &pi](const StateCounts& c) {
            unsigned nstates = 0;
            double homozygosity = 0.0;
            for (auto i : c.counts)
                {
                    if (i)
                        {
                            nstates++;
                            homozygosity += static_cast<double>(i)
                                            * static_cast<double>(i - 1);
                        }
                }
            if (nstates > 1)
                {
                    S += nstates - 1;
                    double d = static_cast<double>(c.n * (c.n - 1));
                    pi += (1.0 - homozygosity / d);
                }
        };
        sstats_algo::aggregate_sites(m, agg, -1);
        if (!S)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        auto a1 = summstats_aux::a_sub_n(static_cast<std::uint32_t>(m.nsam));
        double w = static_cast<double>(S) / a1;
        auto a2 = summstats_aux::b_sub_n(static_cast<std::uint32_t>(m.nsam));
        auto dn = static_cast<double>(m.nsam);
        double b1 = (dn + 1.0) / (3.0 * (dn - 1.0));
        double b2
            = (2.0 * (std::pow(dn, 2.0) + dn + 3.0)) / (9.0 * dn * (dn - 1.0));
        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (dn + 2.0) / (a1 * dn) + a2 / std::pow(a1, 2.0);
        double e1 = c1 / a1;
        double e2 = c2 / (std::pow(a1, 2.0) + a2);
        double denominator = std::pow((e1 * S + e2 * S * (S - 1.0)), 0.5);
        return (pi - w) / denominator;
    }
} // namespace Sequence
