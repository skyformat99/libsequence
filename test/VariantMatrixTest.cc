//! \file ComparisonsTests.cc @brief Tests for Sequence/Comparisons.hpp
#define BOOST_TEST_MODULE VariantMatrixTest

#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <boost/test/included/unit_test.hpp>
#include <algorithm>
#include <numeric>  //for std::iota
#include <iterator>

struct vmatrix_fixture
{
    std::vector<std::int8_t> input_data;
    std::vector<double> input_pos;
    Sequence::VariantMatrix m;

    vmatrix_fixture()
        : input_data(make_input_data()), input_pos(make_intput_pos()),
          m(input_data, input_pos)
    {
    }

    std::vector<std::int8_t>
    make_input_data()
    {
        int nsam = 20;
        int nsites = 5;
        std::vector<std::int8_t> rv;
        for (int i = 0; i < nsites; ++i)
            {
                for (int j = 0; j < nsam; ++j)
                    {
                        std::int8_t state = (j % 2 == 0.) ? 1 : 0;
                        rv.push_back(state);
                    }
            }
        return rv;
    }

    std::vector<double>
    make_intput_pos()
    {
        std::vector<double> rv;
        rv.resize(5);
        std::iota(rv.begin(), rv.end(), 0.);
        return rv;
    }
};

BOOST_FIXTURE_TEST_CASE(test_construction, vmatrix_fixture)
{
    Sequence::VariantMatrix m(std::vector<std::int8_t>(100, 1),
                              std::vector<double>(5, 0.0));
    BOOST_REQUIRE_EQUAL(m.nsites, 5);
    BOOST_REQUIRE_EQUAL(m.nsam, 20);
}

BOOST_FIXTURE_TEST_CASE(test_range_exceptions, vmatrix_fixture)
{
    BOOST_REQUIRE_THROW(m.at(m.nsites + 1, 0), std::out_of_range);
    BOOST_REQUIRE_THROW(m.at(0, m.nsam + 1), std::out_of_range);
}

BOOST_FIXTURE_TEST_CASE(test_iteration, vmatrix_fixture)
{
    for (std::size_t i = 0; i < m.nsam; ++i)
        {
            for (std::size_t j = 0; j < m.nsites; ++j)
                {
                    auto x = m.get(j, i);
                    std::int8_t ex = (i % 2 == 0.) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(x),
                                        static_cast<int>(ex));
                }
        }
}

BOOST_FIXTURE_TEST_CASE(test_row_views, vmatrix_fixture)
{
    for (std::size_t i = 0; i < m.nsites; ++i)
        {
            auto x = Sequence::get_RowView(m, i);
            for (auto j = x.begin(); j != x.end(); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.begin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (auto j = x.cbegin(); j != x.cend(); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.cbegin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (auto j = std::begin(x); j != std::end(x); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.begin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (std::size_t j = 0; j < x.size(); ++j)
                {
                    std::int8_t ex = (j % 2 == 0.) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(x[j]),
                                        static_cast<int>(ex));
                }
            std::size_t j = 0;
            for (auto xj : x)
                {
                    std::int8_t ex = (j % 2 == 0.) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(xj),
                                        static_cast<int>(ex));
                    ++j;
                }
        }
}

BOOST_FIXTURE_TEST_CASE(test_const_row_views, vmatrix_fixture)
{
    for (std::size_t i = 0; i < m.nsites; ++i)
        {
            auto x = Sequence::get_ConstRowView(m, i);

            for (auto j = x.begin(); j != x.end(); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.begin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (auto j = x.cbegin(); j != x.cend(); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.cbegin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (auto j = std::begin(x); j != std::end(x); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.begin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (std::size_t j = 0; j < x.size(); ++j)
                {
                    std::int8_t ex = (j % 2 == 0.) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(x[j]),
                                        static_cast<int>(ex));
                }
            std::size_t j = 0;
            for (auto xj : x)
                {
                    std::int8_t ex = (j % 2 == 0.) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(xj),
                                        static_cast<int>(ex));
                    ++j;
                }
        }
}

BOOST_FIXTURE_TEST_CASE(test_row_view_exceptions, vmatrix_fixture)
{
    BOOST_REQUIRE_THROW(Sequence::get_RowView(m, m.nsites + 1),
                        std::exception);
    BOOST_REQUIRE_THROW(Sequence::get_RowView(m, m.nsites + 1),
                        std::out_of_range);

    auto r = Sequence::get_RowView(m, 0);
    BOOST_REQUIRE_THROW(r.at(m.nsam + 1), std::exception);
    BOOST_REQUIRE_THROW(r.at(m.nsam + 1), std::out_of_range);
}

BOOST_FIXTURE_TEST_CASE(test_const_row_view_exceptions, vmatrix_fixture)
{
    BOOST_REQUIRE_THROW(Sequence::get_ConstRowView(m, m.nsites + 1),
                        std::exception);
    BOOST_REQUIRE_THROW(Sequence::get_ConstRowView(m, m.nsites + 1),
                        std::out_of_range);

    auto r = Sequence::get_ConstRowView(m, 0);
    BOOST_REQUIRE_THROW(r.at(m.nsam + 1), std::exception);
    BOOST_REQUIRE_THROW(r.at(m.nsam + 1), std::out_of_range);
}
