//! \file VariantMatrixTest.cc @brief Tests for Sequence/VariantMatrix.hpp
#define BOOST_TEST_MODULE VariantMatrixTest

#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <boost/test/included/unit_test.hpp>
#include <algorithm>
#include <numeric> //for std::iota
#include <iterator>

struct vmatrix_fixture
{
    std::vector<std::int8_t> input_data;
    std::vector<double> input_pos;
    Sequence::VariantMatrix m, m2;

    vmatrix_fixture()
        : input_data(make_input_data()), input_pos(make_intput_pos()),
          m(input_data, input_pos), m2(input_data, input_pos)
    {
        // The two VariantMatrix objects
        // have same data, but different internal
        // dimensions
        std::swap(m2.nsites, m2.nsam);
        m2.positions.resize(m2.nsites);
        std::iota(std::begin(m2.positions), std::end(m2.positions), 0.);
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

BOOST_FIXTURE_TEST_CASE(test_bad_row_swap, vmatrix_fixture)
{
    auto a = Sequence::get_RowView(m, 0);
    auto b = Sequence::get_RowView(m2, 0);
    BOOST_REQUIRE_THROW(swap(a, b), std::invalid_argument);
}

BOOST_FIXTURE_TEST_CASE(test_bad_column_swap, vmatrix_fixture)
{
    auto a = Sequence::get_ColView(m, 0);
    auto b = Sequence::get_ColView(m2, 0);
    BOOST_REQUIRE_THROW(swap(a, b), std::invalid_argument);
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

BOOST_FIXTURE_TEST_CASE(test_row_view_iterators, vmatrix_fixture)
{
    for (std::size_t i = 0; i < m.nsites; ++i)
        {
            auto row = Sequence::get_RowView(m, i);
            BOOST_REQUIRE_EQUAL(std::distance(row.begin(), row.end()), m.nsam);
            BOOST_REQUIRE_EQUAL(std::distance(row.cbegin(), row.cend()),
                                m.nsam);
        }
}

BOOST_FIXTURE_TEST_CASE(test_column_views, vmatrix_fixture)
{
    for (std::size_t i = 0; i < m.nsam; ++i)
        {
            auto col = Sequence::get_ColView(m, i);
            std::int8_t state = (i % 2 == 0) ? 1 : 0;
            BOOST_REQUIRE_EQUAL(
                std::count(std::begin(col), std::end(col), !state), 0);
            BOOST_REQUIRE_EQUAL(std::count(col.rbegin(), col.rend(), !state),
                                0);
            BOOST_REQUIRE_EQUAL(std::count(col.cbegin(), col.cend(), !state),
                                0);
            BOOST_REQUIRE_EQUAL(std::count(col.crbegin(), col.crend(), !state),
                                0);

            BOOST_REQUIRE_EQUAL(std::distance(std::begin(col), std::end(col)),
                                m.nsites);
            BOOST_REQUIRE_EQUAL(std::distance(col.rbegin(), col.rend()),
                                m.nsites);

            // Check that iterators and reverse iterators have the expected
            // relationships:
            auto fwd = col.begin();
            auto rev = col.rbegin();
            for (; rev < col.rend(); ++rev)
                {
                    auto rf = std::distance(fwd, rev.base());
                    auto rb = std::distance(rev, col.rend());
                    BOOST_REQUIRE_EQUAL(rf, rb);
                }

            auto cfwd = col.cbegin();
            auto crev = col.crbegin();
            for (; crev < col.crend(); ++crev)
                {
                    auto rf = std::distance(cfwd, crev.base());
                    auto rb = std::distance(crev, col.crend());
                    BOOST_REQUIRE_EQUAL(rf, rb);
                }
        }
}

BOOST_FIXTURE_TEST_CASE(tesl_col_view_iterator_increment, vmatrix_fixture)
{
    auto x = Sequence::get_ConstColView(m, 0);
    auto b = x.begin();
    unsigned num_increments = 0;
    while (b < x.end())
        {
            b = b + 2;
            ++num_increments;
        }
    BOOST_REQUIRE_EQUAL(num_increments, 3);
}

BOOST_FIXTURE_TEST_CASE(test_column_view_invalid_compare, vmatrix_fixture)
{
    auto c0 = Sequence::get_ConstColView(m, 0);
    auto c1 = Sequence::get_ConstColView(m, 1);
    BOOST_REQUIRE_NO_THROW(std::distance(c0.begin(), c0.end()));
    BOOST_REQUIRE_THROW(std::distance(c0.begin(), c1.begin()),
                        std::invalid_argument);
}

// The remaining tests apply STL algorithms to column iterators,
// which is a good stress test.  We've already done count above.

BOOST_FIXTURE_TEST_CASE(test_accumulate, vmatrix_fixture)
{
    auto c = Sequence::get_ConstColView(m, 0);
    int sum = static_cast<int>(std::accumulate(c.cbegin(), c.cend(), 0));
    BOOST_REQUIRE_EQUAL(sum, static_cast<int>(m.nsites));
    c = Sequence::get_ConstColView(m, 1);
    sum = static_cast<int>(std::accumulate(c.cbegin(), c.cend(), 0));
    BOOST_REQUIRE_EQUAL(sum, 0);
}