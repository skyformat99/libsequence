//! \file ParallelTransformTest.cc @brief Tests for Sequence/parallel/transform.hpp
#include <boost/test/unit_test.hpp>
#include <Sequence/parallel/transform.hpp>
#include <vector>
#include <list>

BOOST_AUTO_TEST_SUITE(ParallelTransformTest)

BOOST_AUTO_TEST_CASE(test_vector)
{
    std::vector<int> x(100000, 3);
    std::vector<int> y(x.size());
    auto f = [](const int i) { return i; };
    auto e = Sequence::parallel::transform(x.begin(), x.end(), y.begin(), f);
    BOOST_REQUIRE_EQUAL(e == y.end(), true);
    BOOST_REQUIRE_EQUAL(
        std::all_of(y.begin(), y.end(), [](const int i) { return i == 3; }),
        true);
}

BOOST_AUTO_TEST_CASE(test_list)
{
    std::list<int> x(100000, 3);
    std::list<int> y(x.size());
    auto f = [](const int i) { return i; };
    auto e = Sequence::parallel::transform(x.begin(), x.end(), y.begin(), f);
    BOOST_REQUIRE_EQUAL(e == y.end(), true);
    BOOST_REQUIRE_EQUAL(
        std::all_of(y.begin(), y.end(), [](const int i) { return i == 3; }),
        true);
}

BOOST_AUTO_TEST_SUITE_END()
