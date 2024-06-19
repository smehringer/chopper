#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <cstddef>     // for size_t
#include <sstream>     // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>      // for allocator, string
#include <string_view> // for operator<<
#include <vector>      // for vector

#include <chopper/lsh.hpp>

TEST(Multicluster_test, ctor_from_cluster)
{
    chopper::Cluster const cluster1{5};
    chopper::MultiCluster const multi_cluster1{cluster1};

    EXPECT_EQ(multi_cluster1.id(), 5u);
    EXPECT_FALSE(multi_cluster1.empty());
    EXPECT_EQ(multi_cluster1.size(), 1u);
    ASSERT_EQ(multi_cluster1.contained_user_bins().size(), 1u);
    ASSERT_EQ(multi_cluster1.contained_user_bins()[0].size(), 1u);
    EXPECT_EQ(multi_cluster1.contained_user_bins()[0][0], 5u);
    EXPECT_TRUE(multi_cluster1.is_valid(cluster1.id()));
}