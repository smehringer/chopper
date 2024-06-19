#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <cstddef>     // for size_t
#include <sstream>     // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>      // for allocator, string
#include <string_view> // for operator<<
#include <vector>      // for vector

#include <chopper/lsh.hpp>

TEST(Cluster_test, ctor_from_id)
{
    size_t user_bin_idx{5};
    chopper::Cluster const cluster{user_bin_idx};

    EXPECT_EQ(cluster.id(), user_bin_idx);
    EXPECT_FALSE(cluster.empty());
    EXPECT_EQ(cluster.size(), 1u);
    ASSERT_EQ(cluster.contained_user_bins().size(), 1u);
    EXPECT_EQ(cluster.contained_user_bins()[0], user_bin_idx);
    EXPECT_TRUE(cluster.is_valid(user_bin_idx));
}

TEST(Cluster_test, move_to)
{
    size_t user_bin_idx1{5};
    size_t user_bin_idx2{7};
    chopper::Cluster cluster1{user_bin_idx1};
    chopper::Cluster cluster2{user_bin_idx2};

    EXPECT_TRUE(cluster1.is_valid(user_bin_idx1));
    EXPECT_TRUE(cluster2.is_valid(user_bin_idx2));

    cluster2.move_to(cluster1);

    // cluster1 now contains user bins 5 and 7
    EXPECT_EQ(cluster1.size(), 2u);
    ASSERT_EQ(cluster1.contained_user_bins().size(), 2u);
    EXPECT_EQ(cluster1.contained_user_bins()[0], user_bin_idx1);
    EXPECT_EQ(cluster1.contained_user_bins()[1], user_bin_idx2);

    // cluster 2 is empty
    EXPECT_TRUE(cluster2.has_been_moved());
    EXPECT_TRUE(cluster2.empty());
    EXPECT_EQ(cluster2.size(), 0u);
    EXPECT_EQ(cluster2.contained_user_bins().size(), 0u);
    EXPECT_EQ(cluster2.moved_to_cluster_id(), cluster1.id());

    // both should still be valid
    EXPECT_TRUE(cluster1.is_valid(user_bin_idx1));
    EXPECT_TRUE(cluster2.is_valid(user_bin_idx2));
}

TEST(Multicluster_test, ctor_from_cluster)
{
    chopper::Cluster const cluster1{5};
    chopper::MultiCluster const multi_cluster1{cluster1};

    EXPECT_EQ(multi_cluster1.id(), cluster1.id());
    EXPECT_FALSE(multi_cluster1.empty());
    EXPECT_EQ(multi_cluster1.size(), 1u);
    ASSERT_EQ(multi_cluster1.contained_user_bins().size(), 1u);
    ASSERT_EQ(multi_cluster1.contained_user_bins()[0].size(), 1u);
    EXPECT_EQ(multi_cluster1.contained_user_bins()[0][0], 5u);
    EXPECT_TRUE(multi_cluster1.is_valid(cluster1.id()));
}

TEST(Multicluster_test, ctor_from_moved_cluster)
{
    chopper::Cluster cluster1{5};
    chopper::Cluster cluster2{7};
    cluster2.move_to(cluster1);

    chopper::MultiCluster const multi_cluster2{cluster2};

    EXPECT_EQ(multi_cluster2.id(), cluster2.id());
    EXPECT_TRUE(multi_cluster2.empty());
    EXPECT_TRUE(multi_cluster2.has_been_moved());
    EXPECT_EQ(multi_cluster2.moved_to_cluster_id(), cluster2.moved_to_cluster_id());
    EXPECT_EQ(multi_cluster2.size(), 0u);
    ASSERT_EQ(multi_cluster2.contained_user_bins().size(), 0u);
    EXPECT_TRUE(multi_cluster2.is_valid(cluster2.id()));
}