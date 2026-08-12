#include "model_alias_registry.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>

using ntic::lbm::config::ModelAliasRegistry;

int main()
{
    const ModelAliasRegistry r;

    NTIC_EXPECT_EQ(r.resolveLatticeModel("d2q9"), std::string("D2Q9"));
    NTIC_EXPECT_EQ(r.resolveLatticeModel("D3q15"), std::string("D3Q15"));
    NTIC_EXPECT_EQ(r.resolveLatticeModel("d3q19"), std::string("D3Q19"));
    NTIC_EXPECT_EQ(r.resolveLatticeModel("d3q27"), std::string("D3Q27"));

    NTIC_EXPECT_EQ(r.resolveCollisionModel("bgk"), std::string("BGK"));
    NTIC_EXPECT_EQ(r.resolveCollisionModel("mrt"), std::string("MRT"));
    NTIC_EXPECT_EQ(r.resolveCollisionModel("trt"), std::string("TRT"));
    NTIC_EXPECT_EQ(r.resolveCollisionModel("central moment"), std::string("Cascaded"));
    NTIC_EXPECT_EQ(r.resolveCollisionModel("CM"), std::string("Cascaded"));
    NTIC_EXPECT_EQ(r.resolveCollisionModel("cumulant"), std::string("Cumulant"));

    NTIC_EXPECT_EQ(r.resolveWallBoundary("ibb"), std::string("InterpolatedBB"));
    NTIC_EXPECT_EQ(r.resolveWallBoundary("immersed"), std::string("IBM"));

    NTIC_EXPECT_EQ(r.resolveInletBoundary("default"), std::string("Default"));
    NTIC_EXPECT_EQ(r.resolveInletBoundary("zh-vel"), std::string("ZouHeVelocity"));
    NTIC_EXPECT_EQ(r.resolveInletBoundary("zh pre"), std::string("ZouHePressure"));

    NTIC_EXPECT_EQ(r.resolveOutletBoundary("standard"), std::string("Default"));
    NTIC_EXPECT_EQ(r.resolveOutletBoundary("zou-he velocity"), std::string("ZouHeVelocity"));
    NTIC_EXPECT_EQ(r.resolveOutletBoundary("zh-pressure"), std::string("ZouHePressure"));

    NTIC_EXPECT_TRUE(!r.containsCollisionModel("entropic"));
    NTIC_EXPECT_THROW(r.resolveCollisionModel("entropic"), std::runtime_error);
    NTIC_EXPECT_THROW(r.resolveWallBoundary("bounce back"), std::runtime_error);

    ntic::lbm::test::printSummary();
    return 0;
}
