#pragma once

#include <array>
#include <cstddef>


namespace ntic::lbm::lattice
{

enum class LatticeType
{
    D2Q9,
    D3Q15,
    D3Q19,
    D3Q27
};


template <LatticeType Type>
struct Model;


//=============================================================================
// D2Q9
//=============================================================================

template <>
struct Model<LatticeType::D2Q9>
{
    static constexpr std::size_t nStencils = 9;


    inline static constexpr std::array<int, nStencils> ex =
    {
         0,  1,  1,  0, -1, -1, -1,  0,  1
    };

    inline static constexpr std::array<int, nStencils> ey =
    {
         0,  0,  1,  1,  1,  0, -1, -1, -1
    };

    inline static constexpr std::array<double, nStencils> weight =
    {
        4.0 / 9.0,
        1.0 / 9.0, 1.0 / 36.0, 1.0 / 9.0, 1.0 / 36.0,
        1.0 / 9.0, 1.0 / 36.0, 1.0 / 9.0, 1.0 / 36.0
    };
};


//=============================================================================
// D3Q15
//=============================================================================

template <>
struct Model<LatticeType::D3Q15>
{
    static constexpr std::size_t nStencils = 15;


    inline static constexpr std::array<int, nStencils> ex =
    {
         0,  0,  1,  0,  1, -1, -1,  1,
         0, -1,  0, -1,  1,  1, -1
    };

    inline static constexpr std::array<int, nStencils> ey =
    {
         0,  0,  0,  1,  1,  1, -1, -1,
        -1,  0,  0, -1, -1,  1,  1
    };

    inline static constexpr std::array<int, nStencils> ez =
    {
         0,  1,  0,  0,  1,  1,  1,  1,
         0,  0, -1, -1, -1, -1, -1
    };

    inline static constexpr std::array<double, nStencils> weight =
    {
        2.0 / 9.0,

        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,

        1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0,

        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,

        1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0
    };
};


//=============================================================================
// D3Q19
//=============================================================================

template <>
struct Model<LatticeType::D3Q19>
{
    static constexpr std::size_t nStencils = 19;


    inline static constexpr std::array<int, nStencils> ex =
    {
         0,  0,  1,  0,  1,  0, -1,  0,  1, -1,
         0, -1,  0, -1,  0,  1,  0, -1,  1
    };

    inline static constexpr std::array<int, nStencils> ey =
    {
         0,  0,  0,  1,  0,  1,  0, -1,  1,  1,
         0,  0, -1,  0, -1,  0,  1, -1, -1
    };

    inline static constexpr std::array<int, nStencils> ez =
    {
         0,  1,  0,  0,  1,  1,  1,  1,  0,  0,
        -1,  0,  0, -1, -1, -1, -1,  0,  0
    };

    inline static constexpr std::array<double, nStencils> weight =
    {
        1.0 / 3.0,

        1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,

        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,

        1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,

        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };
};


//=============================================================================
// D3Q27
//=============================================================================

template <>
struct Model<LatticeType::D3Q27>
{
    static constexpr std::size_t nStencils = 27;


    inline static constexpr std::array<int, nStencils> ex =
    {
         0,  0,  1,  0,  1,  0, -1,  0,  1, -1,
         1, -1, -1,  1,  0, -1,  0, -1,  0,  1,
         0, -1,  1, -1,  1,  1, -1
    };

    inline static constexpr std::array<int, nStencils> ey =
    {
         0,  0,  0,  1,  0,  1,  0, -1,  1,  1,
         1,  1, -1, -1,  0,  0, -1,  0, -1,  0,
         1, -1, -1, -1, -1,  1,  1
    };

    inline static constexpr std::array<int, nStencils> ez =
    {
         0,  1,  0,  0,  1,  1,  1,  1,  0,  0,
         1,  1,  1,  1, -1,  0,  0, -1, -1, -1,
        -1,  0,  0, -1, -1, -1, -1
    };

    inline static constexpr std::array<double, nStencils> weight =
    {
        8.0 / 27.0,

        2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0,

        1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
        1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,

        1.0 / 216.0, 1.0 / 216.0,
        1.0 / 216.0, 1.0 / 216.0,

        2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0,

        1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
        1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,

        1.0 / 216.0, 1.0 / 216.0,
        1.0 / 216.0, 1.0 / 216.0
    };
};

} // namespace ntic::lbm::lattice