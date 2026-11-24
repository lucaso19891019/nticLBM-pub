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

// Negative velocity components select positive storage anchors.
// Slot zero means self; stored slots follow increasing x=1,y=2,z=4 masks.
template<LatticeType Type>
struct Twist
{
    using Lattice=Model<Type>;
    static constexpr std::size_t dimension=Type==LatticeType::D2Q9 ? 2 : 3;
    static constexpr unsigned mask(std::size_t q)
    {
        unsigned value=(Lattice::ex[q]<0 ? 1U : 0U)|(Lattice::ey[q]<0 ? 2U : 0U);
        if constexpr(dimension==3) value|=Lattice::ez[q]<0 ? 4U : 0U;
        return value;
    }
    inline static constexpr auto masks=[] {
        std::array<bool,8> result{};
        for(std::size_t q=0;q<Lattice::nStencils;++q) result[mask(q)]=true;
        return result;
    }();
    inline static constexpr std::size_t nAnchors=[] {
        std::size_t count=0;
        for(unsigned m=1;m<(1U<<dimension);++m) if(masks[m]) ++count;
        return count;
    }();
    inline static constexpr auto anchorMasks=[] {
        std::array<unsigned,nAnchors> result{}; std::size_t index=0;
        for(unsigned m=1;m<(1U<<dimension);++m) if(masks[m]) result[index++]=m;
        return result;
    }();
    inline static constexpr auto directionAnchor=[] {
        std::array<std::size_t,Lattice::nStencils> result{};
        for(std::size_t q=0;q<Lattice::nStencils;++q)
            for(std::size_t a=0;a<nAnchors;++a)
                if(mask(q)==anchorMasks[a]) result[q]=a+1;
        return result;
    }();
    inline static constexpr auto opposite=[] {
        std::array<std::size_t,Lattice::nStencils> result{};
        for(std::size_t q=0;q<Lattice::nStencils;++q)
        {
            result[q]=Lattice::nStencils;
            for(std::size_t r=0;r<Lattice::nStencils;++r)
                if(Lattice::ex[q]==-Lattice::ex[r] && Lattice::ey[q]==-Lattice::ey[r])
                {
                    if constexpr(dimension==3) { if(Lattice::ez[q]!=-Lattice::ez[r]) continue; }
                    result[q]=r; break;
                }
        }
        return result;
    }();
};

} // namespace ntic::lbm::lattice