#include "neighbor_linking.hpp"
#include "test_multilevel_checks.hpp"

#include <iostream>
#include <set>
#include <stdexcept>

namespace geo=ntic::lbm::geometry;
namespace lat=ntic::lbm::lattice;

namespace
{

template<lat::LatticeType Type>
constexpr bool validTraits()
{
    using T=lat::Twist<Type>; using L=lat::Model<Type>;
    for(std::size_t q=0;q<L::nStencils;++q)
    {
        const auto opposite=T::opposite[q];
        if(opposite>=L::nStencils || T::opposite[opposite]!=q) return false;
        const auto slot=T::directionAnchor[q];
        if(slot>T::nAnchors || (slot==0 ? 0U : T::anchorMasks[slot-1])!=T::mask(q)) return false;
        for(std::size_t d=0;d<T::dimension;++d)
        {
            const int e=d==0 ? L::ex[q] : L::ey[q];
            int velocity=e;
            if constexpr(T::dimension==3) if(d==2) velocity=L::ez[q];
            if(velocity!=static_cast<int>((T::mask(opposite)>>d)&1U)-static_cast<int>((T::mask(q)>>d)&1U)) return false;
        }
    }
    return true;
}
static_assert(lat::Twist<lat::LatticeType::D2Q9>::nAnchors==3);
static_assert(lat::Twist<lat::LatticeType::D3Q15>::nAnchors==7);
static_assert(lat::Twist<lat::LatticeType::D3Q19>::nAnchors==6);
static_assert(lat::Twist<lat::LatticeType::D3Q27>::nAnchors==7);
static_assert(validTraits<lat::LatticeType::D2Q9>() && validTraits<lat::LatticeType::D3Q15>() &&
              validTraits<lat::LatticeType::D3Q19>() && validTraits<lat::LatticeType::D3Q27>());
// D2Q9 q1=(+x,0), q5=(-x,0), q8=(+x,-y).
static_assert(lat::Twist<lat::LatticeType::D2Q9>::directionAnchor[1]==0);
static_assert(lat::Twist<lat::LatticeType::D2Q9>::directionAnchor[5]==1);
static_assert(lat::Twist<lat::LatticeType::D2Q9>::directionAnchor[8]==2);

template<std::size_t D,lat::LatticeType Type>
void checkLinks(const geo::MultilevelGrid<D,Type>& grid,const geo::NeighborLinks<D,Type>& links)
{
    using T=lat::Twist<Type>;
    using Key=typename geo::NeighborLinks<D,Type>::Ghost;
    requireMultilevel(links.originalCount==grid.cells().size(),"Linking changed the original ID prefix");
    std::set<Key> used;
    for(std::size_t c=0;c<grid.cells().size();++c)
        for(std::size_t a=0;a<T::nAnchors;++a)
        {
            Key expected; expected.level=grid.level(c);
            const auto coordinate=grid.coordinate(c);
            for(std::size_t d=0;d<D;++d) expected.coordinate[d]=coordinate[d]+((T::anchorMasks[a]>>d)&1U);
            geo::CellIndex found=geo::invalidIndex;
            // Independent small-fixture oracle: no coveringNode or auxiliary-family lookup.
            for(std::size_t r=0;r<grid.cells().size();++r)
                if(grid.level(r)==expected.level)
                {
                    const auto p=grid.coordinate(r); bool equal=true;
                    for(std::size_t d=0;d<D;++d) equal=equal && p[d]==expected.coordinate[d];
                    if(equal) { requireMultilevel(found==geo::invalidIndex,"Duplicate stored position"); found=r; }
                }
            const auto actual=links.anchors[c][a];
            if(found!=geo::invalidIndex) requireMultilevel(actual==found,"Same-level existing slot was not reused");
            else
            {
                requireMultilevel(links.isGhost(actual),"Missing position has no Ghost endpoint");
                requireMultilevel(links.ghosts[actual-links.originalCount]==expected,"Incorrect Ghost position");
                used.insert(expected);
            }
        }
    requireMultilevel(used.size()==links.ghosts.size(),"Ghosts are duplicated or recursively expanded");
    for(std::size_t level=0;level<links.ghostLevels.size();++level)
        for(auto id=links.ghostLevels[level].begin;id<links.ghostLevels[level].end;++id)
            requireMultilevel(links.ghosts.at(id-links.originalCount).level==level,"Incorrect Ghost level range");
    for(std::size_t c=0;c<grid.cells().size();++c)
        for(std::size_t q=0;q<lat::Model<Type>::nStencils;++q)
            requireMultilevel(links.anchor(c,q)==(T::directionAnchor[q]==0 ? c : links.anchors[c][T::directionAnchor[q]-1]),
                              "Direction-to-storage lookup differs from its trait");
}

template<lat::LatticeType Type>
void uniform3D(std::size_t ghostCount)
{
    const geo::Box shape({0,0,0},{2,3,4});
    geo::MultilevelGrid<3,Type> grid(shape.bounds,1); grid.generate(shape);
    const auto links=geo::linkNeighbors(grid); checkLinks(grid,links);
    requireMultilevel(grid.cells().size()==24 && links.ghosts.size()==ghostCount,"Incorrect uniform 3D Ghost count");
}

void uniform2D()
{
    const geo::Rectangle shape({0,0},{2,3});
    geo::MultilevelGrid<2> grid(shape.bounds,1);
    bool rejected=false;
    try { (void)geo::linkNeighbors(grid); } catch(const std::logic_error&) { rejected=true; }
    requireMultilevel(rejected,"Unfinalized grid was accepted");
    grid.generate(shape); const auto links=geo::linkNeighbors(grid); checkLinks(grid,links);
    requireMultilevel(grid.cells().size()==6 && links.ghosts.size()==6,"Incorrect uniform 2D Ghost count");
    rejected=false;
    try { (void)links.anchor(links.originalCount,0); } catch(const std::out_of_range&) { rejected=true; }
    requireMultilevel(rejected,"Ghost endpoint unexpectedly has its own neighbor table");
}

void adaptive3D()
{
    const geo::Box shape({0,0,0},{4,4,4});
    const geo::Sphere selection({2,2,2},0.4);
    geo::MultilevelOptions options; options.maxLevel=1;
    geo::MultilevelGrid<3,lat::LatticeType::D3Q27> grid(shape.bounds,1,options);
    grid.generate(shape,{geo::RefinementRule<3>::volume(geo::volumeRegion(selection),1)});
    const auto links=geo::linkNeighbors(grid); checkLinks(grid,links);
    bool auxiliaryTarget=false,auxiliarySource=false,upperZ=false;
    for(std::size_t c=0;c<grid.cells().size();++c)
        for(auto id:links.anchors[c]) if(id<links.originalCount)
        {
            auxiliaryTarget=auxiliaryTarget || grid.cells().isAuxiliary(id);
            auxiliarySource=auxiliarySource || (grid.cells().isAuxiliary(c) && !grid.cells().isAuxiliary(id));
            upperZ=upperZ || (grid.cells().isAuxiliary(id) && (grid.cells().auxiliaryChild[id]&4U)!=0);
        }
    requireMultilevel(auxiliaryTarget && auxiliarySource && upperZ,"3D auxiliary reuse does not exercise the z child bit");
}

void sparseAndMultilevel()
{
    const auto shape=geo::Rectangle({0,0},{8,8})-geo::Circle({4,4},1.3);
    const geo::Rectangle selection({2.3,0},{5.7,8});
    geo::MultilevelOptions options; options.maxLevel=1;
    geo::MultilevelGrid<2> grid(shape.bounds,1,options);
    grid.generate(shape,{geo::RefinementRule<2>::volume(geo::volumeRegion(selection),1)});
    const auto oldCells=grid.cells().sourceNode,oldBoundary=grid.boundaries().cells;
    omp_set_num_threads(1); const auto one=geo::linkNeighbors(grid);
    omp_set_num_threads(4); const auto four=geo::linkNeighbors(grid);
    checkLinks(grid,four);
    geo::NeighborLinks<2,lat::LatticeType::D2Q9> nested;
    std::exception_ptr nestedError;
#pragma omp parallel num_threads(4)
    {
#pragma omp single
        {
            try { nested=geo::linkNeighbors(grid); } catch(...) { nestedError=std::current_exception(); }
        }
    }
    if(nestedError) std::rethrow_exception(nestedError);
    requireMultilevel(nested.anchors==four.anchors && nested.ghosts==four.ghosts,"Nested linking differs from standalone");
    requireMultilevel(one.anchors==four.anchors && one.ghosts==four.ghosts,"Linking depends on thread count");
    requireMultilevel(oldCells==grid.cells().sourceNode && oldBoundary==grid.boundaries().cells,"Linking mutated geometry indices");
    bool auxiliaryTarget=false,auxiliarySource=false,hole=false;
    for(std::size_t c=0;c<grid.cells().size();++c)
        for(auto id:four.anchors[c]) if(id<four.originalCount)
        {
            auxiliaryTarget=auxiliaryTarget || grid.cells().isAuxiliary(id);
            auxiliarySource=auxiliarySource || (grid.cells().isAuxiliary(c) && !grid.cells().isAuxiliary(id));
        }
    for(const auto& ghost:four.ghosts)
    {
        const auto h=grid.spacing(ghost.level);
        const double x=(ghost.coordinate[0]+0.5)*h,y=(ghost.coordinate[1]+0.5)*h;
        hole=hole || ((x-4)*(x-4)+(y-4)*(y-4)<1.3*1.3);
    }
    requireMultilevel(auxiliaryTarget && auxiliarySource,"Real/auxiliary reuse was not exercised");
    requireMultilevel(hole,"Dry-hole endpoint was not exercised");
    grid.refine(0);
    bool rejected=false;
    try { (void)geo::linkNeighbors(grid); } catch(const std::logic_error&) { rejected=true; }
    requireMultilevel(rejected,"Invalidated grid was accepted");
    grid.finalize(shape); checkLinks(grid,geo::linkNeighbors(grid));
    grid.derefine(0); grid.finalize(shape); checkLinks(grid,geo::linkNeighbors(grid));
    const geo::Circle absent({20,20},1);
    grid.finalize(absent);
    const auto empty=geo::linkNeighbors(grid);
    requireMultilevel(empty.originalCount==0 && empty.size()==0 && empty.anchors.empty(),"Dry grid has storage endpoints");
}

} // namespace

int main()
{
    try
    {
        omp_set_dynamic(0); omp_set_num_threads(4);
        uniform2D();
        uniform3D<lat::LatticeType::D3Q15>(36);
        uniform3D<lat::LatticeType::D3Q19>(35);
        uniform3D<lat::LatticeType::D3Q27>(36);
        std::cout << "[PASS] twist_traits_and_uniform_ghosts\n";
        sparseAndMultilevel(); adaptive3D();
        std::cout << "[PASS] sparse_multilevel_neighbor_links\n";
        return 0;
    }
    catch(const std::exception& error) { std::cerr << "[FAIL] neighbors: " << error.what() << '\n'; return 1; }
}
