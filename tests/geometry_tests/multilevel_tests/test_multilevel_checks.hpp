#pragma once

#include "multilevel_grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace geo=ntic::lbm::geometry;
namespace lat=ntic::lbm::lattice;

namespace ntic::lbm::geometry::test
{
inline void require(bool condition,const char* message)
{ if(!condition) throw std::runtime_error(message); }

namespace multilevel_check_detail
{
template<std::size_t D,lat::LatticeType Type>
void topology(const geo::MultilevelGrid<D,Type>& grid)
{
    const auto& cells=grid.cells(); const auto& nodes=grid.nodes();
    double volume=0,expected=1;
    for(std::size_t d=0; d<D; ++d) expected*=grid.domain().max[d]-grid.domain().min[d];
    std::set<std::pair<unsigned,std::array<std::uint32_t,D>>> positions;
    for(std::size_t n=0; n<nodes.size(); ++n)
    {
        const auto& node=nodes[n];
        if(!node.active || !node.leaf()) require(node.cell==geo::invalidIndex,"Covered/inactive parent has data storage");
        if(!node.active) continue;
        if(!node.leaf())
            for(std::size_t child=0; child<(std::size_t{1}<<D); ++child)
            {
                const auto& c=nodes[node.firstChild+child];
                require(c.parent==n && c.level==node.level+1 && c.active,"Invalid indexed parent/child relation");
                for(std::size_t d=0; d<D; ++d)
                    require(c.coordinate[d]==2*node.coordinate[d]+((child>>d)&1U),"Invalid child coordinate");
            }
        else
        {
            if(node.cell!=geo::invalidIndex)
                require(node.cell<cells.size() && cells.sourceNode[node.cell]==n,"Leaf has no unique data index");
            volume+=std::pow(grid.spacing(node.level),D);
        }
    }
    require(std::abs(volume-expected)<1e-10*expected,"Leaf volume is missing or double counted");
    for(std::size_t c=0; c<cells.size(); ++c)
    {
        require(positions.emplace(grid.level(c),grid.coordinate(c)).second,"Duplicate cell at the same level");
        if(cells.isAuxiliary(c))
        {
            const auto source=nodes[cells.sourceNode[c]].cell;
            require(source<cells.size() && (cells.role[source]&geo::CoarseInterface),"Auxiliary source is not a coarse interface leaf");
            require(grid.level(c)==grid.level(source)+1,"Auxiliary level mismatch");
            require(!cells.computes(c),"Auxiliary counted as independent solving volume");
        }
        else
        {
            require(cells.cellType[c]!=geo::CellType::Dry,"Dry real leaf occupies a data slot");
            for(const auto neighbour:grid.adjacentRealCells(c))
            {
                require(!cells.isAuxiliary(neighbour),"Partition adjacency contains auxiliary cell");
                require(std::abs(static_cast<int>(grid.level(c))-static_cast<int>(grid.level(neighbour)))<=1,"2:1 balance failed");
                const auto reverse=grid.adjacentRealCells(neighbour);
                require(std::binary_search(reverse.begin(),reverse.end(),c),"Nonreciprocal partition adjacency");
            }
        }
    }
    std::size_t end=0;
    for(unsigned l=0; l<grid.levels().size(); ++l)
    {
        const auto& level=grid.levels()[l]; require(level.all.begin==end,"Noncontiguous level storage");
        for(unsigned role=0; role<5; ++role)
        {
            const auto range=level.roles[role]; require(range.begin==end,"Noncontiguous role storage");
            for(auto c=range.begin; c<range.end; ++c)
                require(grid.level(c)==l && cells.role[c]==role,"Wrong role/level range");
            end=range.end;
        }
        require(end==level.all.end,"Level endpoint mismatch");
        for(unsigned role=0; role<4; ++role)
        {
            auto cursor=level.roles[role].begin;
            const std::array<geo::CellType,3> types{geo::CellType::Interior,geo::CellType::Boundary,geo::CellType::Dry};
            for(unsigned type=0; type<3; ++type)
            {
                const auto range=level.types[role][type]; require(range.begin==cursor,"Noncontiguous geometry blocks");
                for(auto c=range.begin; c<range.end; ++c) require(cells.cellType[c]==types[type],"Wrong geometry block");
                cursor=range.end;
            }
            require(cursor==level.roles[role].end,"Geometry block endpoint mismatch");
        }
    }
    require(end==cells.size(),"Cell pool is not fully covered by ranges");
    std::size_t auxiliaries=0;
    for(const auto& block:grid.interfaces().blocks)
    {
        require(block.coarse<cells.size() && (cells.role[block.coarse]&geo::CoarseInterface),"Invalid interface coarse index");
        require(block.auxiliary.size()==(std::size_t{1}<<D),"Incomplete auxiliary family");
        auxiliaries+=block.auxiliary.size();
        for(auto c=block.auxiliary.begin; c<block.auxiliary.end; ++c)
        {
            require(cells.isAuxiliary(c) && cells.sourceNode[c]==cells.sourceNode[block.coarse],"Auxiliary group lost contiguity");
            require(cells.auxiliaryChild[c]==c-block.auxiliary.begin,"Auxiliary child order changed");
        }
        require(block.fine.size()>0,"Coarse exchange block has no adjacent fine cells");
        for(auto i=block.fine.begin; i<block.fine.end; ++i)
        {
            const auto c=grid.interfaces().fineCells[i];
            require(!cells.isAuxiliary(c) && (cells.role[c]&geo::FineInterface),"Invalid fine interface index");
            require(grid.level(c)==grid.level(block.coarse)+1,"Interface skips a level");
        }
    }
    std::size_t count=0; for(auto role:cells.role) count+=role==geo::Auxiliary;
    require(auxiliaries==count,"Auxiliary cells missing from interface table");
}

template<std::size_t D,lat::LatticeType Type,typename Geometry>
void compareGeometry(const geo::MultilevelGrid<D,Type>& grid,const Geometry& geometry)
{
    using Analysis=std::conditional_t<D==2,geo::GeometryAnalysis2D<Type>,geo::GeometryAnalysis3D<Type>>;
    using Lattice=typename Analysis::Lattice;
    constexpr auto links=Lattice::nStencils-1;
    using Key=std::pair<std::size_t,std::size_t>;
    const auto& cells=grid.cells(); const auto& boundary=grid.boundaries();
    require(boundary.offsets.size()==boundary.nBoundaries+1 && boundary.offsets.front()==0 &&
            boundary.offsets.back()==boundary.cells.size(),"Invalid multilevel boundary offsets");
    require(boundary.q.size()==boundary.cells.size()*links,"Invalid multilevel Q shape");
    std::map<Key,std::array<double,links>> actual;
    std::set<std::pair<std::size_t,std::size_t>> owners;
    for(std::size_t face=0; face<boundary.nBoundaries; ++face)
        for(auto entry=boundary.offsets[face]; entry<boundary.offsets[face+1]; ++entry)
        {
            const auto cell=boundary.cells[entry];
            require(cells.cellType[cell]==geo::CellType::Boundary,"Q references a non-boundary cell");
            std::array<double,links> q{};
            for(std::size_t link=0; link<links; ++link)
            {
                q[link]=boundary.q[entry*links+link];
                require(q[link]==-1 || (q[link]>0 && q[link]<=1),"Invalid Q");
                if(q[link]>0) require(owners.emplace(cell,link).second,"Duplicate physical cell/link Q owner");
            }
            require(actual.emplace(Key{face,cell},q).second,"Duplicate Face/cell entry");
        }
    unsigned treeLevels=1;
    for(const auto& node:grid.nodes()) if(node.leaf()) treeLevels=std::max(treeLevels,node.level+1);
    for(unsigned l=0; l<treeLevels; ++l)
    {
        Analysis reference(grid.spacing(l)); reference.domainBounds=grid.domain();
        reference.gridFrame=geo::AnalysisGridFrame<D>{grid.domain().min,{}};
        geometry.analysis(reference);
        reference.validateBoundaryCells();
        for(const auto& node:grid.nodes())
        {
            if(!node.leaf() || node.level!=l) continue;
            const auto& p=node.coordinate;
            std::size_t id=p[0]+reference.nx*p[1];
            if constexpr(D==3) id+=reference.nx*reference.ny*p[2];
            require((node.cell==geo::invalidIndex)==(reference.scalar[id]==static_cast<double>(geo::CellType::Dry)),
                    "Wet leaf omitted or Dry leaf given a data slot");
        }
        if(l>=grid.levels().size()) continue;
        std::map<Key,std::array<double,links>> expected;
        for(std::size_t face=0; face<reference.nBoundaries; ++face)
            for(auto entry=reference.boundaryOffsets[face]; entry<reference.boundaryOffsets[face+1]; ++entry)
            {
                const auto id=reference.boundaryCells.at(entry);
                std::array<double,links> q{}; std::copy_n(reference.q.begin()+entry*links,links,q.begin());
                expected.emplace(Key{face,id},q);
            }
        for(auto cell=grid.levels()[l].all.begin; cell<grid.levels()[l].all.end; ++cell)
        {
            const auto p=grid.coordinate(cell);
            std::size_t id=p[0]+reference.nx*p[1];
            if constexpr(D==3) id+=reference.nx*reference.ny*p[2];
            require(static_cast<double>(cells.cellType[cell])==reference.scalar[id],
                    "Sparse tile scalar differs from uniform analysis");
            for(std::size_t face=0; face<boundary.nBoundaries; ++face)
            {
                const auto a=actual.find({face,cell}),b=expected.find({face,id});
                require((a==actual.end())==(b==expected.end()),"Sparse tile Face membership differs from uniform analysis");
                if(a!=actual.end())
                    for(std::size_t link=0; link<links; ++link)
                        require(a->second[link]==b->second[link],
                                "Sparse tile Q differs from uniform analysis");
            }
        }
    }
}

template<typename Analysis>
struct AnalysisCallback
{
    const std::function<void(Analysis&)>& analyze;
    void analysis(Analysis& result) const { analyze(result); }
};
} // namespace multilevel_check_detail

inline void checkTopology(const geo::MultilevelGrid<2>& grid)
{ multilevel_check_detail::topology(grid); }

inline void checkGeometry(const geo::MultilevelGrid<2>& grid, const std::function<void(geo::GeometryAnalysis2D<lat::LatticeType::D2Q9>&)>& analyze)
{ multilevel_check_detail::compareGeometry(grid,multilevel_check_detail::AnalysisCallback<geo::GeometryAnalysis2D<lat::LatticeType::D2Q9>>{analyze}); }

inline void checkTopology(const geo::MultilevelGrid<3,lat::LatticeType::D3Q15>& grid)
{ multilevel_check_detail::topology(grid); }

inline void checkGeometry(const geo::MultilevelGrid<3,lat::LatticeType::D3Q15>& grid, const std::function<void(geo::GeometryAnalysis3D<lat::LatticeType::D3Q15>&)>& analyze)
{ multilevel_check_detail::compareGeometry(grid,multilevel_check_detail::AnalysisCallback<geo::GeometryAnalysis3D<lat::LatticeType::D3Q15>>{analyze}); }

inline void checkTopology(const geo::MultilevelGrid<3>& grid)
{ multilevel_check_detail::topology(grid); }

inline void checkGeometry(const geo::MultilevelGrid<3>& grid, const std::function<void(geo::GeometryAnalysis3D<lat::LatticeType::D3Q19>&)>& analyze)
{ multilevel_check_detail::compareGeometry(grid,multilevel_check_detail::AnalysisCallback<geo::GeometryAnalysis3D<lat::LatticeType::D3Q19>>{analyze}); }

inline void checkTopology(const geo::MultilevelGrid<3,lat::LatticeType::D3Q27>& grid)
{ multilevel_check_detail::topology(grid); }

inline void checkGeometry(const geo::MultilevelGrid<3,lat::LatticeType::D3Q27>& grid, const std::function<void(geo::GeometryAnalysis3D<lat::LatticeType::D3Q27>&)>& analyze)
{ multilevel_check_detail::compareGeometry(grid,multilevel_check_detail::AnalysisCallback<geo::GeometryAnalysis3D<lat::LatticeType::D3Q27>>{analyze}); }

inline geo::MultilevelGrid<2> makeRegionGrid(int threads)
{
    const geo::BoundingBox2D domain{{0,0},{8,8}};
    const geo::Rectangle a({0,0},{5,8}),b({3,0},{8,8}); const auto physical=a|b;
    const auto left=geo::volumeRegion(geo::Circle({1.5,4},0.6));
    const auto right=geo::volumeRegion(geo::Circle({6.5,4},0.6));
    const auto hole=geo::volumeRegion(geo::Circle({1.5,4},0.2));
    geo::MultilevelOptions options; options.maxLevel=3; options.analysisTileWidth=4;
    const std::vector<geo::RefinementRule<2>> rules{
        geo::RefinementRule<2>::volume(left|right,1),
        geo::RefinementRule<2>::volume(left-hole,3)};
    omp_set_num_threads(threads);
    geo::MultilevelGrid<2> grid(domain,1,options); grid.generate(physical,rules);
    return grid;
}

} // namespace ntic::lbm::geometry::test
