#include "mpi/partition/partition.hpp"
#include "vtk_output_3d.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace geo=ntic::lbm::geometry;
namespace part=ntic::lbm::mpi;

namespace
{

using PlanePoint=std::array<double,2>;

struct PortDimensions
{
    std::size_t faceID=0;
    double areaEquivalentDiameter=0;
    double minimumFeretWidth=0;
};

geo::Point cross(const geo::Point& a,const geo::Point& b)
{
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}

geo::Point unit(const geo::Point& vector)
{
    const double length=std::hypot(vector[0],vector[1],vector[2]);
    if(!(length>0) || !std::isfinite(length))
        throw std::runtime_error("Port has an invalid plane normal or projection basis.");
    return {vector[0]/length,vector[1]/length,vector[2]/length};
}

double minimumFeretWidth(std::vector<PlanePoint> points)
{
    for(const auto& point:points)
        if(!std::isfinite(point[0]) || !std::isfinite(point[1]))
            throw std::runtime_error("Port projection contains a nonfinite coordinate.");
    std::sort(points.begin(),points.end());
    points.erase(std::unique(points.begin(),points.end()),points.end());
    if(points.size()<3) throw std::runtime_error("Port has fewer than three distinct projected vertices.");

    PlanePoint low=points.front(),high=points.front();
    for(const auto& point:points)
        for(std::size_t d=0; d<2; ++d)
        { low[d]=std::min(low[d],point[d]); high[d]=std::max(high[d],point[d]); }
    const double scale=std::max(high[0]-low[0],high[1]-low[1]);
    const double coordinateTolerance=128*std::numeric_limits<double>::epsilon()*scale;
    const double turnTolerance=coordinateTolerance*scale;
    if(!(scale>0) || !std::isfinite(scale) || !std::isfinite(turnTolerance))
        throw std::runtime_error("Port projection has invalid extent.");
    const auto turn=[](const PlanePoint& a,const PlanePoint& b,const PlanePoint& c) {
        return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
    };
    const auto append=[&](std::vector<PlanePoint>& chain,const PlanePoint& point) {
        while(chain.size()>1 && turn(chain[chain.size()-2],chain.back(),point)<=turnTolerance)
            chain.pop_back();
        chain.push_back(point);
    };
    std::vector<PlanePoint> lower,upper;
    lower.reserve(points.size()); upper.reserve(points.size());
    for(const auto& point:points) append(lower,point);
    for(auto point=points.rbegin(); point!=points.rend(); ++point) append(upper,*point);
    lower.pop_back(); upper.pop_back();
    lower.insert(lower.end(),upper.begin(),upper.end());
    const auto& hull=lower;
    if(hull.size()<3) throw std::runtime_error("Port projection is collinear or degenerate.");

    // A convex polygon's minimum supporting-line distance occurs normal to
    // an edge. Use the actual projected contour, not a bounding-box diagonal.
    double width=std::numeric_limits<double>::infinity();
    for(std::size_t edge=0; edge<hull.size(); ++edge)
    {
        const auto& a=hull[edge]; const auto& b=hull[(edge+1)%hull.size()];
        const double dx=b[0]-a[0],dy=b[1]-a[1],length=std::hypot(dx,dy);
        if(!(length>coordinateTolerance) || !std::isfinite(length))
            throw std::runtime_error("Port convex hull has a degenerate edge.");
        const double nx=-dy/length,ny=dx/length;
        double minimum=std::numeric_limits<double>::infinity(),maximum=-minimum;
        for(const auto& point:hull)
        {
            const double projection=nx*point[0]+ny*point[1];
            minimum=std::min(minimum,projection); maximum=std::max(maximum,projection);
        }
        width=std::min(width,maximum-minimum);
    }
    if(!(width>coordinateTolerance) || !std::isfinite(width))
        throw std::runtime_error("Port minimum Feret width is not positive and finite.");
    return width;
}

std::vector<PortDimensions> measurePorts(const geo::STLGeometry& geometry)
{
    const auto& features=geometry.boundaryFeatures;
    const auto& topology=geometry.topology;
    if(features.nPlanarFaces!=7 || features.planarFaces.size()!=7)
        throw std::runtime_error("Expected seven planar ports in the branched-channel STL.");
    if(features.facetFaceIDs.size()!=topology.geometry.facetVertexIDs.size() ||
       features.facetFaceIDs.size()!=topology.facetGeometry.size())
        throw std::runtime_error("Incomplete port facet metadata.");
    std::vector<std::vector<std::size_t>> vertices(features.nPlanarFaces+1);
    std::vector<geo::Point> normals(features.nPlanarFaces+1);
    for(std::size_t f=0; f<features.facetFaceIDs.size(); ++f)
    {
        const auto face=features.facetFaceIDs[f];
        if(face==0) continue;
        if(face>=vertices.size()) throw std::runtime_error("Invalid planar port Face ID.");
        const auto& facet=topology.facetGeometry[f];
        for(std::size_t d=0; d<3; ++d) normals[face][d]+=facet.area*facet.normal[d];
        const auto& ids=topology.geometry.facetVertexIDs[f];
        vertices[face].insert(vertices[face].end(),ids.begin(),ids.end());
    }
    std::vector<PortDimensions> result;
    result.reserve(features.planarFaces.size());
    for(const auto& face:features.planarFaces)
    {
        auto& ids=vertices.at(face.faceID);
        std::sort(ids.begin(),ids.end()); ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
        if(ids.size()<3 || !(face.area>0) || !std::isfinite(face.area))
            throw std::runtime_error("Port has insufficient vertices or invalid area.");
        const auto normal=unit(normals[face.faceID]);
        std::size_t axis=0;
        for(std::size_t d=1; d<3; ++d)
            if(std::abs(normal[d])<std::abs(normal[axis])) axis=d;
        geo::Point reference{}; reference[axis]=1;
        const auto u=unit(cross(normal,reference)),v=unit(cross(normal,u));
        // Subtract a nearby origin before projection to reduce cancellation.
        const auto& origin=topology.geometry.vertices.at(ids.front());
        std::vector<PlanePoint> projected;
        projected.reserve(ids.size());
        for(auto id:ids)
        {
            const auto& point=topology.geometry.vertices.at(id);
            PlanePoint position{};
            for(std::size_t d=0; d<3; ++d)
            {
                const double delta=point[d]-origin[d];
                position[0]+=delta*u[d]; position[1]+=delta*v[d];
            }
            projected.push_back(position);
        }
        const double areaDiameter=2*std::sqrt(face.area/std::acos(-1.0));
        if(!(areaDiameter>0) || !std::isfinite(areaDiameter))
            throw std::runtime_error("Port area-equivalent diameter is invalid.");
        result.push_back({face.faceID,areaDiameter,minimumFeretWidth(std::move(projected))});
    }
    return result;
}

} // namespace

int main()
{
    try
    {
        std::cout << "STL branched channel: multilevel partition metadata\n"
                  << "Four topology levels (0,1,2,3), eight equal-capacity logical ranks.\n"
                  << "Ranks 0..3 on node 0, 4..7 on node 1; no actual MPI processes are launched.\n"
                  << "Resolution is fixed at 25 cells across the minimum port Feret width.\n"
                  << "Only Other Surface (Face 0) seeds wall bands. Port-rim cells may still refine.\n\n";

        const int threads=std::max(1,std::min(128,omp_get_max_threads()));
        omp_set_num_threads(threads);
        const auto file=std::filesystem::path(MPI_TEST_SOURCE_DIR)/
            "stls"/"smooth_irregular_branched_channel.stl";
        std::cout << "Reading and validating STL: " << file << '\n' << std::flush;
        auto data=ntic::lbm::stl::read(file.string());
        ntic::lbm::stl::FacetTopology topology;
        ntic::lbm::stl::validate(data,"full",topology);
        const geo::STLGeometry geometry(std::move(topology),geo::FlowType::Internal);
        const auto ports=measurePorts(geometry);
        const auto smallest=std::min_element(ports.begin(),ports.end(),[](const auto& a,const auto& b) {
            return a.minimumFeretWidth<b.minimumFeretWidth;
        });
        constexpr unsigned cellsAcrossMinimumPort=25;
        const double minimumDiameter=smallest->minimumFeretWidth;
        const double finestSpacing=minimumDiameter/cellsAcrossMinimumPort,coarseSpacing=8*finestSpacing;
        if(!(finestSpacing>0) || !std::isfinite(coarseSpacing))
            throw std::runtime_error("Derived port-based spacing is invalid.");
        std::cout << std::setprecision(12);
        for(const auto& port:ports)
            std::cout << "Port Face " << port.faceID << ": area-equivalent diameter=" << port.areaEquivalentDiameter
                      << ", minimum Feret width=" << port.minimumFeretWidth << '\n';
        std::cout << "Selected port Face " << smallest->faceID << ", Dmin=" << minimumDiameter
                  << ", fixed cells across Dmin=" << cellsAcrossMinimumPort << ", hmin=" << finestSpacing << '\n';
        for(unsigned level=0; level<=3; ++level)
            std::cout << "Level " << level << " spacing: " << std::ldexp(coarseSpacing,-static_cast<int>(level)) << '\n';
        double bandCells=0;
        std::cout << "Narrow wall-band thickness in finest cells (recommended 2): "; std::cin >> bandCells;
        const double narrowBand=bandCells*finestSpacing,middleBand=2*narrowBand,wideBand=4*narrowBand;
        if(!std::cin || !(bandCells>0) || !std::isfinite(bandCells) ||
           !(narrowBand>0) || !std::isfinite(wideBand))
            throw std::invalid_argument("Wall-band cell count must be finite and positive.");
        if(wideBand>=0.5*minimumDiameter)
            std::cout << "Warning: the widest band reaches half the minimum port width; thin branches may be fully refined.\n";
        const auto walls=geo::fluidRegion(geometry);
        const std::vector<geo::RefinementRule<3>> rules{
            geo::RefinementRule<3>::surface(walls,1,wideBand,{0}),
            geo::RefinementRule<3>::surface(walls,2,middleBand,{0}),
            geo::RefinementRule<3>::surface(walls,3,narrowBand,{0})};
        geo::MultilevelOptions gridOptions;
        gridOptions.maxLevel=3; gridOptions.analysisThreads=static_cast<unsigned>(threads);
        gridOptions.maxNodes=6000000; gridOptions.maxCells=4000000;
        std::cout << "Generating the internal channel grid.\n"
                  << "Limits: " << gridOptions.maxNodes << " tree nodes, " << gridOptions.maxCells << " retained slots.\n"
                  << "These allocation guards are not preallocations or a measured memory budget.\n"
                  << std::flush;
        geo::MultilevelGrid<3> grid(geometry.bounds,coarseSpacing,gridOptions);
        std::cout << "Initial level-0 topology cells: " << grid.nodes().size() << '\n' << std::flush;
        grid.generate(geometry,rules);

        part::PartitionOptions options;
        options.ranks=8; options.blockLevel=0; options.graphThreads=std::min(8,threads);
        options.hardware.computeRates.assign(8,1.0);
        options.hardware.nodeIds={0,0,0,0,1,1,1,1};
        // Retain the initial cost model defaults; these are not hardware measurements.
        std::cout << "Partitioning with METIS...\n" << std::flush;
        const auto plan=part::partition(grid,options);
        std::ostringstream description;
        description << std::setprecision(17) << "Internal smooth_irregular_branched_channel.stl.\n"
                    << "Minimum port Face=" << smallest->faceID << ", minimum Feret Dmin=" << minimumDiameter
                    << ", fixed cells across Dmin=" << cellsAcrossMinimumPort << ".\n"
                    << "hmin=" << finestSpacing << ", h0=" << coarseSpacing << ", topology levels 0/1/2/3.\n"
                    << "Level-1 wall band=" << wideBand << ", level-2 wall band=" << middleBand
                    << ", level-3 wall band=" << narrowBand
                    << " (" << bandCells << " finest cells).\n"
                    << "Limits: " << gridOptions.maxNodes << " tree nodes and " << gridOptions.maxCells
                    << " retained slots; these guards are not preallocations or a measured memory budget.\n"
                    << "Face 0 is the only wall-band seed; planar port caps are excluded as seeds.\n"
                    << "Cells at port rims may refine through side-wall proximity and 2:1 balancing.\n"
                    << "The entire port cap is not guaranteed to remain coarse.";
        part::test::writePartitionOutputs(MPI_TEST_OUTPUT_DIR,"test_partition_branched_channel_vtks",
            grid,plan,options,description.str(),gridOptions.analysisThreads);
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
