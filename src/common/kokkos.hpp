#pragma once

#include <Kokkos_Core.hpp>

#include <cstddef>


namespace ntic::lbm::kokkos {

//=============================================================================
// Backend configuration
//=============================================================================
//
// nticLBM uses one execution backend for the entire code base.
//
// Host and device memory are always treated as separate memory spaces.
// Unified/shared memory is intentionally not used.
//
// Current backend: CUDA
//
// -----------------------------------------------------------------------------
// CUDA
// -----------------------------------------------------------------------------
//
// using ExecSpace      = Kokkos::Cuda;
// using DeviceMemSpace = Kokkos::CudaSpace;
//
// using HostExecSpace  = Kokkos::Serial;
// using HostMemSpace   = Kokkos::CudaHostPinnedSpace;
//
// -----------------------------------------------------------------------------
// HIP
// -----------------------------------------------------------------------------
//
// using ExecSpace      = Kokkos::HIP;
// using DeviceMemSpace = Kokkos::HIPSpace;
//
// using HostExecSpace  = Kokkos::Serial;
// using HostMemSpace   = Kokkos::HIPHostPinnedSpace;
//
// -----------------------------------------------------------------------------
// SYCL
// -----------------------------------------------------------------------------
//
// using ExecSpace      = Kokkos::SYCL;
// using DeviceMemSpace = Kokkos::SYCLDeviceUSMSpace;
//
// using HostExecSpace  = Kokkos::Serial;
// using HostMemSpace   = Kokkos::SYCLHostUSMSpace;
//
// NOTE:
// SYCLDeviceUSMSpace and SYCLHostUSMSpace are separate device/host
// allocations. nticLBM does NOT use SYCLSharedUSMSpace.
//
// -----------------------------------------------------------------------------
// OpenMP
// -----------------------------------------------------------------------------
//
// using ExecSpace      = Kokkos::OpenMP;
// using DeviceMemSpace = Kokkos::HostSpace;
//
// using HostExecSpace  = Kokkos::Serial;
// using HostMemSpace   = Kokkos::HostSpace;
//
// -----------------------------------------------------------------------------
// OpenMP Target
// -----------------------------------------------------------------------------
//
// The exact OpenMP Target execution/memory-space type names depend on the
// Kokkos version in use. When enabling this backend, set:
//
//     ExecSpace
//     DeviceMemSpace
//
// to the corresponding OpenMP Target types provided by that Kokkos build,
// while keeping:
//
// using HostExecSpace = Kokkos::Serial;
// using HostMemSpace  = Kokkos::HostSpace;
//
//=============================================================================


// Current backend: CUDA

using ExecSpace      = Kokkos::Cuda;
using DeviceMemSpace = Kokkos::CudaSpace;

using HostExecSpace  = Kokkos::Serial;
using HostMemSpace   = Kokkos::CudaHostPinnedSpace;


//=============================================================================
// Device types
//=============================================================================

using DeviceType =
    Kokkos::Device<ExecSpace, DeviceMemSpace>;

using HostDeviceType =
    Kokkos::Device<HostExecSpace, HostMemSpace>;


//=============================================================================
// Index type
//=============================================================================

using Index = std::size_t;


//=============================================================================
// Execution policies
//=============================================================================

using RangePolicy =
    Kokkos::RangePolicy<
        ExecSpace,
        Kokkos::IndexType<Index>>;

using HostRangePolicy =
    Kokkos::RangePolicy<
        HostExecSpace,
        Kokkos::IndexType<Index>>;

using MDRangePolicy2D =
    Kokkos::MDRangePolicy<
        ExecSpace,
        Kokkos::Rank<2>,
        Kokkos::IndexType<Index>>;

using MDRangePolicy3D =
    Kokkos::MDRangePolicy<
        ExecSpace,
        Kokkos::Rank<3>,
        Kokkos::IndexType<Index>>;

using HostMDRangePolicy2D =
    Kokkos::MDRangePolicy<
        HostExecSpace,
        Kokkos::Rank<2>,
        Kokkos::IndexType<Index>>;

using HostMDRangePolicy3D =
    Kokkos::MDRangePolicy<
        HostExecSpace,
        Kokkos::Rank<3>,
        Kokkos::IndexType<Index>>;

using TeamPolicy =
    Kokkos::TeamPolicy<ExecSpace>;

using TeamMember =
    TeamPolicy::member_type;


//=============================================================================
// Device managed views
//=============================================================================

template <typename T>
using DeviceView1D =
    Kokkos::View<T*, DeviceType>;

template <typename T>
using DeviceView2D =
    Kokkos::View<T**, DeviceType>;

template <typename T>
using DeviceView3D =
    Kokkos::View<T***, DeviceType>;


//=============================================================================
// Device unmanaged views
//=============================================================================

template <typename T>
using DeviceUnmanagedView1D =
    Kokkos::View<
        T*,
        DeviceType,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template <typename T>
using DeviceUnmanagedView2D =
    Kokkos::View<
        T**,
        DeviceType,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template <typename T>
using DeviceUnmanagedView3D =
    Kokkos::View<
        T***,
        DeviceType,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;


//=============================================================================
// Host managed views
//=============================================================================

template <typename T>
using HostView1D =
    Kokkos::View<T*, HostDeviceType>;

template <typename T>
using HostView2D =
    Kokkos::View<T**, HostDeviceType>;

template <typename T>
using HostView3D =
    Kokkos::View<T***, HostDeviceType>;


//=============================================================================
// Host unmanaged views
//=============================================================================

template <typename T>
using HostUnmanagedView1D =
    Kokkos::View<
        T*,
        HostDeviceType,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template <typename T>
using HostUnmanagedView2D =
    Kokkos::View<
        T**,
        HostDeviceType,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template <typename T>
using HostUnmanagedView3D =
    Kokkos::View<
        T***,
        HostDeviceType,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;


//=============================================================================
// Team scratch memory
//=============================================================================

using TeamScratchSpace =
    TeamMember::scratch_memory_space;


template <typename T>
using TeamScratchView1D =
    Kokkos::View<
        T*,
        TeamScratchSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template <typename T>
using TeamScratchView2D =
    Kokkos::View<
        T**,
        TeamScratchSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template <typename T>
using TeamScratchView3D =
    Kokkos::View<
        T***,
        TeamScratchSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;


//=============================================================================
// Team scratch allocation helpers
//=============================================================================

template <typename T>
inline auto perTeamScratch1D(const std::size_t n0)
{
    return Kokkos::PerTeam(
        TeamScratchView1D<T>::shmem_size(n0));
}


template <typename T>
inline auto perTeamScratch2D(
    const std::size_t n0,
    const std::size_t n1)
{
    return Kokkos::PerTeam(
        TeamScratchView2D<T>::shmem_size(n0, n1));
}


template <typename T>
inline auto perTeamScratch3D(
    const std::size_t n0,
    const std::size_t n1,
    const std::size_t n2)
{
    return Kokkos::PerTeam(
        TeamScratchView3D<T>::shmem_size(n0, n1, n2));
}


} // namespace ntic::lbm::kokkos