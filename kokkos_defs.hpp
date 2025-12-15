#pragma once
// kokkos_defs.hpp
//
// Unified Kokkos setup for Ludwig
// - Selects execution & memory spaces from compile-time backend macros
// - 1D View aliases (extend locally if you need 2D/3D)
// - C++20 templated-lambda `static_for` (device-friendly)
//
// Backends (define exactly one):
//   -DLUDWIG_BACKEND_CUDA
//   -DLUDWIG_BACKEND_HIP
//   -DLUDWIG_BACKEND_SYCL
//   -DLUDWIG_BACKEND_OPENMP
//   -DLUDWIG_BACKEND_SERIAL
//
// Kokkos must be built with the same backend.

#include <Kokkos_Core.hpp>

//------------------------------------------------------------------------------
// Backend selection -> default execution & memory spaces
//------------------------------------------------------------------------------
#if defined(LUDWIG_BACKEND_CUDA)
using ExecSpace = Kokkos::Cuda;
using MemSpace = Kokkos::CudaSpace;
#elif defined(LUDWIG_BACKEND_HIP)
using ExecSpace = Kokkos::HIP;
using MemSpace = Kokkos::HIPSpace;
#elif defined(LUDWIG_BACKEND_SYCL)
using ExecSpace = Kokkos::SYCL;
using MemSpace = Kokkos::SYCLDeviceUSMSpace;
#elif defined(LUDWIG_BACKEND_OPENMP)
using ExecSpace = Kokkos::OpenMP;
using MemSpace = Kokkos::HostSpace;
#elif defined(LUDWIG_BACKEND_SERIAL)
using ExecSpace = Kokkos::Serial;
using MemSpace = Kokkos::HostSpace;
#else
using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = typename ExecSpace::memory_space;
#endif

using Device = Kokkos::Device<ExecSpace, MemSpace>;
using RangePolicy = Kokkos::RangePolicy<ExecSpace>;
using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
using Member = typename TeamPolicy::member_type;

//------------------------------------------------------------------------------
// 1D View aliases
//------------------------------------------------------------------------------
template <class T> using View = Kokkos::View<T*, MemSpace>;
template <class T> using CView = Kokkos::View<const T*, MemSpace>;
template <class T> using HostView = Kokkos::View<T*, Kokkos::HostSpace>;

