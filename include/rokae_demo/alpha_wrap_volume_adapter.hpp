#pragma once

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_3.h>
#include <CGAL/Alpha_wrap_3/internal/Point_set_oracle.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/version.h>

#include <array>
#include <chrono>
#include <cstddef>
#include <map>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace rokae_demo
{

// This is the only project component allowed to depend on CGAL Alpha_wrap_3
// internal APIs. CGAL 5.6.3 has no public API for retaining the volumetric
// Delaunay triangulation after wrapping. Keep this adapter version-pinned and
// re-audit it before upgrading CGAL.
class AlphaWrapVolumeAdapter
{
public:
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Point = Kernel::Point_3;
  using SurfaceMesh = CGAL::Surface_mesh<Point>;

  struct Tetrahedron
  {
    std::size_t id = 0;
    std::array<std::size_t, 4> vertices{};
    // Neighbor IDs refer only to retained (non-outside) finite cells.
    // -1 means that this face is on the wrap boundary.
    std::array<int, 4> neighbors{{-1, -1, -1, -1}};
  };

  struct Statistics
  {
    std::size_t finite_cells = 0;
    std::size_t outside_cells = 0;
    std::size_t non_outside_cells = 0;
    double alpha_wrap_ms = 0.0;
    double tetra_extract_ms = 0.0;
  };

private:
  using Oracle = CGAL::Alpha_wraps_3::internal::Point_set_oracle<Kernel>;
  using Wrapper = CGAL::Alpha_wraps_3::internal::Alpha_wrap_3<Oracle>;

  std::unique_ptr<Oracle> oracle_;
  std::unique_ptr<Wrapper> wrapper_;
  SurfaceMesh surface_mesh_;
  std::vector<Point> vertices_;
  std::vector<Tetrahedron> tetrahedra_;
  Statistics statistics_;

public:
  void run(const std::vector<Point>& points, double alpha, double offset)
  {
    if(points.empty() || alpha <= 0.0 || offset <= 0.0)
      throw std::invalid_argument(
          "Alpha wrap requires points and positive alpha/offset");

    // Destroy the wrapper before its referenced oracle.
    wrapper_.reset();
    oracle_.reset();
    surface_mesh_.clear();
    vertices_.clear();
    tetrahedra_.clear();
    statistics_ = {};

    oracle_ = std::make_unique<Oracle>(Kernel{});
    oracle_->add_point_set(points, CGAL::parameters::default_values());
    wrapper_ = std::make_unique<Wrapper>(*oracle_);
    const auto wrap_start = std::chrono::steady_clock::now();
    (*wrapper_)(alpha, offset, surface_mesh_);
    statistics_.alpha_wrap_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - wrap_start)
            .count();

    const auto extract_start = std::chrono::steady_clock::now();
    extractInteriorCells();
    statistics_.tetra_extract_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - extract_start)
            .count();
  }

  const SurfaceMesh& surfaceMesh() const { return surface_mesh_; }
  SurfaceMesh& surfaceMesh() { return surface_mesh_; }
  const std::vector<Point>& vertices() const { return vertices_; }
  const std::vector<Tetrahedron>& interiorCells() const
  {
    return tetrahedra_;
  }
  const Statistics& statistics() const { return statistics_; }

  // Point-location in the original Alpha Wrap triangulation makes validation
  // independent of a potentially expensive linear scan over tetrahedra.
  bool containsInRetainedCells(const Point& point) const
  {
    if(!wrapper_)
      throw std::logic_error("AlphaWrapVolumeAdapter has not been run");
    const auto& triangulation = wrapper_->triangulation();
    const auto cell = triangulation.locate(point);
    return !triangulation.is_infinite(cell) &&
           !cell->info().is_outside;
  }

private:
  void extractInteriorCells()
  {
    const auto& triangulation = wrapper_->triangulation();
    using Triangulation =
        std::remove_cv_t<std::remove_reference_t<decltype(triangulation)>>;
    using CellHandle = typename Triangulation::Cell_handle;
    using VertexHandle = typename Triangulation::Vertex_handle;

    std::map<CellHandle, std::size_t> retained_ids;
    for(auto cell = triangulation.finite_cells_begin();
        cell != triangulation.finite_cells_end(); ++cell)
    {
      ++statistics_.finite_cells;
      // In CGAL 5.6.3, the final classification flag is stored in the cell
      // info. Deliberately retain every non-outside cell, including cells that
      // Alpha Wrap may have added while enforcing manifoldness.
      if(cell->info().is_outside)
      {
        ++statistics_.outside_cells;
        continue;
      }
      const std::size_t id = retained_ids.size();
      retained_ids.emplace(cell, id);
      ++statistics_.non_outside_cells;
    }

    std::map<VertexHandle, std::size_t> vertex_ids;
    tetrahedra_.resize(retained_ids.size());
    for(const auto& entry : retained_ids)
    {
      const CellHandle cell = entry.first;
      Tetrahedron& tetrahedron = tetrahedra_[entry.second];
      tetrahedron.id = entry.second;
      for(int i = 0; i < 4; ++i)
      {
        const VertexHandle vertex = cell->vertex(i);
        auto inserted = vertex_ids.emplace(vertex, vertex_ids.size());
        if(inserted.second)
          vertices_.push_back(vertex->point());
        tetrahedron.vertices[static_cast<std::size_t>(i)] =
            inserted.first->second;

        const auto neighbor = cell->neighbor(i);
        const auto neighbor_id = retained_ids.find(neighbor);
        if(neighbor_id != retained_ids.end())
          tetrahedron.neighbors[static_cast<std::size_t>(i)] =
              static_cast<int>(neighbor_id->second);
      }
    }
  }
};

}  // namespace rokae_demo
