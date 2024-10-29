compute_cell_area!(output,cells::Cells{false}) = compute_polygon_area_periodic!(output, cells.info.diagram.vertices, cells.vertices, cells.x_period, cells.y_period)

compute_cell_area(cells::Cells{false}) = compute_polygon_area_periodic(cells.info.diagram.vertices, cells.vertices,cells.x_period, cells.y_period)

compute_cell_area!(output,cells::Cells{true}) = compute_polygon_area_spherical!(output, cells.info.diagram.vertices, cells.vertices, cells.sphere_radius)

compute_cell_area(cells::Cells{true}) = compute_polygon_area_spherical(cells.info.diagram.vertices, cells.vertices, cells.sphere_radius)

compute_cell_centroid!(output,cells::Cells{false}) = compute_polygon_centroid_periodic!(output, cells.info.diagram.vertices, cells.vertices, cells.x_period, cells.y_period)

compute_cell_centroid(cells::Cells{false}) = compute_polygon_centroid_periodic(cells.info.diagram.vertices, cells.vertices,cells.x_period, cells.y_period)

compute_cell_centroid!(output,cells::Cells{true}) = compute_polygon_centroid_spherical!(output, cells.info.diagram.vertices, cells.vertices, cells.sphere_radius)

compute_cell_centroid(cells::Cells{true}) = compute_polygon_centroid_spherical(cells.info.diagram.vertices, cells.vertices, cells.sphere_radius)

compute_cell_longitude!(output,cells::Cells{false}) = fill!(output, zero(float_type(cells)))
compute_cell_longitude(cells::Cells{false}) = compute_longitude_periodic(cells.position)

compute_cell_longitude!(output,cells::Cells{true}) = compute_longitude!(output, cells.position)
compute_cell_longitude(cells::Cells{true}) = compute_longitude_spherical(cells.position)

compute_cell_latitude!(output,cells::Cells{false}) = fill!(output, zero(float_type(cells)))
compute_cell_latitude(cells::Cells{false}) = compute_latitude_periodic(cells.position)

compute_cell_latitude!(output,cells::Cells{true}) = compute_latitude!(output, cells.position)
compute_cell_latitude(cells::Cells{true}) = compute_latitude_spherical(cells.position)

compute_cell_normal!(output,cells::Cells{false}) = fill!(output, 𝐤)
compute_cell_normal(cells::Cells{false}) = VecArray(z = ones(float_type(cells), cells.n))

compute_cell_normal!(output,cells::Cells{true}) = tmap!(output, normalize, cells.position)
compute_cell_normal(cells::Cells{true}) = compute_cell_normal!(similar(cells.position), cells)

compute_cell_zonalVector!(output,cells::Cells{false}) = fill!(output, 𝐢)
compute_cell_zonalVector(cells::Cells{false}) = VecArray(x = ones(float_type(cells), cells.n))

compute_cell_zonalVector!(output,cells::Cells{true}) = compute_zonalVector!(output, cells.position)
compute_cell_zonalVector(cells::Cells{true}) = compute_zonalVector_spherical(cells.position)

compute_cell_meridionalVector!(output,cells::Cells{false}) = fill!(output, 𝐣)
compute_cell_meridionalVector(cells::Cells{false}) = VecArray(y = ones(float_type(cells), cells.n))

compute_cell_meridionalVector!(output,cells::Cells{true}) = compute_meridionalVector!(output, cells.position)
compute_cell_meridionalVector(cells::Cells{true}) = compute_meridionalVector_spherical(cells.position)

