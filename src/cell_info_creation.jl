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

