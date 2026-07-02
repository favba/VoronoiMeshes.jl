using VoronoiMeshes, DelaunayTriangulation, NCDatasets, TensorsLite


function main(xpt, ypt, dc)

    #Create a homogeneous regular hexagon mesh with x_period = xpt, y_period ≈ ypt, and cell distance "dc"
    perfect_mesh = create_planar_hex_mesh(xpt, ypt, dc)

    pm1 = perfect_mesh
    dc = dc

    #xpt and ypt are just "target" periods, here we retrieve the true x_period and y_period of the mesh
    xp::Float64 = pm1.x_period
    yp::Float64 = pm1.y_period

    #This is the density function we'll be using for the variable resolution mesh.
    f = function(p)
        a = 2pi/xp
        b = 2pi/yp
        return 3.5 + 2.5*cos(a*p.x)*sin(b*p.y)
    end

    #Relative tolerance for Lloyd's iteration
    rtol = 1e-7
    #Maximum time to spend doing Lloyd's iteration in minutes
    mt = 120.0

    #Create unitialized array to store the new distorted meshes
    more_dist_meshes = Vector{VoronoiMesh{false, N, Int32, Float64, TensorsLite.Zeros.Zero} where N}(undef, 6)

    #Our first distorted mesh adds some generator points on top of those from the homogeneous mesh, so as to create pentagons and heptagons, and also
    #perform LLoyd's iteration with the density function "f" to create distortion due to variable resolution.
    more_dist_meshes[1] = VoronoiMesh(
        vcat(
            pm1.cells.position,
            [(xp/2)𝐢 + (2.5yp/4)𝐣, (1xp/4)𝐢 + (1yp/4)𝐣, (3xp/4)𝐢 + (1yp/4)𝐣]
        ),
        xp, yp, density = f, rtol = rtol, max_time = mt
    )

    #Name of the file to store the first mesh.
    base_name_more_dist = string("mesh_xp",xpt,"_yp",ypt,"_dc",dc,"_more-dist_L1.nc")

    save(base_name_more_dist, more_dist_meshes[1])

    #This loop will create new meshes with roughly double the resolution of the previous one.
    for i in 2:6
        #Use both the cell position and edge positions of the previous mesh as the initial generator points of the new mesh.
        #At each iteration the time taken for the Lloyd's iteration to converge will increae, hence the max_time being a function of "i".
        more_dist_meshes[i] = VoronoiMesh(vcat(more_dist_meshes[i-1].cells.position, more_dist_meshes[i-1].edges.position), xp, yp, max_time = 4*i*mt, density=f, rtol=rtol)
        #Save the new mesh into a file with a similar name to the first one, replacing only "L1" to "Li" where i is the current loop iteration.
        save(replace(base_name_more_dist, "L1" => "L$i"), more_dist_meshes[i])
    end


    return nothing
end

main(parse(Float64, ARGS[1]), parse(Float64, ARGS[2]), parse(Float64, ARGS[3]))

