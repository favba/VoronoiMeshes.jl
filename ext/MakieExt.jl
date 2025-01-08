module MakieExt

using VoronoiMeshes
import VoronoiMeshes: meshplot, meshplot!, diagramplot, diagramplot!
import Makie as plt

plt.@recipe(MeshPlot, mesh) do scene
    plt.Theme(
        color = :deepskyblue3,
        edgelinestyle = :solid,
        ghostedgelinestyle = :dash,
        linewidth = 1,
    )
end

const PeriodicMesh = MeshPlot{<:Tuple{<:AbstractVoronoiMesh{false}}}

function plt.plot!(plot::PeriodicMesh)
    mesh = plot.mesh
    both_edges = plt.lift(VoronoiMeshes.create_cell_linesegments, mesh)
    x_edges = plt.lift(both_edges) do both_edges
        both_edges[1][1]
    end
    y_edges = plt.lift(both_edges) do both_edges
        both_edges[1][2]
    end
    x_ghost_edges = plt.lift(both_edges) do both_edges
        both_edges[2][1]
    end
    y_ghost_edges = plt.lift(both_edges) do both_edges
        both_edges[2][2]
    end
    plt.linesegments!(plot, x_edges, y_edges, color = plot.color, linestyle = plot.edgelinestyle, linewidth = plot.linewidth)
    plt.linesegments!(plot, x_ghost_edges, y_ghost_edges, color = plot.color, linestyle = plot.ghostedgelinestyle, linewidth = plot.linewidth)
    plot
end

plt.@recipe(DiagramPlot, diagram) do scene
    plt.Theme(
        color = :deepskyblue3,
        linestyle = :solid,
        linewidth = 1,
    )
end

const PeriodicDiagram = DiagramPlot{<:Tuple{<:VoronoiDiagram{false}}}

function plt.plot!(plot::PeriodicDiagram)
    diagram = plot.diagram
    edges = plt.lift(VoronoiMeshes.create_diagram_linesegments, diagram)
    x_edges = plt.lift(edges) do edges
        edges[1]
    end
    y_edges = plt.lift(edges) do edges
        edges[2]
    end
    plt.linesegments!(plot, x_edges, y_edges, color = plot.color, linestyle = plot.linestyle, linewidth = plot.linewidth)
    plot
end

end
