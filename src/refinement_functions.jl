function circular_refinement_function(lx, ly; center = (lx/2)𝐢 + (ly/2)𝐣, inner_radius = min(lx/2, ly/2)/2, buffer_radius = min(lx/2, ly/2), inner_density = 3.0, outer_density = 1.0)
    f = let center = center, inner_radius = inner_radius, buffer_radius = buffer_radius,
            inner_density = inner_density, outer_density = outer_density, lx = lx, ly = ly

        function(x)
            xp = closest(center, x, lx, ly)
            d = norm(xp - center)
            if d <= inner_radius
                return inner_density
            elseif d <= buffer_radius
                l = buffer_radius - inner_radius
                p = (d - inner_radius) / l
                return inner_density*(1 - p) + outer_density*p
            else
                return outer_density
            end
        end

    end
    return f
end

function y_refinement_function(lx, ly; center_y = ly/2, length = ly/4, buffer_length = ly/2, inner_density = 3.0, outer_density = 1.0)
    f = let center_y = center_y, length = length, buffer_length = buffer_length,
            inner_density = inner_density, outer_density = outer_density

        function(x)
            yp = closest(center_y, x.y, ly)
            d = norm(yp - center_y)
            if d <= length
                return inner_density
            elseif d <= buffer_length
                l = buffer_length - length
                p = (d - length) / l
                return inner_density*(1 - p) + outer_density*p
            else
                return outer_density
            end
        end
    end
    return f
end

function x_refinement_function(lx, ly; center_x = lx/2, length = lx/4, buffer_length = lx/2, inner_density = 3.0, outer_density = 1.0)
    f = let center_x = center_x, length = length, buffer_length = buffer_length,
            inner_density = inner_density, outer_density = outer_density

        function(x)
            xp = closest(center_x, x.x, lx)
            d = norm(xp - center_x)
            if d <= length
                return inner_density
            elseif d <= buffer_length
                l = buffer_length - length
                p = (d - length) / l
                return inner_density*(1 - p) + outer_density*p
            else
                return outer_density
            end
        end
    end
    return f
end
