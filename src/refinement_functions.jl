function circular_refinement_function(center, inner_radius, buffer_radius, inner_density, outer_density, lx, ly)
    f = let center = center, inner_radius = inner_radius, buffer_radius = buffer_radius,
            inner_density = inner_density, outer_density = outer_density, lx = lx, ly = ly

        function(x)
            xp = periodic_to_base_point(x, lx, ly)
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

function y_refinement_function(center_y, length, buffer_length, inner_density, outer_density, lx, ly)
    f = let center_y = center_y, length = length, buffer_length = buffer_length,
            inner_density = inner_density, outer_density = outer_density

        function(x)
            xp = periodic_to_base_point(x, lx, ly)
            d = norm(xp.y - center_y)
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

function x_refinement_function(center_x, length, buffer_length, inner_density, outer_density, lx, ly)
    f = let center_x = center_x, length = length, buffer_length = buffer_length,
            inner_density = inner_density, outer_density = outer_density

        function(x)
            xp = periodic_to_base_point(x, lx, ly)
            d = norm(xp.x - center_x)
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
