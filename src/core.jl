const BivariateLike = Union{<:Tuple,<:AbstractVector}


function rotate_point(dx, dy, theta)
    # generate rotation matrix
    # (theta is degrees CCW from x-axis)
    R = RotMatrix{2}(-deg2rad(theta))
    # rotate points
    return R * SA[dx, dy]
end
