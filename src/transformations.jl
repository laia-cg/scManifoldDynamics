
# ===========================================
# Functions related to the transformation of the data using vector fields
# ===========================================

"""
    mutable struct FunctionParams{T, P}
A container to store a transformation function (`tfunction`) and its parameters (`params`).

# Fields
- `tfunction`: The transformation function.
- `params`: Parameters for the transformation function.
"""
mutable struct FunctionParams{T, P}
    tfunction::T
    params::P
end

# -------------------------------------------------
# Transformation functions
# -------------------------------------------------

"""
    fgeneral_transformation(x, y; kwargs...)

Applies a general transformation to the coordinates `(x, y)`.

# Keyword Arguments
- `α`: Rotation angle (default: 0.0).
- `cx, cy`: Center of rotation (default: 0.0).
- `xlinshift, ylinshift`: Linear shifts in `x` and `y` (default: 0.0).
- `xpush_a, xpush_b, ypush_a, ypush_b`: Parameters for coordinate pushing (by default they all are set to 0.0).
- `xmoveto, ymoveto`: Target coordinates to move towards (by default they all are set to 0.0).
- `xmoveto_steps, ymoveto_steps`: Number of steps to reach the target coordinates (by default they all are set to 0.0).

# Returns
- A 2-element vector with the transformed `[x, y]` coordinates.
"""
function fgeneral_transformation(x, y;
    α::Number = 0.0, 
    cx::Number = 0.0, 
    cy::Number = 0.0, 
    xlinshift::Number = 0.0, 
    ylinshift::Number = 0.0,
    xpush_a::Number = 0.0,
    xpush_b::Number = 0.0,
    ypush_a::Number = 0.0,
    ypush_b::Number = 0.0,
    xmoveto::Number = 0.0, 
    xmoveto_steps::Number = 0.0,
    ymoveto::Number = 0.0, 
    ymoveto_steps::Number = 0.0
    )
    x_coord = cos(α) * (x - cx) - sin(α) * (y - cy) + cx
    y_coord = sin(α) * (x - cx) + cos(α) * (y - cy) + cy

    if xpush_a != 0
        x_coord += xpush_a / (x + xpush_b)
    end

    if ypush_a != 0
        y_coord += ypush_a / (y + ypush_b)
    end

    if xmoveto_steps != 0
        x_coord += (xmoveto - x) / xmoveto_steps
    end

    if ymoveto_steps != 0
        y_coord += (ymoveto - y) / ymoveto_steps
    end

    x_coord += xlinshift
    y_coord += ylinshift

    return [x_coord, y_coord]
end


"""
    params_information(f)

Generates a string summarizing the parameters of a given transformation function `f`.

# Arguments
- `f`: A `FunctionParams` object that contains the transformation function (`tfunction`) and its parameters (`params`).

# Returns
- A string describing the parameters of the transformation.
  - If the transformation is a `divide_cluster` function, it includes:
    - Parameters of the dividing line.
    - Parameters of the transformation functions above and below the line.
  - Otherwise, it lists the relevant parameters for the function.
"""
function params_information(f::FunctionParams)
    if typeof(f.tfunction) == typeof(divide_cluster)
        main_params_str = join(["$key=$(repr(value))" for (key, value) in pairs(f.params) if key ∉ [:fa, :fb]], ", ")
        subtitle1 = join(["$key=$(repr(value))" for (key, value) in pairs(f.params.fa.params)], ", ")
        subtitle2 = join(["$key=$(repr(value))" for (key, value) in pairs(f.params.fb.params)], ", ")
        return "Parameters sep. line: $main_params_str \n Function above: $subtitle1 \n Function below: $subtitle2"
    else
        excluded_keys = [:fa, :fb, :movetobool, :pushbool]
        params_str = join(["$key=$(repr(value))" for (key, value) in pairs(f.params) if key ∉ excluded_keys], ", ")
        return "$params_str"
    end
end


"""
    fdivide_countinuous(x, y; a, b, A, k, y0, B, C)

Applies a continuous division transformation to the coordinates `(x, y)`.

# Arguments
- `x, y`: Input coordinates.
- `a, b`: Parameters for the `x` transformation: `x' = x + a / (x + b)`.
- `A, k, y0, B, C`: Parameters for the `y` transformation:
  - `y' = y + (A / (1 + exp(-k * (y + y0))) - B) / C`.

# Returns
- A 2-element vector `[x', y']` of transformed coordinates.
"""
function fdivide_countinuous(x,y;
    a::Number=0.0, b::Number=0.0,
    A::Number=0.0, k::Number=1.0, y0::Number=-2.0, 
    B::Number=0.0, C::Number=1.0)
    x_coord = x + a/(x + b)
    y_coord = y + (A/(1+exp(-k*(y+y0)))-B)/C  
    return [x_coord, y_coord] 
end





"""
    divide_cluster(x, y; m, n, fa, fb, sigma=0.5)

Divides the points `(x, y)` of a group of points (e.g a cell cluster) using a line and applies one of two transformations (`fa` or `fb`) based on the point's position relative to the line.

# Arguments
- `x, y`: Input coordinates.
- `m, n`: Parameters for the dividing line: `y = m * x + n`.
- `fa`: Transformation function to apply above the line.
- `fb`: Transformation function to apply below the line.
- `sigma`: Standard deviation of random noise added to `y` for probabilistic splitting (default: `0.5`).

# Returns
- Transformed coordinates `[x', y']` using `fa` if the point is above the line, or `fb` if below.
"""
function divide_cluster(x, y; 
    m::Number, 
    n::Number, 
    fa::FunctionParams, 
    fb::FunctionParams, 
    sigma::Number=0.5)

    d = rand(Normal(0.0, sigma))
    if (y + d) >= m * x + n
        return fa.tfunction(x, y; fa.params...)
    else
        return fb.tfunction(x, y; fb.params...)
    end
end



# -------------------------------------------------
# Data transformation utilities
# -------------------------------------------------

"""
    Vector_Field(lowdata, f)

Applies the vector field transformation defined by `f` on `lowdata`.

# Arguments
- `lowdata`: A matrix where each row is a 2D point `[x, y]`.
- `f`: A `FunctionParams` object defining the transformation function and parameters.

# Returns
- A matrix of transformed data, where each row corresponds to a transformed `[x, y]` point.
"""
function Vector_Field(lowdata::AbstractMatrix{T}, f::FunctionParams) where T <: Real
    if f.params === nothing
        transdata = hcat(f.tfunction.(lowdata[:,1], lowdata[:,2])...)'
    else
        transdata = hcat(f.tfunction.(lowdata[:,1], lowdata[:,2]; f.params...)...)'
    end
    return transdata
end




"""
    transform_data!(ntimepoints, ind_clusterchange, lowdata, c, f)

Applies a series of vector field transformations to a dataset across multiple time points, updating the data in place.

# Arguments
- `ntimepoints`: Number of time points for the transformation.
- `ind_clusterchange`: Indices of clusters to be transformed at each time step.
- `lowdata`: Initial 2D data matrix where each row is `[x, y]`.
- `c`: Cluster assignments for each data point.
- `f`: A list of `FunctionParams` objects, one per cluster, specifying the transformations.

# Returns
- An updated `latent_transformed_data` object in place with the transformed data for all time points.
"""
function transform_data!(ntimepoints::Int, ind_clusterchange::Vector{Int}, lowdata::AbstractMatrix{T}, ::Vector{Vector{Int}}, f::Vector{FunctionParams}) where T <: Real

    latent_transformed_data = fill(0.0f0, size(lowdata, 1), 2, ntimepoints) 
    latent_transformed_data[:,:,1] = lowdata

    trans_seed = 1234
    vec_clusterchange = [c[ind_clusterchange[1]]] 

    for tp in 2:ntimepoints
        for ct in ind_clusterchange
            latent_transformed_data[c[ct],:,tp] = Vector_Field(latent_transformed_data[c[ct],:,tp-1], f[ct])  
            vec_clusterchange = union(vec_clusterchange, c[ct])
        end
        σ = 0.1
        g = Normal(0, σ)
        re_lowdata = lowdata[setdiff(1:end, vec_clusterchange),:] 
        Random.seed!(trans_seed)
        latent_transformed_data[setdiff(1:end, vec_clusterchange),:,tp] = re_lowdata .+ rand(g,size(re_lowdata,1), size(re_lowdata,2))  
    end
end


"""
    transform_data(ntimepoints, ind_clusterchange, lowdata, c, f)

Generates a series of vector field transformations to a dataset across multiple time points, returning the results.

# Arguments
- `ntimepoints`: Number of time points for the transformation.
- `ind_clusterchange`: Indices of clusters to be transformed at each time step.
- `lowdata`: Initial 2D data matrix where each row is `[x, y]`.
- `c`: Cluster assignments for each data point.
- `f`: A list of `FunctionParams` objects, one per cluster, specifying the transformations.

# Returns
- `latent_transformed_data`: A 3D array where `[:,:,t]` gives the transformed data at time point `t`.
"""
function transform_data(ntimepoints::Int, ind_clusterchange::Vector{Int}, lowdata::AbstractMatrix{T}, c::Vector{Vector{Int}}, f::Vector{FunctionParams}) where T <: Real
    latent_transformed_data = fill(0.0f0, size(lowdata, 1), 2, ntimepoints) 
    latent_transformed_data[:,:,1] = lowdata

    trans_seed = 1234
    vec_clusterchange = [c[ind_clusterchange[1]]] 

    for tp in 2:ntimepoints
        for ct in ind_clusterchange
            latent_transformed_data[c[ct],:,tp] = Vector_Field(latent_transformed_data[c[ct],:,tp-1], f[ct])  
            vec_clusterchange = union(vec_clusterchange, c[ct])
        end
        σ = 0.1
        g = Normal(0, σ)
        re_lowdata = lowdata[setdiff(1:end, vec_clusterchange),:] 
        Random.seed!(trans_seed)
        latent_transformed_data[setdiff(1:end, vec_clusterchange),:,tp] = re_lowdata .+ rand(g,size(re_lowdata,1), size(re_lowdata,2))  
    end
    return latent_transformed_data
end





# function _transform_data!(ntimepoints::Int, ind_clusterchange::Vector{Int}, latent_transformed_data, lowdata, c, f)
#     trans_seed = 1234
#     vec_clusterchange = [c[ind_clusterchange[1]]] 

#     for tp in 2:ntimepoints
#         for ct in ind_clusterchange
#             latent_transformed_data[c[ct],:,tp] = Vector_Field(latent_transformed_data[c[ct],:,tp-1], f[ct])  
#             vec_clusterchange = union(vec_clusterchange, c[ct])
#         end
#         σ = 0.1
#         g = Normal(0, σ)
#         re_lowdata = lowdata[setdiff(1:end, vec_clusterchange),:] 
#         Random.seed!(trans_seed)
#         latent_transformed_data[setdiff(1:end, vec_clusterchange),:,tp] = re_lowdata .+ rand(g,size(re_lowdata,1), size(re_lowdata,2))  
#     end
# end

# function transform_data!(ntimepoints::Int, ind_clusterchange::Vector{Int}, lowdata, c, f)
#     latent_transformed_data = fill(0.0f0, size(lowdata, 1), 2, ntimepoints) 
#     latent_transformed_data[:,:,1] = lowdata
#     _transform_data!(ntimepoints, ind_clusterchange, latent_transformed_data, lowdata, c, f)
#     return nothing
# end

# function transform_data(ntimepoints::Int, ind_clusterchange::Vector{Int}, lowdata, c, f)
#     latent_transformed_data = fill(0.0f0, size(lowdata, 1), 2, ntimepoints) 
#     latent_transformed_data[:,:,1] = lowdata
#     _transform_data!(ntimepoints, ind_clusterchange, latent_transformed_data, lowdata, c, f)
#     return latent_transformed_data
# end
