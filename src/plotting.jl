# ===========================================
# Functions used for plotting 
# ===========================================

# -------------------------------------------
# Functions related to the process of data transformation
# -------------------------------------------

"""
    plot_vector_field(f; grid_xmin=-5, grid_xmax=5, grid_ymin=-5, grid_ymax=5, grid_step=1.0, 
                      show_title=true, lengthscale=1)

Plots a vector field based on the transformation function `f`.

# Arguments
- `f::FunctionParams`: Transformation function and its parameters.
- `grid_xmin::Number`, `grid_xmax::Number`: Grid limits for the x-axis (default: `-5`, `5`).
- `grid_ymin::Number`, `grid_ymax::Number`: Grid limits for the y-axis (default: `-5`, `5`).
- `grid_step::Number`: Step size for the grid (default: `1.0`).
- `show_title::Bool`: Whether to display the title on the plot (default: `true`).
- `lengthscale::Number`: Scaling factor for arrow lengths (default: `1`).

# Returns
- A plot of the vector field.
"""
function plot_vector_field(
    f::FunctionParams,
    grid_xmin::Number=-5, 
    grid_xmax::Number=5, 
    grid_ymin::Number=-5, 
    grid_ymax::Number=5, 
    grid_step::Number=1.0; 
    show_title::Bool=true, 
    lengthscale::Number=1
    )
    # personalised_lims::Bool=false,
    # xmin::Number=-5,
    # xmax::Number=5,
    # ymin::Number=-5,
    # ymax::Number=5
    x, y = meshgrid(grid_xmin:grid_step:grid_xmax, grid_ymin:grid_step:grid_ymax)

    us = hcat(f.tfunction.(x, y; f.params...)...)'[:,1]
    vs = hcat(f.tfunction.(x, y; f.params...)...)'[:,2]

    us = us .- x
    vs = vs .- y

    fig = Figure(resolution = (400, 400))
    ax = Axis(fig[1, 1], backgroundcolor = "white")

    arrows!(ax, x, y, us, vs, arrowsize = 10, lengthscale = lengthscale)
    if show_title
        if typeof(f.tfunction) == typeof(divide_cluster)
            main_params_str = join(["$key=$(repr(value))" for (key, value) in pairs(f.params) if key ∉ [:fa, :fb]], ", ")
            subtitle1 = join(["$key=$(repr(value))" for (key, value) in pairs(f.params.fa.params)], ", ")
            subtitle2 = join(["$key=$(repr(value))" for (key, value) in pairs(f.params.fb.params)], ", ")
            ax.title = "Parameters sep. line: $main_params_str \n Function above: $subtitle1 \n Function below: $subtitle2"
        else
            excluded_keys = [:fa, :fb, :movetobool, :pushbool] 
            params_str = join(["$key=$(repr(value))" for (key, value) in pairs(f.params) if key ∉ excluded_keys], ", ")  # Convert the parameters into a string so we can plot it in the title
            ax.title = "$params_str"
        end
    end
    # if personalised_lims
    #     xlims!(ax, xmin, xmax)  
    #     ylims!(ax, ymin, ymax) 
    # end
    return fig
end

"""
    meshgrid(x, y)

Creates a 2D meshgrid from vectors `x` and `y`.

# Arguments
- `x::Vector`: Vector representing the x-axis.
- `y::Vector`: Vector representing the y-axis.

# Returns
- A tuple `(X, Y)` where `X` and `Y` are the 2D meshgrid arrays.
"""
function meshgrid(x::AbstractVector{<:Number}, y::AbstractVector{<:Number})
    return (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
end


"""
    plot_vector_field!(ax, f; grid_xmin=-5, grid_xmax=5, grid_ymin=-5, grid_ymax=5, grid_step=1.0, show_title=true)

Adds a vector field plot to an existing axis `ax` based on the transformation function `f`.

# Arguments
- `ax`: The existing axis object to which the vector field will be added.
- `f::FunctionParams`: Transformation function and its parameters.
- `grid_xmin::Int`, `grid_xmax::Int`: Grid limits for the x-axis (default: `-5`, `5`).
- `grid_ymin::Int`, `grid_ymax::Int`: Grid limits for the y-axis (default: `-5`, `5`).
- `grid_step::Float64`: Step size for the grid (default: `1.0`).
- `show_title::Bool`: Whether to display the title on the plot (default: `true`).

# Returns
- The previous base plot with the new vector field plotted in the assigned ax position.
"""
function plot_vector_field!(
    ax;
    f, 
    grid_xmin::Int=-5, 
    grid_xmax::Int=5, 
    grid_ymin::Int=-5,
    grid_ymax::Int=5, 
    grid_step::Float64=1.0, 
    show_title::Bool=true)

    x, y = meshgrid(grid_xmin:grid_step:grid_xmax, grid_ymin:grid_step:grid_ymax)

    us = hcat(f.tfunction.(x, y; f.params...)...)'[:,1]
    vs = hcat(f.tfunction.(x, y; f.params...)...)'[:,2]

    us = us .- x
    vs = vs .- y

    arrows!(ax, x, y, us, vs, arrowsize = 10, lengthscale = 1)
    if show_title
        if typeof(f.tfunction) == typeof(divide_cluster)
            main_params_str = join(["$key=$(repr(value))" for (key, value) in pairs(f.params) if key ∉ [:fa, :fb]], ", ")
            subtitle1 = join(["$key=$(repr(value))" for (key, value) in pairs(f.params.fa.params)], ", ")
            subtitle2 = join(["$key=$(repr(value))" for (key, value) in pairs(f.params.fb.params)], ", ")
            ax.title = "Parameters sep. line: $main_params_str \n Function above: $subtitle1 \n Function below: $subtitle2"
        else
            excluded_keys = [:fa, :fb, :movetobool, :pushbool] 
            params_str = join(["$key=$(repr(value))" for (key, value) in pairs(f.params) if key ∉ excluded_keys], ", ")
            ax.title = "$params_str"
        end
    end
    return
end


"""
    plot_tf_data(latent_transformed_data; plot_width=300, plot_height=200)

Plots scatter plots of transformed data at different time points.

# Arguments
- `latent_transformed_data::Array`: 3D array of transformed data, where `[:,:,t]` is the data at time point `t`.
- `plot_width::Number`, `plot_height::Number`: Dimensions of the plot (default: `300`, `200`).

# Returns
- A scatter plot showing transformed data across time points.
"""
function plot_tf_data(
    latent_transformed_data::AbstractArray{<:Number}; 
    plot_width::Number=300, 
    plot_height::Number=200
)
    latent_tf_data = cat(collect(latent_transformed_data[:,:, tp] for tp in 1:ntimepoints)..., dims=1)

    df_scatter_transformed = DataFrames.DataFrame(z1 = latent_tf_data[:,1], 
                                            z2 = latent_tf_data[:,2], 
                                            timepoint = cat(collect(fill(i, size(latent_transformed_data[:,1,1], 1)) for i in 1:ntimepoints)..., dims=1),
                                            cloud  = repeat(df_toy_data.label, outer=ntimepoints), 
    )

    plot_transformed_data = df_scatter_transformed |>  
    @vlplot(width=plot_width,
            height=plot_height,
            :circle, 
            x={:z1, title="z1", axis={titleFontSize=20, labelFontSize=20, tickCount=5}}, 
            y={:z2, title="z2", axis={titleFontSize=20, labelFontSize=20, tickCount=5}}, 
            column = {:"timepoint:n", axis={title="timepoint"}}, 
            color={"cloud:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
            #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
            size={value=35},
            labelFontSize=50,
            config={legend={titleFontSize=20, labelFontSize=20}}
    )
    return plot_transformed_data
end



"""
    plot_tf_data_t(latent_transformed_data, transformation_function)

Plots scatter plots of transformed data over time, with a title derived from the parameters of the transformation function `transformation_function`.

# Arguments
- `latent_transformed_data::Array`: 3D array of transformed data, where `[:,:,t]` is the data at time point `t`.
- `transformation_function::FunctionParams`: Transformation function and its parameters.

# Returns
- A scatter plot showing transformed data across time points, with a title describing the transformation parameters.
"""


function plot_tf_data_t(
    latent_transformed_data::AbstractArray{<:Number};
    transformation_function::FunctionParams, 
    cell_annotation::AbstractVector,
    ntimepoints::Int=4
    )

    latent_tf_data = cat(collect(latent_transformed_data[:,:, tp] for tp in 1:ntimepoints)..., dims=1)

    df_scatter_transformed = DataFrames.DataFrame(z1 = latent_tf_data[:,1], 
                                            z2 = latent_tf_data[:,2], 
                                            timepoint = cat(collect(fill(i, size(latent_transformed_data[:,1,1], 1)) for i in 1:ntimepoints)..., dims=1),
                                            cloud  = repeat(cell_annotation, outer=ntimepoints), 
    )

    
    title_str = ""
    if transformation_function !== nothing
        if typeof(transformation_function.tfunction) == typeof(divide_cluster)
            main_params_str = join(["$key=$(repr(value))" for (key, value) in pairs(transformation_function.params) if key ∉ [:fa, :fb]], ", ")
            subtitle1 = join(["$key=$(repr(value))" for (key, value) in pairs(transformation_function.params.fa.params)], ", ")
            subtitle2 = join(["$key=$(repr(value))" for (key, value) in pairs(transformation_function.params.fb.params)], ", ")
            title_str = "Parameters sep. line: $main_params_str \n Function above: $subtitle1 \n Function below: $subtitle2"
        else
            excluded_keys = [:fa, :fb, :movetobool, :pushbool]
            params_str = join(["$key=$(repr(value))" for (key, value) in pairs(transformation_function.params) if key ∉ excluded_keys], ", ")
            title_str = params_str
        end
    end

    plot_transformed_data = df_scatter_transformed |>  
    @vlplot(title=(title_str == "" ? nothing : (text=title_str, fontSize=20)), #title={"text"="$title", "fontSize"=20},
            width=220,
            height=220,
            :circle, 
            x={:z1, title="z1", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
            y={:z2, title="z2", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
            column = {:"timepoint:n", title=nothing}, 
            color={"cloud:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
            #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
            size={value=35},
            labelFontSize=50,
            config={legend={titleFontSize=20, labelFontSize=20}}
    )
    return plot_transformed_data
end



"""
    plot_loesses(loss_registry)

Plots loss trends over epochs for total loss, KL divergence, and reconstruction loss.

# Arguments
- `loss_registry::Dict`: A dictionary with keys `"total_loss"`, `"kl_z"`, and `"reconstruction"`, each containing a vector of loss values over epochs.

# Returns
- A line plot showing the three components of the VAE loss over epochs.
"""
function plot_loesses(loss_registry::Dict)
    epochs = length(loss_registry["total_loss"])
    df_loss = DataFrame(
        epoch = repeat(1:epochs, 3),
        loss_value = vcat(loss_registry["total_loss"], loss_registry["kl_z"], loss_registry["reconstruction"]),
        loss_type = repeat(["Total Loss", "KL Divergence", "Reconstruction"], inner=epochs)
    )

    plot_loss = df_loss |>
    @vlplot(
        :line,
        width=400,
        height=300,
        x={field=:epoch, scale = {zero = false}, type="quantitative", title="Epoch", axis={titleFontSize=12, labelFontSize=12}}, # 
        y={field=:loss_value, scale = {zero = false}, type="quantitative", title="Loss Value", axis={titleFontSize=12, labelFontSize=12}}, #
        column = {"loss_type:n", header={title = nothing, labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}},
        resolve={scale={x="independent",y="independent"}}
    )
    return plot_loss
end

