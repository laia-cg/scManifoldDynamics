using Pkg
Pkg.activate(".")

using CSV, DataFrames, VegaLite, Plots, Distributions, Random
using UMAP, TSne
using scVI  # Pkg.add(url="https://github.com/maren-ha/scVI.jl", rev="v0.1.0")
using Flux
using JLD2

using OptimalTransport
using KernelDensity
using Distances
using LinearAlgebra



# -------- GET ALL BASE VARIABLES READY --------------
picked_dataset = "Lymph"
picked_model = "tsne_sae"
transnum = "t12" 
old_transformations = false

if old_transformations
    folderdata = "/Users/laiacanal/Documents/TimeSeries_project/new-results-time-transformation-project/$(picked_dataset)/$(picked_model)/Data/"
    folderresults = "/Users/laiacanal/Documents/TimeSeries_project/new-results-time-transformation-project/$(picked_dataset)/$(picked_model)/Results/"
else
    folderdata = "./Data/Synthetic_data/$(picked_dataset)/$(picked_model)/"
    folderresults = "./Data/Synthetic_data/$(picked_dataset)/$(picked_model)/Results/"
end


cell_annotation = CSV.read("Data/Datasets/$(picked_dataset)_cell_annotation.csv", DataFrames.DataFrame)
cell_annotation = cell_annotation[:,:x]


save_folder = "/Users/laiacanal/Documents/TimeSeries_project/Revisions/hyperparameter_justification/timeseries/$(picked_dataset)_$(transnum)_$(picked_model)/quantitative_study"

if !isdir(save_folder)
    mkpath(save_folder)
    println("Created folder: $save_folder")
else
    println("Folder already exists: $save_folder")
end

# -------- READ TRANSFORMED DATA -----------------

data_dict = load(folderdata*"matrices_$(transnum).jld2")["$(transnum)"]
println(data_dict["transformation_info"][1])
tf_data = data_dict["generated"]
tf_ncells, tf_ngenes = size(tf_data)
latent_tf_data = data_dict["latent"]

latent_tf_data = min_max_scale(latent_tf_data) 

ntimepoints = 4

# tf_adata = AnnData(countmatrix=tf_data, 
#             celltypes = repeat(cell_annotation, outer=ntimepoints),
#             obs = DataFrame(timepoint =  cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1))
#     )

translognormcounts = log.(normalization(tf_data) .+ 1)  #logNormCounts(tf_data, 2)


# ------ PLOT DATA ----------
df_latent_tf_data = DataFrames.DataFrame(z1 = latent_tf_data[:,1], 
                                                z2 = latent_tf_data[:,2], 
                                                timepoint_name = cat(collect(fill("t" * string(i), length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
                                                timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
                                                Celltype  = repeat(cell_annotation, outer=ntimepoints), 
        )


df_t1 = df_latent_tf_data[df_latent_tf_data.timepoint .== 1, :]  # all cells from timepoint 1
df_t2 = df_latent_tf_data[df_latent_tf_data.timepoint .== 2, :]
df_t3 = df_latent_tf_data[df_latent_tf_data.timepoint .== 3, :]
df_t4 = df_latent_tf_data[df_latent_tf_data.timepoint .== 4, :]
        
plot_transformed_data = df_latent_tf_data |>  
@vlplot(width=500,
        height=400,
        :circle, 
        x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        column = {:"timepoint_name:n",  header={title = " ", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
        # color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
        color={"Celltype:n", scale={range=tableau20_extended}, legend={disable=true, title="Clusters",orient="right"}},
        #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
        size={value=10},
        labelFontSize=50,
        config={legend={titleFontSize=20, labelFontSize=20}}
)

plot_file = joinpath(save_folder, "original_transformation.png")
VegaLite.save(plot_file, plot_transformed_data)


# ------ KERNEL ESTIMATION ------
timepoint1_data = [df_t1.z1 df_t1.z2]
timepoint2_data = [df_t2.z1 df_t2.z2]
timepoint3_data = [df_t3.z1 df_t3.z2]
timepoint4_data = [df_t4.z1 df_t4.z2]

kde_npoints = 64 #100 #64
kde_t1 = kde(timepoint1_data; npoints=(kde_npoints,kde_npoints))  # or kde((timepoint1_data[:,1], timepoint1_data[:,2]))
kde_t2 = kde(timepoint2_data; npoints=(kde_npoints,kde_npoints))  
kde_t3 = kde(timepoint3_data; npoints=(kde_npoints,kde_npoints))  
kde_t4 = kde(timepoint4_data; npoints=(kde_npoints,kde_npoints))  

plot_kde_t1 = plot(kde_t1,  linewidth=2, legend = false)
plot_kde_t2 = plot(kde_t2,  linewidth=2, legend = false)
plot_kde_t3 = plot(kde_t3, linewidth=2, legend = false)
plot_kde_t4 = plot(kde_t4, linewidth=2, legend = false)

all_kde_plots = plot(plot_kde_t1, plot_kde_t2, plot_kde_t3, plot_kde_t4, layout=(1,4), size=(2146,519))

Plots.savefig(all_kde_plots, save_folder*"/kde_estimation.png")

# ===============================================
# Compute Earth Mover’s (Wasserstein) Distance
# ===============================================

wd_t1_t2 = round(compute_wasserstein(kde_t1, kde_t2), digits = 4)
wd_t2_t3 = round(compute_wasserstein(kde_t2, kde_t3), digits = 4)
wd_t3_t4 = round(compute_wasserstein(kde_t3, kde_t4), digits = 4)
wd_t1_t4 = round(compute_wasserstein(kde_t1, kde_t4), digits = 4)

wd_original_df = DataFrame(
    Pair = ["t1-t2", "t2-t3", "t3-t4", "t1-t4"],
    WassersteinDistance = [wd_t1_t2, wd_t2_t3, wd_t3_t4, wd_t1_t4]
)

println("wd_original_df with kde_npoints = $(kde_npoints)")
println(wd_original_df)



# ------ option 1: KDE + Wasserstein ------- 



# ========================================
# Methods applied to the synthetic dataset
# ========================================


results = CSV.read(folderresults*"$(transnum)_df_reduced_dimensions_tf_data.csv", DataFrame)

study_method = "scVI" #tSNE
df_selected_results = results[results.Method .== study_method, :] 
df_selected_results.timepoint_name = cat(collect(fill("t" * string(i), length(cell_annotation)) for i in 1:ntimepoints)..., dims=1)

df_selected_results.z1 = min_max_scale(df_selected_results.z1)
df_selected_results.z2 = min_max_scale(df_selected_results.z2)

df_res_t1 = df_selected_results[df_selected_results.timepoint .== 1, :]  # all cells from timepoint 1
df_res_t2 = df_selected_results[df_selected_results.timepoint .== 2, :]
df_res_t3 = df_selected_results[df_selected_results.timepoint .== 3, :]
df_res_t4 = df_selected_results[df_selected_results.timepoint .== 4, :]
        
plot_transformed_data = df_selected_results |>  
@vlplot(width=500,
        height=400,
        :circle, 
        x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        column = {:"timepoint_name:n",  header={title = " ", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
        # color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
        color={"Celltype:n", scale={range=tableau20_extended}, legend={disable=true, title="Clusters",orient="right"}},
        #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
        size={value=10},
        labelFontSize=50,
        config={legend={titleFontSize=20, labelFontSize=20}}
)

plot_file = joinpath(save_folder, "$(study_method)_original_result.png")
VegaLite.save(plot_file, plot_transformed_data)


# ------ KERNEL ESTIMATION ------

timepoint1_res_data = [df_res_t1.z1 df_res_t1.z2]
timepoint2_res_data = [df_res_t2.z1 df_res_t2.z2]
timepoint3_res_data = [df_res_t3.z1 df_res_t3.z2]
timepoint4_res_data = [df_res_t4.z1 df_res_t4.z2]


kde_res_t1 = kde(timepoint1_res_data; npoints=(kde_npoints,kde_npoints))  # or kde((timepoint1_data[:,1], timepoint1_data[:,2]))
kde_res_t2 = kde(timepoint2_res_data; npoints=(kde_npoints,kde_npoints))  
kde_res_t3 = kde(timepoint3_res_data; npoints=(kde_npoints,kde_npoints))  
kde_res_t4 = kde(timepoint4_res_data; npoints=(kde_npoints,kde_npoints))  

plot_kde_res_t1 = plot(kde_res_t1,  linewidth=2, legend = false)
plot_kde_res_t2 = plot(kde_res_t2,  linewidth=2, legend = false)
plot_kde_res_t3 = plot(kde_res_t3, linewidth=2, legend = false)
plot_kde_res_t4 = plot(kde_res_t4, linewidth=2, legend = false)

all_kde_plots = plot(plot_kde_res_t1, plot_kde_res_t2, plot_kde_res_t3, plot_kde_res_t4, layout=(1,4), size=(2146,519))

Plots.savefig(all_kde_plots, save_folder*"/study_method_$(study_method)_kde_estimation.png")

# ===============================================
# Compute Earth Mover’s (Wasserstein) Distance
# ===============================================


wd_t1_t2 = round(compute_wasserstein(kde_res_t1, kde_res_t2), digits = 4)
wd_t2_t3 = round(compute_wasserstein(kde_res_t2, kde_res_t3), digits = 4)
wd_t3_t4 = round(compute_wasserstein(kde_res_t3, kde_res_t4), digits = 4)
wd_t1_t4 = round(compute_wasserstein(kde_res_t1, kde_res_t4), digits = 4)

results_df = DataFrame(
    Pair = ["t1-t2", "t2-t3", "t3-t4", "t1-t4"],
    WassersteinDistance = [wd_t1_t2, wd_t2_t3, wd_t3_t4, wd_t1_t4]
)

println("wd_results_df $(study_method) with kde_npoints = $(kde_npoints)")
println(results_df)



# ====================
# new util functions
# ====================


function min_max_scale(data)
    mins = minimum(data, dims=1)
    maxs = maximum(data, dims=1)
    # Avoid division by zero if maxs equals mins
    return (data .- mins) ./ (maxs .- mins .+ eps())
end


function create_grid(x, y)
    nx = length(x)
    ny = length(y)
    grid = zeros(nx * ny, 2)
    k = 1
    for j in 1:ny  # y coordinate
        for i in 1:nx  # x coordinate
            grid[k, :] = [x[i], y[j]]
            k += 1
        end
    end
    return grid
end


function compute_wasserstein(kde1, kde2, ε=0.1)
    xA, yA = kde1.x, kde1.y
    ZA_norm = kde1.density ./ sum(kde1.density)
    
    xB, yB = kde2.x, kde2.y
    ZB_norm = kde2.density ./ sum(kde2.density)
    
    gridA = create_grid(xA, yA)
    gridB = create_grid(xB, yB)
    
    pA = vec(ZA_norm)
    pB = vec(ZB_norm)
    
    C = pairwise(SqEuclidean(), gridA, gridB; dims=1)
    
    π = sinkhorn(pA, pB, C, ε)
    W = sum(π .* C)
    
    return W
end


# # ---- deprecated code -------

# # option 1: 
# kde_result = kde(latent_tf_data)
# plot(kde_result, title="Kernel Density Estimation", xlabel="Value", ylabel="Density", linewidth=2)


# # option 2: 
# kde_result = kde((latent_tf_data[:, 1], latent_tf_data[:, 2]))
# xrange = kde_result.x
# yrange = kde_result.y
# Z      = kde_result.density

# contourf(xrange, yrange, Z')

# contourf(x, y, z', title="2D Kernel Density Estimation", xlabel="X", ylabel="Y", color=:viridis)

# df_latent_tf_data = DataFrames.DataFrame(z1 = latent_tf_data[:,1], 
#     z2 = latent_tf_data[:,2], 
#     timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
#     Celltype  = repeat(cell_annotation, outer=ntimepoints), 
# )



xA = kde_t1.x
yA = kde_t1.y
ZA = kde_t1.density
# Usually sum(ZA) is not exactly 1; you can renormalize:
ZA_normed = ZA ./ sum(ZA)

xB = kde_t2.x
yB = kde_t2.y
ZB = kde_t2.density
ZB_normed = ZB ./ sum(ZB)


gridA = create_grid(xA, yA)
gridB = create_grid(xB, yB)


# Flatten
pA = vec(ZA_normed)
pB = vec(ZB_normed)

C_grid = pairwise(SqEuclidean(), gridA, gridB; dims=1)

ε = 0.1

# Option 1: Solve the entropically regularized optimal transport problem.
π_grid = sinkhorn(pA, pB, C_grid, ε)
W_grid = sum(π_grid .* C_grid)

# # Option 2: Calcuate the Wasserstein distance directly
# distance_grid = OptimalTransport.emd(pA, pB, C_grid, Tulip.Optimizer())
# W = sum(distance_grid .* C_grid)

# ------ option 2: raw data ------- (doesn't work so far)

# Suppose you want the distance between timepoint1_data and timepoint2_data
X = timepoint1_data  # shape (n1, 2)
Y = timepoint2_data  # shape (n2, 2)

μ = fill(1.0 / size(X, 1), size(X, 1))
ν = fill(1.0 / size(Y, 1), size(Y, 1))

X_scaled = min_max_scale(X)
Y_scaled = min_max_scale(Y)

C_scaled = pairwise(SqEuclidean(), X_scaled, Y_scaled; dims=1)

# Set the regularization parameter
ε =0.1

# Solve the entropically regularized optimal transport problem
π = sinkhorn(μ, ν, C_scaled, ε)

# Compute the regularized Wasserstein distance as the sum of π .* C
W = sum(π .* C_scaled)
