
using Pkg
Pkg.activate(".")

using CSV, DataFrames, VegaLite, Plots, Distributions, Random
using UMAP, TSne
using scVI  # Pkg.add(url="https://github.com/maren-ha/scVI.jl", rev="v0.1.0")
using Flux
using JLD2

using OptimalTransport
using KernelDensity

using Random, Distributions, StatsBase, Plots
using LinearAlgebra


include("../src/scManifoldDynamics.jl") 
using .scManifoldDynamics  

# -------- GET ALL BASE VARIABLES READY --------------
picked_dataset = "PBMC"
picked_model = "tsne_sae"
transnum = "t1" 

folderdata = "./Data/Synthetic_data/$(picked_dataset)/$(picked_model)/"
folderresults = "./Data/Synthetic_data/$(picked_dataset)/$(picked_model)/Results/"

cell_annotation = CSV.read("Data/Datasets/$(picked_dataset)_cell_annotation.csv", DataFrames.DataFrame)
cell_annotation = cell_annotation[:,:x]


save_folder = "./Data/figures/sample_timepoints/$(picked_dataset)/"
if !isdir(save_folder)
    mkpath(save_folder)  
    println("Directory created: ", save_folder)
else
    println("Directory already exists: ", save_folder)
end

# -------- READ TRANSFORMED DATA -----------------

data_dict = load(folderdata*"matrices_$(transnum).jld2")["$(transnum)"]
println(data_dict["transformation_info"][1])
tf_data = data_dict["generated"]
tf_ncells, tf_ngenes = size(tf_data)
latent_tf_data = data_dict["latent"]


ntimepoints = 4

# tf_adata = AnnData(countmatrix=tf_data, 
#             celltypes = repeat(cell_annotation, outer=ntimepoints),
#             obs = DataFrame(timepoint =  cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1))
#     )

translognormcounts = log.(normalization(tf_data) .+ 1)  #logNormCounts(tf_data, 2)



if picked_dataset == "PBMC"
    big_TSNE_args = TSNEArgs(seed = 1013, perplexity = 200, pca_dims = 20)
    big_UMAP_args = UMAPArgs(seed = 87, pca_dims = 10) #105
elseif picked_dataset == "Zeisel"
    big_TSNE_args = TSNEArgs(seed = 1013, perplexity = 100, pca_dims = 20)
    big_UMAP_args = UMAPArgs(seed = 9711, pca_dims = 100)
elseif picked_dataset == "Lymph"
    big_TSNE_args = TSNEArgs(seed = 1013, perplexity = 50, pca_dims = 50)
    big_UMAP_args = UMAPArgs(seed = 9711, pca_dims = 40,  n_nbrs = 70)
end


# ------ PLOT DATA ----------


unique_clusters = unique(cell_annotation)
cluster_name_to_index = Dict(c => i for (i, c) in enumerate(unique_clusters))
cell_cluster_indices = [cluster_name_to_index[c] for c in cell_annotation]

# Now repeat for all timepoints
all_cell_cluster_indices = repeat(cell_cluster_indices, ntimepoints)

df_latent_tf_data = DataFrames.DataFrame(z1 = latent_tf_data[:,1], 
                                                z2 = latent_tf_data[:,2], 
                                                timepoint_name = cat(collect(fill("t" * string(i), length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
                                                timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
                                                Celltype  = repeat(cell_annotation, outer=ntimepoints), 
        )
df_latent_tf_data.cluster_index = all_cell_cluster_indices

plot_transformed_data = df_latent_tf_data |>  
@vlplot(width=500,
        height=400,
        :circle, 
        x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        column = {:"timepoint_name:n",  header={title = " ", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
        color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
        #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
        size={value=30},
        labelFontSize=50,
        config={legend={titleFontSize=20, labelFontSize=20}}
)
plot_transformed_data_file = joinpath(save_folder, "original_transformation_$(picked_model)_$(transnum).png")
VegaLite.save(plot_transformed_data_file, plot_transformed_data)



function sample_cells_by_index_v2(data::DataFrame,
                                   spatial::String;
                                   keep_ratio=[1.0, 0.75, 0.5, 0.25],
                                   cell_indices::Vector{Int},
                                   cluster_modes::Dict{Int, String})

    sampled_indices = Int[]
    subdata = data[cell_indices, :]

    timepoints = unique(subdata.timepoint)
    clusters = unique(subdata.cluster_index)

    for cl in clusters
        for (i, t) in enumerate(timepoints)
            cell_inds = findall(row -> row.cluster_index == cl && row.timepoint == t, eachrow(subdata))
            if isempty(cell_inds)
                continue
            end

            full_inds = cell_indices[cell_inds]

            # Should we modify this cluster?
            if cl ∉ keys(cluster_modes)
                append!(sampled_indices, full_inds)
                continue
            end

            mode = cluster_modes[cl]
            ratio = mode == "equal"    ? 0.6 : #rand(0.5:0.05:1.0) :
                    mode == "decrease" ? keep_ratio[i] :
                    mode == "increase" ? keep_ratio[end - i + 1] :
                    error("Unknown mode")

            n_keep = Int(round(length(full_inds) * ratio))

            if spatial == "uniform"
                selected = sample(full_inds, n_keep; replace=false)
            elseif spatial == "shaped"
                coords = reduce(hcat, [[data[i, :z1], data[i, :z2]] for i in full_inds])
                μ = mean(coords; dims=2)
                dists = sqrt.(sum((coords .- μ).^2, dims=1))

                # Exponential decay for higher prob near center
                weights = exp.(-vec(dists) ./ std(vec(dists)))
                probabilities = weights ./ sum(weights)
                selected = sample(full_inds, Weights(probabilities), n_keep; replace=false)
            else
                error("Unknown spatial mode")
            end

            append!(sampled_indices, selected)
        end
    end

    return sort(sampled_indices)
end


function select_cells_per_timepoint(timepoint::Vector{Int}, ncells_per_timepoint::Vector{Int})
    unique_timepoints = unique(timepoint)
    selected_indices = Int[]
    for (i, tp) in enumerate(unique_timepoints)
            tp_indices = findall(x -> x == tp, timepoint)
            ncells = ncells_per_timepoint[i]
            if length(tp_indices) >= ncells
                    selected = sample(tp_indices, ncells, replace=false)
            else
                    selected = tp_indices # If there are not enough cells, we take all of them
            end
            append!(selected_indices, selected)     
    end
    return selected_indices
end



# cluster_modes = Dict(3 => "increase", 1 => "increase")

# sampled_inds = sample_cells_by_index_v2(df_latent_tf_data,
#                                         "shaped";
#                                         cell_indices=cell_indices = collect(1:nrow(df_latent_tf_data)),
#                                         cluster_modes=cluster_modes)

# sampled_inds = sample_by_timepoint(df_latent_tf_data, 0.60)
ncells = length(cell_annotation)
timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1)
ncells_per_timepoint = [Int(0.7*ncells), Int(0.4*ncells), Int(0.6*ncells), Int(0.3*ncells)]
sampled_inds = select_cells_per_timepoint(timepoint, ncells_per_timepoint)


df_sampled = df_latent_tf_data[sampled_inds, :]

# ncells = length(cell_annotation)

plot_transformed_data = df_sampled |>  
@vlplot(width=500,
        height=400,
        :circle, 
        x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        column = {:"timepoint_name:n",  header={title = " ", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
        color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
        #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
        size={value=30},
        labelFontSize=50,
        config={legend={titleFontSize=20, labelFontSize=20}}
)
plot_transformed_data_file = joinpath(save_folder, "sampled_cells_$(picked_model)_$(transnum).png")
VegaLite.save(plot_transformed_data_file, plot_transformed_data)





tf_sample_data = tf_data[sampled_inds,:]

# elapsed_time_seconds = @elapsed begin
#     Random.seed!(big_TSNE_args.seed) 
#     tSNE = tsne(rescale(tf_sample_data, dims=1), 2, big_TSNE_args.pca_dims, 1000, big_TSNE_args.perplexity, progress=true); 
# end


translognormcounts_sample = log.(normalization(tf_sample_data) .+ 1)  #log
# pcs = scVI.prcomps(log.(normalization(tf_adata.countmatrix) .+ 1))
elapsed_time_seconds = @elapsed begin
        pcs = scVI.prcomps(translognormcounts_sample)
end
elapsed_time_seconds = @elapsed begin
    Random.seed!(big_UMAP_args.seed)
    umapout = umap(pcs[:,1:big_UMAP_args.pca_dims]', min_dist=0.5, n_neighbors = big_UMAP_args.n_nbrs) 
    # umapout_pcs_10 = umap(pcs[:,1:10]', min_dist=0.5) 
    # umapout_pcs_15 = umap(pcs[:,1:15]', min_dist=0.5) 
    # umapout = umap(Transpose(translognormcounts), min_dist=0.5)
    # umapout_origdata = umap(tf_data')
end


# df_sampled.tSNE1 = tSNE[:,1]
# df_sampled.tSNE2 = tSNE[:,2]


df_sampled.UMAP1 = umapout[1,:]
df_sampled.UMAP2 = umapout[2,:]

plot_transformed_data = df_sampled |>  
@vlplot(width=500,
        height=400,
        :circle, 
        x={:UMAP1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        y={:UMAP2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        column = {:"timepoint_name:n",  header={title = " ", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
        color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
        #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
        size={value=30},
        labelFontSize=50,
        config={legend={titleFontSize=20, labelFontSize=20}}
)
plot_transformed_data_file = joinpath(save_folder, "sampled_umap_$(picked_model)_$(transnum).png")
VegaLite.save(plot_transformed_data_file, plot_transformed_data)





translognormcounts = log.(normalization(tf_data) .+ 1)  #log
# pcs, explained_variance2 = prcomps2(translognormcounts)
# pcs = scVI.prcomps(log.(normalization(tf_data) .+ 1))
elapsed_time_seconds = @elapsed begin
        pcs = scVI.prcomps(translognormcounts)
end

# firstpcs = pcs[:,1:100]

# elapsed_time_seconds = @elapsed begin
#         Random.seed!(big_TSNE_args.seed) 
#         tSNE = tsne(rescale(tf_data, dims=1), 2, big_TSNE_args.pca_dims, 1000, big_TSNE_args.perplexity, progress=true); 
# end
# minutes, seconds = divrem(elapsed_time_seconds, 60)
# println("Execution time tSNE with $(big_TSNE_args.pca_dims) pcs: ", minutes, " minutes, ", seconds, " seconds")


elapsed_time_seconds = @elapsed begin
        Random.seed!(big_UMAP_args.seed)
        umapout = umap(pcs[:,1:big_UMAP_args.pca_dims]', min_dist=0.5, n_neighbors = big_UMAP_args.n_nbrs) 
        # umapout_pcs_10 = umap(pcs[:,1:10]', min_dist=0.5) 
        # umapout_pcs_15 = umap(pcs[:,1:15]', min_dist=0.5) 
        # umapout = umap(Transpose(translognormcounts), min_dist=0.5)
        # umapout_origdata = umap(tf_data')
end
minutes, seconds = divrem(elapsed_time_seconds, 60)
println("Execution time UMAP with $(big_UMAP_args.pca_dims) pcs: ", minutes, " minutes, ", seconds, " seconds")

df_latent_tf_data.UMAP1 = umapout[1,:]
df_latent_tf_data.UMAP2 = umapout[2,:]

plot_transformed_data = df_latent_tf_data |>  
@vlplot(width=500,
        height=400,
        :circle, 
        x={:UMAP1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        y={:UMAP2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        column = {:"timepoint_name:n",  header={title = " ", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
        color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
        #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
        size={value=30},
        labelFontSize=50,
        config={legend={titleFontSize=20, labelFontSize=20}}
)
plot_transformed_data_file = joinpath(save_folder, "original_umap_$(picked_model)_$(transnum).png")
VegaLite.save(plot_transformed_data_file, plot_transformed_data)
