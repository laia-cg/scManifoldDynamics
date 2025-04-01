using Pkg
Pkg.activate(".")

using CSV, DataFrames, VegaLite, Plots, Distributions, Random
using UMAP, TSne
using Flux
using JLD2
using Clustering
using ProgressMeter


include("../src/scManifoldDynamics.jl") 
using .scManifoldDynamics  

save_folder = "figures/hyperparameter_evaluation/"
mkpath(save_folder)

# =======================
#   One time point 
# =======================

dataset = "Lymph"
adata = get_data(dataset)

reduction_args = get(adata.uns, "reduction_args", ReductionArgs())

#PCA
pcs = scVI.prcomps(log.(normalization(adata.countmatrix) .+ 1))

X_lognorm = log.(normalization(adata.countmatrix) .+ 1)
cell_annotation = get_celltypes(adata)

println(pwd())

# -----  study on UMAP -------

vec_nbrs = [15, 20, 30, 50]  # 15 default
vec_min_dist = [0.01, 0.1, 0.3, 0.5]  # 0.1 default
vec_pcs = [10, 20, 30, 50, 100]

total_iterations = length(vec_pcs) * length(vec_nbrs) * length(vec_min_dist)
p = Progress(total_iterations, desc="Processing UMAP parameters...")

results_umap = Dict() # Dictionary to store the results for further analysis

data_results_folder = normpath(joinpath(@__DIR__, "..", "Data", "hyperparameter_evaluation/UMAP_$(dataset)/"))
mkpath(data_results_folder)

figures_results_folder = normpath(joinpath(@__DIR__, "..", "figures", "hyperparameter_evaluation/UMAP_$(dataset)/"))
mkpath(figures_results_folder)

for n_pcs in vec_pcs
    for n_nbrs in vec_nbrs
        for min_dist in vec_min_dist
            data = n_pcs == 0 ? X_lognorm' : pcs[:, 1:n_pcs]'

            @info "Running UMAP with pcs = $n_pcs, n_neighbors=$n_nbrs, min_dist=$min_dist "
            
            # Run UMAP with different parameter configurations
            Random.seed!(reduction_args.umap.seed) 
            umapout = umap(data, 
                        min_dist=min_dist, 
                        n_neighbors=n_nbrs)


            df_stacked_umap = DataFrame(z1 = umapout[1, :], 
                                        z2 = umapout[2, :],
                                        Celltype = cell_annotation)

            # Plot UMAP embedding
            plot_reduced = df_stacked_umap |>
            @vlplot(width=350,
                    height=250,
                    :circle, 
                    x={:z1, title=" ", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
                    y={:z2, title=" ", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
                    # color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Cell annotation",orient="right"}},
                    color={
                        "Celltype:n", 
                        # scale={scheme="tableau10"}, 
                        scale={range=tableau20_extended}, 
                        legend={disable = true, title="Cell annotation", orient="right", direction="vertical"}
                    },
                    size={value=5},
                    config={legend={titleFontSize=20, labelFontSize=15}},
                    resolve={scale={x="independent", y="independent"}})

  
            plot_file = joinpath(figures_results_folder, "UMAP_n_pcs_$(n_pcs)_nbrs_$(n_nbrs)_min_dist_$(min_dist).svg")
            VegaLite.save(plot_file, plot_reduced)

            umap_coords = transpose(umapout)  # Transpose for compatibility with Clustering.jl

            # # Compute Silhouette score
            # dist_matrix = pairwise(Euclidean(), umap_coords)
            # silhouette_score = mean(silhouettes(predicted_labels, dist_matrix))

            # Store results
            param_key = "n_pcs_$(n_pcs)_nbrs_$(n_nbrs)_min_dist_$(min_dist)"
            results_umap[param_key] = Dict( "embedding" => umap_coords)               
            next!(p)                
        end
    end
end

@save joinpath(data_results_folder, "$(dataset)_umap_results.jld2") results_umap
println("UMAP analysis completed. Results saved to $data_results_folder.")


# ----- study on tSNE ---------

data_results_folder = base_dir = normpath(joinpath(@__DIR__, "..", "Data", "hyperparameter_evaluation/tSNE_$(dataset)/"))
mkpath(data_results_folder)

figures_results_folder = normpath(joinpath(@__DIR__, "..", "figures", "hyperparameter_evaluation/tSNE_$(dataset)/"))
mkpath(figures_results_folder)

# Define parameter ranges for tSNE
vec_perplexity = [15, 30, 50]  # 15 default
vec_pcs = [10, 20, 30, 40, 100]

total_iterations = length(vec_perplexity) * length(vec_pcs) 
p = Progress(total_iterations, desc="Processing tSNE parameters...")

results_tsne = Dict() # Dictionary to store the results for further analysis

# Iterate over combinations of parameters
for n_pcs in vec_pcs
    for per in vec_perplexity
        next!(p)                
        @info "Running tSNE with pcs = $n_pcs, perplexity=$per "
        
        # Run tSNE with different parameter configurations
        Random.seed!(reduction_args.tsne.seed) 
        tSNE = tsne(rescale(adata.countmatrix, dims=1), 2, n_pcs, 1000, per, progress=true) 

        df_stacked_tsne = DataFrame(z1 = tSNE[:, 1], 
                                    z2 = tSNE[:, 2],
                                    Celltype = cell_annotation)

        # Plot tSNE embedding
        plot_reduced = df_stacked_tsne |>
        @vlplot(width=350,
                height=250,
                :circle, 
                x={:z1, title=" ", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
                y={:z2, title=" ", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
                # color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Cell annotation",orient="right"}},
                color={"Celltype:n", scale={range=tableau20_extended}, legend={disable=true, title="Cell annotation",orient="right"}},
                size={value=5},
                config={legend={titleFontSize=20, labelFontSize=15}},
                resolve={scale={x="independent", y="independent"}})

        plot_file = joinpath(figures_results_folder, "n_pcs_$(n_pcs)_tsne_perplexity_$(per).svg")
        VegaLite.save(plot_file, plot_reduced)

        # # Compute Silhouette score
        # dist_matrix = pairwise(Euclidean(), umap_coords)
        # silhouette_score = mean(silhouettes(predicted_labels, dist_matrix))

        # Store results
        param_key = "perplexity_$(per)_n_pcs_$(n_pcs)"
        results_tsne[param_key] = Dict("embedding" => tSNE)               
    end
end

@save joinpath(data_results_folder, "$(dataset)_tsne_results.jld2") results_tsne
println("t-SNE analysis completed. Results saved to $data_results_folder.")




# ==============================
# Further analysis 
# ==============================

results_dict = results_tsne
# for (param_key, param_dict) in results_dict
#     println("$(param_key)")
#     println(param_dict["embedding"][1:2,:])
# end

num_points = size(embedding, 1)

best_k = 0
best_ari = -1
k_values = [10]
ari_scores = []
resolution_leiden = 0.0005

for k in k_values
    # Construct kNN graph as adjacency matrix
    dist_matrix = pairwise(Euclidean(), embedding, dims=1)
    A = spzeros(num_points, num_points)

    for i in 1:num_points
        neighbors = sortperm(dist_matrix[i, :])[2:k+1]
        for j in neighbors
            A[i, j] = 1
            A[j, i] = 1
        end
    end

    # Run Leiden clustering
    Random.seed!(1234)
    leiden_result = Leiden.leiden(A, resolution=resolution_leiden)
    
    # Extract labels
    predicted_labels = zeros(Int, num_points)
    for (cluster_id, cluster_points) in enumerate(leiden_result.partition)
        for point in cluster_points
            predicted_labels[point] = cluster_id
        end
    end

    # Compute ARI
    ari_score = randindex(integer_labels, predicted_labels)
    push!(ari_scores, ari_score)

    df_scatter_transformed[!,"predicted_labels"] = predicted_labels

    plot_transformed_data = df_scatter_transformed |>
    @vlplot(
        width=400,
        height=300,
        :circle,
        x = { :z1, title = " ", axis = { titleFontSize = 10, labelFontSize = 10, tickCount = 5 } },
        y = { :z2, title = " ", axis = { titleFontSize = 10, labelFontSize = 10, tickCount = 5 } },
        # color = { "Celltype:n", scale = { scheme = "tableau10" }, legend = { disable = true, title = "Cell annotation", orient = "bottom", direction = "horizontal" } },
        color = {"predicted_labels:n", scale = { scheme = "tableau10" }, legend = { disable = false, title = "Cell annotation", orient = "right", direction = "vertical" } },
        size = { value = 20 },
        labelFontSize = 50,
        title = "$(param_key)_$(silhouette_score)",
        config = { legend = { titleFontSize = 20, labelFontSize = 20 } }
    )

    # Track best k
    if ari_score > best_ari
        best_ari = ari_score
        best_k = k
    end
end

# Plot ARI vs. k
plot(k_values, ari_scores, xlabel="k", ylabel="ARI Score", title="Choosing k for Leiden", lw=2)
println("Best k based on ARI: ", best_k)

# =======================
#  Multiple time points
# =======================


# -------- GET ALL BASE VARIABLES READY --------------
picked_dataset = "Lymph"
picked_model = "tsne_sae"
transnum = "t1" 


folderdata = "Data/Synthetic_data/$(picked_dataset)/$(picked_model)/"

folderresults = "Data/Synthetic_data/$(picked_dataset)/$(picked_model)/Results/hyperparameter_evaluation/"
mkpath(folderresults)

cell_annotation = CSV.read("Data/Datasets/$(picked_dataset)_cell_annotation.csv", DataFrames.DataFrame)
cell_annotation = cell_annotation[:,:x]



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

# --------- Look at the loaded transformation (if needed) ---------------

check_data = latent_tf_data  #umapout' #umapout' # latent_tf_data #umapout_pcs_10'

plot_original_latent = true
if plot_original_latent
        df_scatter_transformed = DataFrames.DataFrame(z1 = check_data[:,1], 
                                                z2 = check_data[:,2], 
                                                timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
                                                Celltype  = repeat(cell_annotation, outer=ntimepoints), 
        )


        plot_transformed_data = df_scatter_transformed |>  
        @vlplot(width=500,
                height=400,
                :circle, 
                x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
                y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
                column = {:"timepoint:n",  header={title = "Timepoint", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
                color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
                #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
                size={value=30},
                labelFontSize=50,
                config={legend={titleFontSize=20, labelFontSize=20}}
        )
        # VegaLite.save(folderexperiments*"$(picked_dataset)_TIME_$(check).png", plot_transformed_data)
end


# ---------- Explained variance -----------

pcs = scVI.prcomps(translognormcounts)

df_pcs = DataFrames.DataFrame(z1 = pcs[:,1], 
                                                z2 = pcs[:,2], 
                                                timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
                                                Celltype  = repeat(cell_annotation, outer=ntimepoints), 
        )


plot_transformed_data = df_pcs |>  
@vlplot(width=500,
        height=400,
        :circle, 
        x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        column = {:"timepoint:n",  header={title = "Timepoint", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
        color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
        #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
        size={value=30},
        labelFontSize=50,
        config={legend={titleFontSize=20, labelFontSize=20}}
)


# ---------- Parameters for the models -----------


if picked_dataset == "PBMC"
    big_TSNE_args = TSNEArgs(seed = 1013, perplexity = 200, pca_dims = 20)
    big_UMAP_args = UMAPArgs(seed = 87, pca_dims = 10) #105
elseif picked_dataset == "Zeisel"
    big_TSNE_args = TSNEArgs(seed = 1013, perplexity = 100, pca_dims = 20)
    big_UMAP_args = UMAPArgs(seed = 9711, pca_dims = 100)
elseif picked_dataset == "Lymph"
    big_TSNE_args = TSNEArgs(seed = 1013, perplexity = 50, pca_dims = 50)
    big_UMAP_args = UMAPArgs(seed = 9711, pca_dims = 50,  n_nbrs = 30)
end


# -----  study on UMAP -------
using ProgressMeter

vec_nbrs = [15, 30, 50, 70]  # 15 default
vec_min_dist = [0.1, 0.3, 0.5]  # 0.1 default
vec_pcs = [10, 20, 30, 40, 100]

total_iterations = length(vec_pcs) * length(vec_nbrs) * length(vec_min_dist)
p = Progress(total_iterations, desc="Processing UMAP parameters...")

results_umap = Dict()

# Iterate over combinations of parameters
for n_pcs in vec_pcs
    for n_nbrs in vec_nbrs
        for min_dist in vec_min_dist
            data = n_pcs == 0 ? translognormcounts' : pcs[:, 1:n_pcs]'

            @info "Running UMAP with pcs = $n_pcs, n_neighbors=$n_nbrs, min_dist=$min_dist "
            
            # Run UMAP with specified parameters
            Random.seed!(big_UMAP_args.seed)
            umapout = umap(data, 
                        min_dist=min_dist, 
                        n_neighbors=n_nbrs)

            df_stacked_umap = DataFrame(z1 = umapout[1, :], 
                                        z2 = umapout[2, :],
                                        timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
                                        Celltype  = repeat(cell_annotation, outer=ntimepoints))

            plot_reduced = df_stacked_umap |>
            @vlplot(width=500,
                    height=400,
                    :circle, 
                    x={:z1, title=" ", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
                    y={:z2, title=" ", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
                    column = {:"timepoint:n",  header={title = " ", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}},
                    color={
                        "Celltype:n", 
                        scale={range=tableau20_extended}, 
                        legend={disable = true, title="Cell annotation", orient="bottom", direction="horizontal"}
                    },
                    size={value=10},
                    config={legend={titleFontSize=20, labelFontSize=15}},
                    resolve={scale={x="independent", y="independent"}})


            plot_file = joinpath(save_folder, "UMAP_n_pcs_$(n_pcs)_nbrs_$(n_nbrs)_min_dist_$(min_dist).png")
            # save(plot_file, plot_reduced)
            VegaLite.save(plot_file, plot_reduced)

            param_key = "n_nbrs_$(n_nbrs)_min_dist_$(min_dist)_npcs_$(n_pcs)"
            results_umap[param_key] = Dict( "embedding" => umapout)               
            next!(p)                
        end
    end
end

@save joinpath(save_folder, "$(picked_dataset)_umap_results.jld2") results_umap



# ----- study on tSNE ---------
save_folder = "/Users/laiacanal/Documents/TimeSeries_project/Revisions/hyperparameter_justification/$(picked_dataset)_$(transnum)_$(picked_model)_tSNE"

# Define parameter ranges for tSNE
vec_perplexity = [30, 100, 200]  # 15 default
vec_pcs = [10, 20, 50, 100]

total_iterations = length(vec_pcs) * length(vec_perplexity) 
p = Progress(total_iterations, desc="Processing tSNE parameters...")

# Initialize random seed
Random.seed!(1013)



# Prepare to store results
results_tsne = Dict()

# Iterate over combinations of parameters
for n_pcs in vec_pcs
    for per in vec_perplexity
        next!(p)                
        @info "Running tSNE with pcs = $n_pcs, perplexity=$per "
        
        # Run tSNE with specified parameters
        Random.seed!(1013) 
        tSNE = tsne(rescale(tf_data, dims=1), 2, n_pcs, 1000, per, progress=true) 

        df_scatter_transformed = DataFrames.DataFrame(z1 = tSNE[:,1], 
            z2 = tSNE[:,2], 
            timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
            Celltype  = repeat(cell_annotation, outer=ntimepoints), 
        )


        plot_transformed_data = df_scatter_transformed |>  
        @vlplot(width=500,
            height=400,
            :circle, 
            x={:z1, title=" ", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
            y={:z2, title=" ", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
            column = {:"timepoint:n",  header={title = " ", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
            # color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
            color={
                "Celltype:n", 
                scale={scheme="tableau10"}, 
                legend={disable = true, title="Cell annotation", orient="bottom", direction="horizontal"}
            },
            #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
            size={value=30},
            labelFontSize=50,
            config={legend={titleFontSize=20, labelFontSize=20}}
        )

        plot_file = joinpath(save_folder, "n_pcs_$(n_pcs)_tsne_perplexity_$(per).svg")
        VegaLite.save(plot_file, plot_transformed_data)

        # # Compute Silhouette score
        # dist_matrix = pairwise(Euclidean(), umap_coords)
        # silhouette_score = mean(silhouettes(predicted_labels, dist_matrix))

        # Store results
        param_key = "perplexity_$(per)_n_pcs_$(n_pcs)"
        results_tsne[param_key] = Dict("embedding" => tSNE)               
    end
end

@save joinpath(save_folder, "$(picked_dataset)_tsne_results.jld2") results_tsne
