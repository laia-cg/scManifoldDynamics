
# ===========================================
# Utility functions
# ===========================================


# -------------------------------------------
# Preprocessing functions
# -------------------------------------------

"""
    standardize(mat::AbstractMatrix{<:Number})

Standardizes the columns of a matrix to have zero mean and unit variance.

# Arguments
- `mat::AbstractMatrix{<:Number}`: Input data matrix.

# Returns
- A standardized matrix with zero mean and unit variance for each column.
"""
function standardize(mat::AbstractMatrix{<:Number})
    (mat .- mean(mat, dims = 1)) ./ std(mat, dims = 1)
end


"""
    rescale(mat::AbstractMatrix{<:Number}; dims::Int=1)

Rescales the rows or columns of a matrix to range between 0 and 1.

# Arguments
- `mat::AbstractMatrix{<:Number}`: Input data matrix.
- `dims::Int`: Dimension along which to rescale (1 for columns, 2 for rows, default: `1`).

# Returns
- A rescaled matrix with values ranging between 0 and 1.
"""
function rescale(mat::AbstractMatrix{<:Number}; dims::Int=1)
    (mat .- mean(mat, dims=dims)) ./ max.(std(mat, dims=dims), eps())
end



"""
    normalization(data::AbstractMatrix{<:Number}, scaleval::Number=1)

Normalizes the input data using library size normalization.

# Arguments
- `data::AbstractMatrix{<:Number}`: Input data matrix.
- `scaleval::Number`: Scaling factor for normalization (default: `1`).

# Returns
- A normalized data matrix.
"""
function normalization(data::AbstractMatrix{<:Number}, scaleval::Number=1)
    libsizes = vec(sum(data, dims=2))
    sizefactors = libsizes./mean(libsizes)
    datan = (data./sizefactors)*scaleval
    return datan
end



"""
    logNormCounts(data::AbstractMatrix{<:Number}, logbase::Number)

Applies library size normalization and log transformation to the data.

# Arguments
- `data::AbstractMatrix{<:Number}`: Input data matrix.
- `logbase::Number`: Base for the logarithm.

# Returns
- A log-normalized data matrix.
"""
function logNormCounts(data::AbstractMatrix{<:Number}, logbase::Number)
    datan = normalization(data, 1)
    datan = log.(logbase, datan .+1)
    return datan
end


int(x) = floor(Int, x)

# -------------------------------------------
# Data management functions
# -------------------------------------------


"""
    get_data(dataset::String)

Loads a dataset and initializes reduction arguments for dimensionality reduction.

# Arguments
- `dataset::String`: Name of the dataset to load.

# Returns
- An `AnnData` object containing the dataset and its metadata.
"""
function get_data(dataset::String)
    countmatrix_tmp = CSV.read("../Data/Datasets/$(dataset)_counts.csv", DataFrames.DataFrame)
    data = Array{Float32,2}(Array{Float32,2}(countmatrix_tmp[:,2:end])')  

    cell_annotation = CSV.read("../Data/Datasets/$(dataset)_cell_annotation.csv", DataFrames.DataFrame)
    cell_annotation = cell_annotation[:,:x]
    
    adata = AnnData(countmatrix=data, 
            celltypes = cell_annotation
    )
    adata.uns = Dict()
    if dataset == "PBMC"
        TSNE_args = TSNEArgs(seed = 1013, perplexity = 50, pca_dims = 35)
        UMAP_args = UMAPArgs(seed = 105, pca_dims = 35)
    elseif dataset == "Zeisel"
        TSNE_args = TSNEArgs(seed = 1013, perplexity = 20, pca_dims = 100)
        UMAP_args = UMAPArgs(seed = 9711, pca_dims = 100)
    end
    reduction_args = ReductionArgs(tsne = TSNE_args, umap = UMAP_args) 
    adata.uns["reduction_args"] = reduction_args
    return adata
end






"""
    find_clusters(adata::AnnData)

Identifies clusters based on cell annotations.

# Arguments
- `adata::AnnData`: Input data object with cell annotations.

# Returns
- `c::Vector{Vector{Int}}`: Indices of cells grouped by cluster.
- `combined_labels::Vector{String}`: Combined cell and cluster labels.
"""
function find_clusters(adata)
    cell_annotation = get_celltypes(adata)
    unique_clusters = unique(cell_annotation)
    c = Vector{Vector{Int}}(undef, length(unique_clusters))

    # Group cells by cluster
    for (i, cluster) in enumerate(unique_clusters)
        c[i] = findall(x -> x == cluster, cell_annotation)
    end
    display("$(length(unique_clusters)) different cell types have been found:")

    # We want to map cell index to cluster number, so we have the correspondence between cluster number and name
    cell_cluster_map = Dict{Int, Int}() 

    for (cluster_idx, cells) in enumerate(c)
        for cell in cells
            cell_cluster_map[cell] = cluster_idx
        end
    end
    combined_labels = [string(get_celltypes(adata)[i], " - c", cell_cluster_map[i]) for i in 1:length(get_celltypes(adata))]
    return c, combined_labels
end



"""
    create_clusters(selected_indices_vectors::Vector{Vector{Int}}, adata::AnnData)

Creates clusters based on selected indices and assigns remaining cells to an additional cluster.

# Arguments
- `selected_indices_vectors::Vector{Vector{Int}}`: Selected cell indices for initial clusters.
- `adata::AnnData`: Input data object.

# Returns
- `c::Vector{Vector{Int}}`: Indices of cells grouped by cluster.
- `combined_labels::Vector{String}`: Combined cell and cluster labels.
"""
function create_clusters(selected_indices_vectors::Vector{Vector{Int}}, adata)
    num_clusters = length(selected_indices_vectors)
    c = Vector{Vector{Int}}(undef, num_clusters + 1)  # +1 for the additional cluster
    
    for i in 1:num_clusters
        c[i] = selected_indices_vectors[i]
    end

    # Initialize the additional cluster for the remaining cells
    c[end] = []  # This is the additional cluster (last cluster in c)

    cell_cluster_map = Dict{Int, Int}()
    
    for (cluster_idx, cells) in enumerate(c[1:end-1])  # Exclude the additional cluster for this step
        for cell in cells
            cell_cluster_map[cell] = cluster_idx
        end
    end

    # Assign remaining cells to the additional cluster
    for cell in 1:size(adata.countmatrix,1)
        if !haskey(cell_cluster_map, cell)
            push!(c[end], cell)  # Add the cell to the additional cluster
            cell_cluster_map[cell] = num_clusters + 1  # Map the cell to the additional cluster
        end
    end
    
    combined_labels = [string("c", cell_cluster_map[cell]) for cell in 1:size(adata.countmatrix,1)]
    
    return c, combined_labels
end




"""
    get_celltypes(a::AnnData)

Retrieves cell type annotations from an `AnnData` object.

# Arguments
- `a::AnnData`: Input data object.

# Returns
- A vector of cell type annotations, or `nothing` if none are found.
"""
function get_celltypes(a::AnnData)
    celltypes = nothing
    if !isnothing(a.obs)
        if hasproperty(a.obs, :cell_type)
            celltypes = a.obs.cell_type
        elseif hasproperty(a.obs, :celltype)
            celltypes = a.obs.celltype
        elseif hasproperty(a.obs, :celltypes)
            celltypes = a.obs.celltypes
        elseif hasproperty(a.obs, :cell_types)
            celltypes = a.obs.cell_types
        end
    end

    if isnothing(celltypes)
        if hasproperty(a, :celltypes)
            celltypes = a.celltypes
        elseif hasproperty(a, :cell_types)
            celltypes = a.cell_types
        end
    end

    if isnothing(celltypes)
        println("No celltypes found")
    end

    return celltypes
end


# -------------------------------------------
# Models utility functions
# -------------------------------------------

"""
    load_trained_model(dataset::String, model::String)

Loads a trained model from file.

# Arguments
- `dataset::String`: Dataset name.
- `model::String`: Model name.

# Returns
- `m`: The trained model.
- `trained_data`: The associated training data.
"""
function load_trained_model(dataset::String, model::String)
    trained_model = load("../Data/Trained_models/$(dataset)_$(model).jld2")["$(dataset)_$(model)"]
    m = trained_model["trained_model"]
    trained_data = trained_model["trained_data"]
    return m, trained_data
end


"""
    save_trained_model(dataset::String, model::String, m, latent_data, latent_from_generated; original_labels=nothing)

Saves a trained model and its associated data to file.

# Arguments
- `dataset::String`: Dataset name.
- `model::String`: Model name.
- `m`: Trained model.
- `latent_data`: Latent representation of the data.
- `generated_data`: Generated data by the decoder.
- `original_labels`: Optional original labels for the data.

# Returns
- Saves the model to the specified file path.
"""
function save_trained_model(
    dataset::String, 
    model::String, 
    m::Any, 
    latent_data::AbstractArray{<:Number}, 
    generated_data::AbstractArray{<:Number}; 
    original_labels::Union{AbstractArray{<:Number}, Nothing}=nothing
    )
    
    save_directory = "../Data/Trained_models/"
    if !isdir(save_directory)
        mkpath(save_directory)
        @info "Created directory: $save_directory"
    else
        @info "Directory already exists: $save_directory"
    end

    if model == "scvi"
        trained_data = Dict(
            "latent" => latent_data,
            "generated" => generated_data
        )
    else
        trained_data = Dict(
        "original_labels" => original_labels,
        "latent" => latent_data,
        "generated" => generated_data
    )
    end
    save_model = Dict("trained_model" => m, "trained_data"=> trained_data)
    save_path = joinpath(save_directory, "$(dataset)_$(model).jld2")
    
    JLD2.save(save_path, "$(dataset)_$(model)", save_model)
    @info "Model saved at: $save_path"
end






# function get_celltypes(a::AnnData)
#     return a.celltypes
# end

# function save_trained_model(dataset::String, model::String, m, original_labels, latent_data, latent_from_generated)
#     trained_data = Dict(
#         "original_labels" => original_labels,
#         "latent" => latent_data,
#         "generated" => latent_from_generated
#     )
    
#     save_model = Dict("trained_model" => m, "trained_data"=> trained_data)
#     save_path = "../Data/Trained_models/$(dataset)_$(model).jld2"
    
#     JLD2.save(save_path, "$(dataset)_$(model)", save_model)
# end