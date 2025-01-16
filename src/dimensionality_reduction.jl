
# ===========================================
# Functions related to dimensionality reduction methods
# ===========================================



# -------------------------------------------
# Structures to encode parameters related to dimensionality reduction methods
# -------------------------------------------

"""
    @with_kw struct TSNEArgs

Arguments for configuring t-SNE dimensionality reduction.

# Fields
- `seed::Int`: Random seed for reproducibility (default: `1013`).
- `perplexity::Int`: Perplexity parameter for t-SNE (default: `20`).
- `pca_dims::Int`: Number of PCA dimensions to reduce before applying t-SNE (default: `100`).
"""
@with_kw struct TSNEArgs
    seed::Int = 1013
    perplexity::Int = 20
    pca_dims::Int = 100
end

"""
    @with_kw struct UMAPArgs

Arguments for configuring UMAP dimensionality reduction.

# Fields
- `seed::Int`: Random seed for reproducibility (default: `9711`).
- `pca_dims::Int`: Number of PCA dimensions to reduce before applying UMAP (default: `100`).
"""
@with_kw struct UMAPArgs
    seed::Int = 9711
    pca_dims::Int = 100
    min_dist::Number = 0.5
end

"""
    @with_kw struct ReductionArgs

A container for arguments related to dimensionality reduction methods.

# Fields
- `tsne::TSNEArgs`: Arguments for t-SNE.
- `umap::UMAPArgs`: Arguments for UMAP.
"""
@with_kw struct ReductionArgs
    tsne::TSNEArgs = TSNEArgs()
    umap::UMAPArgs = UMAPArgs()
end

# -------------------------------------------
# Dimensionality reduction functions
# -------------------------------------------

"""
    prcomps(mat::AbstractMatrix{<:Number}, standardizeinput::Bool=true)

Performs Principal Component Analysis (PCA) on the given matrix.

# Arguments
- `mat::AbstractMatrix{<:Number}`: Input data matrix.
- `standardizeinput::Bool`: Whether to standardize the input data before PCA (default: `true`).

# Returns
- `prcomps`: Principal components of the input matrix.
- `explained_variance_ratio`: Vector of variance ratios explained by each principal component.
"""
function prcomps(mat::AbstractMatrix{<:Number}, standardizeinput::Bool=true)
    if standardizeinput
        mat = standardize(mat)
    end
    u,s,v = svd(mat)
    prcomps = u * Diagonal(s)

    # Obtain percentage of variance explained
    explained_variance = (s.^2) / (size(mat, 1) - 1)
    explained_variance_ratio = explained_variance / sum(explained_variance)

    return prcomps, explained_variance_ratio
end


"""
    get_reduced_dimensions(adata; reduction_args=nothing, tSNE_seed::Int=1013, UMAP_seed::Int=9711)

Performs dimensionality reduction on the dataset using PCA, t-SNE, and UMAP.

# Arguments
- `adata`: An AnnData object containing the dataset.
- `reduction_args::Union{ReductionArgs, Nothing}`: Arguments for dimensionality reduction methods (default: `nothing`, which takes the default parameters).


# Returns
- `df_stacked_all_reduced::DataFrame`: DataFrame containing reduced dimensions for PCA, t-SNE, and UMAP.
- `plot_reduced`: Plot visualizing reduced dimensions.
"""
function get_reduced_dimensions(
    adata; 
    reduction_args::Union{ReductionArgs, Nothing}=nothing, 
)
    
    if isnothing(reduction_args)
        reduction_args = get(adata.uns, "reduction_args", ReductionArgs())
    end
    
    #PCA
    pcs = scVI.prcomps(log.(normalization(adata.countmatrix) .+ 1))

    #UMAP
    Random.seed!(reduction_args.umap.seed) 
    umapout = umap(pcs[:,1:reduction_args.umap.pca_dims]', min_dist=reduction_args.umap.min_dist)
    ncells = size(umapout, 2)

    #tSNE
    Random.seed!(reduction_args.tsne.seed) 
    tSNE = tsne(rescale(adata.countmatrix, dims=1), 2, reduction_args.tsne.pca_dims, 1000, reduction_args.tsne.perplexity, progress=false) 

    
    df_stacked_all_reduced = DataFrame(z1 = vcat(pcs[:,1], tSNE[:,1], umapout[1,:]), 
                            z2 = vcat(pcs[:,2], tSNE[:,2], umapout[2,:]),
                            Celltype = repeat(get_celltypes(adata), outer = 3),
                            Method = vcat(fill("PCS", ncells), fill("tSNE", ncells), fill("UMAP", ncells)), 
    )

    # Plot
    plot_reduced = df_stacked_all_reduced |>
    @vlplot(width=450,
            height=300,
            :circle, 
            x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
            y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
            column = {"Method:n", header={title = nothing, labelFontSize=20, titleFontSize=15, labelFontWeight="bold"}},
            color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=false, title="Cell annotation",orient="right"}},
            size={value=25},
            config={legend={titleFontSize=20, labelFontSize=15}},
            resolve={scale={x="independent",y="independent"}}
    )
    return df_stacked_all_reduced, plot_reduced
end
