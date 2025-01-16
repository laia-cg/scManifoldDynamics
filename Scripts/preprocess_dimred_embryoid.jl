#-------------------------------------------------------------------------------------------------
# setup 
#-------------------------------------------------------------------------------------------------

using Pkg;
Pkg.activate(".")

using DataFrames
using HDF5
using SparseArrays
using StatsBase
using LinearAlgebra
#Pkg.add(url="https://github.com/maren-ha/scVI.jl")
using scVI 

#-------------------------------------------------------------------------------------------------
# functions
#-------------------------------------------------------------------------------------------------

function standardize(x)
    (x .- mean(x, dims = 1)) ./ std(x, dims = 1)
end

function prcomps(mat, standardizeinput = true)
    if standardizeinput
        mat = standardize(mat)
    end
    u,s,v = svd(mat)
    prcomps = u * Diagonal(s)
    return prcomps
end

open_h5_data(filename::String; mode::String="r+") = h5open(filename, mode)

function Base.sort!(idx::AbstractArray, vals::AbstractArray...)
    if !issorted(idx)
        ordering = sortperm(idx)
        permute!(idx, ordering)
        for v in vals
            permute!(v, ordering)
        end
    end
end

# from Muon.jl: https://github.com/scverse/Muon.jl/blob/9ed2b4cadd1173ba90f8bc2932dfe8f421972872/src/hdf5_io.jl#L58
function read_matrix(f::HDF5.Group; kwargs...)
    enctype = read_attribute(f, "encoding-type")

    if enctype == "csc_matrix" || enctype == "csr_matrix"
        shape = read_attribute(f, "shape")
        iscsr = enctype[1:3] == "csr"

        indptr = read(f, "indptr")
        indices = read(f, "indices")
        data = read(f, "data")

        indptr .+= eltype(indptr)(1)
        indices .+= eltype(indptr)(1)

        # the row indices in every column need to be sorted
        @views for (colstart, colend) in zip(indptr[1:(end - 1)], indptr[2:end])
            sort!(indices[colstart:(colend - 1)], data[colstart:(colend - 1)])
        end

        if iscsr
            reverse!(shape)
        end
        mat = SparseMatrixCSC(shape..., indptr, indices, data)
        return iscsr ? mat' : mat
    elseif enctype == "categorical"
        ordered = read_attribute(f, "ordered") > 0
        categories = read(f, "categories")
        codes = read(f, "codes") .+ 1

        T = any(codes .== 0) ? Union{Missing, eltype(categories)} : eltype(categories)
        mat = CategoricalVector{T}(
            undef, length(codes); levels=categories, ordered=ordered)
        copy!(mat.refs, codes)

        return mat
    else
        error("unknown encoding $enctype")
    end
end

function rescale(A; dims=1)
    (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())
end

function normalization(data, scaleval)
    libsizes = vec(sum(data, dims=2))
    #datan = (data ./libsizes)*scaleval 
    sizefactors = libsizes./mean(libsizes)
    datan = (data./sizefactors)*scaleval
    return datan
end

function logNormCounts(data, logbase) # log-transformed normalized (taking into consideration the sizefactors ) expression values from a count matrix. It's the same as logNormCounts in R-scuttle package
    datan = normalization(data, 1)
    datan = log.(logbase, datan .+1)
    return datan
end

#-------------------------------------------------------------------------------------------------
# get the data 
#-------------------------------------------------------------------------------------------------

filename = "embryoid_anndata_complete.h5ad"
file = open_h5_data(filename)

X = read(file, "X")' # shape: cell x gene 
counts = read_matrix(file["layers"]["counts"])
norm = read_matrix(file["layers"]["norm"])
sqrt = read_matrix(file["layers"]["sqrt"])

layers = read(file, "layers")
layers = file["layers"]
# transform obs to dataframe
obs = read(file, "obs")
obs_df = DataFrame(obs)
# sample_labels 
sample_labels = obs["sample_labels"]["categories"][Int.(obs["sample_labels"]["codes"]).+1 ]
obs_df[!,:sample_labels] = sample_labels
obs_df[!,:sample_labels_codes] = Int.(obs["sample_labels"]["codes"]).+1
# check for leftover Dicts
any(x -> x == Dict{String, Any}, eltype.(eachcol(obs_df))) && @warn "There are still dictionary types in the columns of the obs df, please check!"
# uns
uns = read(file, "uns")
# transform var to dataframe 
var = read(file, "var")
var_df = DataFrame(var)
# eltype.(eachcol(var_df)) # seems fine 

adata = AnnData(countmatrix = counts, 
                layers = Dict("counts" => counts, 
                            "norm" => norm, 
                            "sqrt" => sqrt),
                obs = obs_df,
                var = var_df
)

rename(adata.var, "highly_variable" => "highly_variable_python")
highly_variable_genes!(adata; n_top_genes=1000, replace_hvgs=true, verbose=true)
subset_to_hvg!(adata)

#hvgs = adata.var[!,:highly_variable_python]
#sum(hvgs)
#@assert size(adata.countmatrix,2) == length(hvgs)
#adata.countmatrix = adata.countmatrix[:,hvgs]
#adata.var = adata.var[hvgs,:]
#@assert sum(adata.var[!,:highly_variable]) == size(adata.countmatrix,2)
#@assert !any(isnan.(adata.var[!,:highly_variable_rank]))

X = Float32.(adata.countmatrix)
lognormcounts = logNormCounts(X, 2)

using TSne
using UMAP
using Random 
using CSV

# PCA
pcs = scVI.prcomps(lognormcounts)
CSV.write("dimension_reduction_results/PCs.csv", DataFrame(pcs, :auto), delim='\t')

# tSNE on rescaled only 
Random.seed!(1013)  #1234   #3281 (not the best)
tSNE_rescaledonly = tsne(rescale(X, dims=1), 2, 50, 1000, 20.0)
CSV.write("dimension_reduction_results/tSNE_rescaleddata.csv", DataFrame(tSNE_rescaledonly, :auto), delim = '\t')

# tSNE on rescaled lognormcounts 
Random.seed!(1013)  #1234   #3281 (not the best)
tSNE_rescaledata = tsne(rescale(lognormcounts, dims=1), 2, 50, 1000, 20.0)
CSV.write("dimension_reduction_results/tSNE_rescaledlognormcounts.csv", DataFrame(tSNE_rescaledata, :auto), delim = '\t')

# tSNE on lognormcounts 
Random.seed!(1013)  #1234   #3281 (not the best)
tSNE_lognormcounts = tsne(lognormcounts, 2, 50, 1000, 20.0)
CSV.write("dimension_reduction_results/tSNE_lognormcounts.csv", DataFrame(tSNE_lognormcounts, :auto), delim = '\t')

# tSNE on first 100 PCs # takes forever
Random.seed!(1013)  #1234   #3281 (not the best)
tSNE_100pcs = tsne(pcs[:,1:100], 2, 50, 1000, 20.0)
CSV.write("dimension_reduction_results/tSNE_100pcs.csv", DataFrame(tSNE_100pcs, :auto), delim = '\t')

# UMAP 
Random.seed!(3821) #1234 #74
umapout_lognormcounts = umap(Transpose(lognormcounts))
CSV.write("dimension_reduction_results/umap_lognormcounts.csv", DataFrame(umapout_lognormcounts, :auto), delim = '\t')

Random.seed!(3821) #1234 #74
umapout_rescaled = umap(Transpose(rescale(adata.countmatrix)))
CSV.write("dimension_reduction_results/umap_rescaled.csv", DataFrame(umapout_rescaled, :auto), delim = '\t')

Random.seed!(3821) #1234 #74
umapout_100pcs = umap(Transpose(pcs[:,1:100]))
CSV.write("dimension_reduction_results/umap_100pcs.csv", DataFrame(umapout_100pcs, :auto), delim = '\t')

# scVI 

adata.countmatrix = Float32.(adata.countmatrix)
m = scVAE(size(adata.countmatrix,2);
    n_latent=2,
    #n_layers=2
)
training_args = TrainingArgs(
    max_epochs=100, # 50 for 10-dim 
    weight_decay=Float32(1e-6)
)
train_model!(m, adata, training_args)
latent_scvi = get_latent_representation(m, adata.countmatrix)
CSV.write("dimension_reduction_results/scvi.csv", DataFrame(latent_scvi, :auto), delim ='\t')

#
#
#

df_all = DataFrame(PC1 = pcs[:,1], PC2 = pcs[:,2], 
                tSNE_resc1 = tSNE_rescaledonly[:,1], tSNE_resc2 = tSNE_rescaledonly[:,2],
                tSNE_lognorm1 = tSNE_lognormcounts[:,1], tSNE_lognorm2 = tSNE_lognormcounts[:,2],
                tSNE_resclognorm1 = tSNE_rescaledata[:,1], tSNE_resclognorm2 = tSNE_rescaledata[:,2],
                tSNE_100pcs1 = tSNE_100pcs[:,1], tSNE_100pcs2 = tSNE_100pcs[:,2],
                UMAP_lognorm1 = umapout_lognormcounts[1,:], UMAP_lognorm2 = umapout_lognormcounts[2,:],
                UMAP_resc1 = umapout_rescaled[1,:], UMAP_resc2 = umapout_rescaled[2,:],
                UMAP_100pcs1 = umapout_100pcs[1,:], UMAP_100pcs2 = umapout_100pcs[2,:],
                scVI_z1 = latent_scvi[1,:], scVI_z2 = latent_scvi[2,:]
)
CSV.write("dimension_reduction_results/df_all_methods.csv", df_all, delim='\t')

dict_all = Dict("pcs" => pcs, 
                "tSNE_rescaled" => tSNE_rescaledonly,
                "tSNE_rescaledlognormcounts" => tSNE_rescaledata, 
                "tSNE_lognormcounts" => tSNE_lognormcounts,
                "tSNE_100pcs" => tSNE_100pcs,
                "umapout_lognormcounts" => umapout_lognormcounts,
                "umapout_rescaled" => umapout_rescaled,
                "UMAP_100pcs" => umapout_100pcs,
                "latent_scvi" => latent_scvi
)

using JLD2
@save "dimension_reduction_results/dictionary_with_all_methods.jld2" dict_all

# Plot all together
ncells, ngenes = size(adata.countmatrix)
df_stacked_all_reduced = DataFrame(z1 = vcat(pcs[:,1], tSNE_rescaledonly[:,1], umapout_lognormcounts[1,:], latent_scvi[1,:]), 
                        z2 = vcat(pcs[:,2], tSNE_rescaledonly[:,2], umapout_lognormcounts[2,:], latent_scvi[1,:]),
                        clusters_spectral = repeat(adata.obs[!,"clusters_spectral"], outer = 4),
                        clusters_louvain = repeat(adata.obs[!,"clusters_louvain"], outer = 4),
                        timepoint = repeat(adata.obs[!,"timepoint"], outer = 4),
                        Method = vcat(fill("PCA", ncells), fill("TSNE(rescaled)", ncells), fill("UMAP", ncells), fill("scVI", ncells)), 
)
CSV.write("dimension_reduction_results/stacked_df_for_plotting.csv", df_stacked_all_reduced, delim='\t')

#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
# repeat with 500 HVGs 
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------

adata = AnnData(countmatrix = counts, 
                layers = Dict("counts" => counts, 
                            "norm" => norm, 
                            "sqrt" => sqrt),
                obs = obs_df,
                var = var_df
)

rename(adata.var, "highly_variable" => "highly_variable_python")
highly_variable_genes!(adata; n_top_genes=500, replace_hvgs=true, verbose=true)
subset_to_hvg!(adata)

#hvgs = adata.var[!,:highly_variable_python]
#sum(hvgs)
#@assert size(adata.countmatrix,2) == length(hvgs)
#adata.countmatrix = adata.countmatrix[:,hvgs]
#adata.var = adata.var[hvgs,:]
#@assert sum(adata.var[!,:highly_variable]) == size(adata.countmatrix,2)
#@assert !any(isnan.(adata.var[!,:highly_variable_rank]))

X = Float32.(adata.countmatrix)
lognormcounts = logNormCounts(X, 2)

using TSne
using UMAP
using Random 
using CSV

# PCA
pcs = scVI.prcomps(lognormcounts)
CSV.write("dimension_reduction_results_500hvgs/PCs.csv", DataFrame(pcs, :auto), delim='\t')

# tSNE on rescaled only 
Random.seed!(1013)  #1234   #3281 (not the best)
tSNE_rescaledonly = tsne(rescale(X, dims=1), 2, 50, 1000, 20.0)
CSV.write("dimension_reduction_results_500hvgs/tSNE_rescaleddata.csv", DataFrame(tSNE_rescaledonly, :auto), delim = '\t')

# tSNE on rescaled lognormcounts 
Random.seed!(1013)  #1234   #3281 (not the best)
tSNE_rescaledata = tsne(rescale(lognormcounts, dims=1), 2, 50, 1000, 20.0)
CSV.write("dimension_reduction_results_500hvgs/tSNE_rescaledlognormcounts.csv", DataFrame(tSNE_rescaledata, :auto), delim = '\t')

# tSNE on lognormcounts 
Random.seed!(1013)  #1234   #3281 (not the best)
tSNE_lognormcounts = tsne(lognormcounts, 2, 50, 1000, 20.0)
CSV.write("dimension_reduction_results_500hvgs/tSNE_lognormcounts.csv", DataFrame(tSNE_lognormcounts, :auto), delim = '\t')

# tSNE on first 100 PCs # takes forever
Random.seed!(1013)  #1234   #3281 (not the best)
tSNE_100pcs = tsne(pcs[:,1:100], 2, 50, 1000, 20.0)
CSV.write("dimension_reduction_results_500hvgs/tSNE_100pcs.csv", DataFrame(tSNE_100pcs, :auto), delim = '\t')

# UMAP 
Random.seed!(3821) #1234 #74
umapout_lognormcounts = umap(Transpose(lognormcounts))
CSV.write("dimension_reduction_results_500hvgs/umap_lognormcounts.csv", DataFrame(umapout_lognormcounts, :auto), delim = '\t')

Random.seed!(3821) #1234 #74
umapout_rescaled = umap(Transpose(rescale(adata.countmatrix)))
CSV.write("dimension_reduction_results_500hvgs/umap_rescaled.csv", DataFrame(umapout_rescaled, :auto), delim = '\t')

Random.seed!(3821) #1234 #74
umapout_100pcs = umap(Transpose(pcs[:,1:100]))
CSV.write("dimension_reduction_results_500hvgs/umap_100pcs.csv", DataFrame(umapout_100pcs, :auto), delim = '\t')

# scVI 

adata.countmatrix = Float32.(adata.countmatrix)
m = scVAE(size(adata.countmatrix,2);
    n_latent=2,
    #n_layers=2
)
training_args = TrainingArgs(
    max_epochs=100, # 50 for 10-dim 
    weight_decay=Float32(1e-6)
)
train_model!(m, adata, training_args)
latent_scvi = get_latent_representation(m, adata.countmatrix)
CSV.write("dimension_reduction_results_500hvgs/scvi.csv", DataFrame(latent_scvi, :auto), delim ='\t')

#
#
#

df_all = DataFrame(PC1 = pcs[:,1], PC2 = pcs[:,2], 
                tSNE_resc1 = tSNE_rescaledonly[:,1], tSNE_resc2 = tSNE_rescaledonly[:,2],
                tSNE_lognorm1 = tSNE_lognormcounts[:,1], tSNE_lognorm2 = tSNE_lognormcounts[:,2],
                tSNE_resclognorm1 = tSNE_rescaledata[:,1], tSNE_resclognorm2 = tSNE_rescaledata[:,2],
                tSNE_100pcs1 = tSNE_100pcs[:,1], tSNE_100pcs2 = tSNE_100pcs[:,2],
                UMAP_lognorm1 = umapout_lognormcounts[1,:], UMAP_lognorm2 = umapout_lognormcounts[2,:],
                UMAP_resc1 = umapout_rescaled[1,:], UMAP_resc2 = umapout_rescaled[2,:],
                UMAP_100pcs1 = umapout_100pcs[1,:], UMAP_100pcs2 = umapout_100pcs[2,:],
                scVI_z1 = latent_scvi[1,:], scVI_z2 = latent_scvi[2,:]
)
CSV.write("dimension_reduction_results_500hvgs/df_all_methods.csv", df_all, delim='\t')

dict_all = Dict("pcs" => pcs, 
                "tSNE_rescaled" => tSNE_rescaledonly,
                "tSNE_rescaledlognormcounts" => tSNE_rescaledata, 
                "tSNE_lognormcounts" => tSNE_lognormcounts,
                "tSNE_100pcs" => tSNE_100pcs,
                "umapout_lognormcounts" => umapout_lognormcounts,
                "umapout_rescaled" => umapout_rescaled,
                "UMAP_100pcs" => umapout_100pcs,
                "latent_scvi" => latent_scvi
)

using JLD2
@save "dimension_reduction_results_500hvgs/dictionary_with_all_methods.jld2" dict_all

# Plot all together
ncells, ngenes = size(adata.countmatrix)
df_stacked_all_reduced = DataFrame(z1 = vcat(pcs[:,1], tSNE_rescaledonly[:,1], umapout_lognormcounts[1,:], latent_scvi[:,1]), 
                        z2 = vcat(pcs[:,2], tSNE_rescaledonly[:,2], umapout_lognormcounts[2,:], latent_scvi[:,2]),
                        Celltype = repeat(adata.celltypes, outer = 4),
                        timepoint = repeat(adata.obs[!,"timepoints"], outer = 4),
                        Method = vcat(fill("PCA", ncells), fill("TSNE(rescaled)", ncells), fill("UMAP", ncells), fill("scVI", ncells)), 
)
CSV.write("dimension_reduction_results_500hvgs/stacked_df_for_plotting.csv", df_stacked_all_reduced, delim='\t')