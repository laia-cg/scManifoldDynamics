module scManifoldDynamics

using TSne
using UMAP
using VegaLite
using Statistics
using LinearAlgebra
using Random
using ProgressMeter
using CSV, DataFrames, Makie  
using scVI 
using Base: @kwdef
using Parameters: @with_kw
using Distributions
using Flux
using JLD2

include("utils.jl")
include("dimensionality_reduction.jl")
include("transformations.jl")
include("plotting.jl")
include("train_models_functions.jl")


export
    # util functions
    standardize, rescale, normalization, logNormCounts, get_data, find_clusters, create_clusters, get_celltypes, load_trained_model, save_trained_model,
    # dimensionality reduction related functions
    prcomps, get_reduced_dimensions,
    # transformation related functions
    fgeneral_transformation, params_information, fdivide_countinuous, divide_cluster, Vector_Field, transform_data!, transform_data,
    # plotting functions 
    plot_vector_field, meshgrid, plot_vector_field!, plot_tf_data, plot_tf_data_t, plot_loesses,
    # train model functions
    get_scvi_model, get_sae_model, decodefromlatent,
    # structures 
    TSNEArgs, UMAPArgs, ReductionArgs, FunctionParams
end
