using Pkg
Pkg.activate(".")

using CSV, DataFrames, VegaLite, Plots, Distributions, Random
using UMAP, TSne
using scVI  # Pkg.add(url="https://github.com/maren-ha/scVI.jl", rev="v0.1.0")
using Flux
using JLD2

using OptimalTransport
using KernelDensity
using StatsPlots


include("../src/scManifoldDynamics.jl") 
using .scManifoldDynamics  
using .scManifoldDynamics: tableau20_extended


# -------- GET ALL BASE VARIABLES READY --------------
picked_dataset = "Lymph"
picked_model = "tsne_sae"

folderdata = "./Data/Synthetic_data/$(picked_dataset)/$(picked_model)/"
folderresults = "./Data/Synthetic_data/$(picked_dataset)/$(picked_model)/Results/"

if !isdir(folderdata)
        println("Folder '$(folderdata)' does not exist.")
else
        println("Folder '$(folderdata)' exists.")
end

if !isdir(folderresults)
        println("Folder '$(folderresults)' does not exist. Creating it now.")
        mkpath(folderresults)
else
        println("Folder '$(folderresults)' already exists.")
end

cell_annotation = CSV.read("./Data/Datasets/$(picked_dataset)_cell_annotation.csv", DataFrames.DataFrame)
cell_annotation = cell_annotation[:,:x]



# -------- READ TRANSFORMED DATA -----------------

vec_transnum = ["t11", "t21"] 

for transnum in vec_transnum

        data_dict = load(folderdata*"matrices_$(transnum).jld2")["$(transnum)"]
        println(data_dict["transformation_info"])
        tf_data = data_dict["generated"]
        tf_ncells, tf_ngenes = size(tf_data)
        latent_tf_data = data_dict["latent"]

        ntimepoints = 4

        tf_adata = AnnData(countmatrix=tf_data, 
                celltypes = repeat(cell_annotation, outer=ntimepoints),
                obs = DataFrame(timepoint =  cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1))
        )

        translognormcounts = log.(normalization(tf_data) .+ 1)  #logNormCounts(tf_data, 2)

        # # --------- Look at the loaded transformation (if needed) ---------------
        # experiment = "UMAP_parameters"
        # check = "latent_data"
        # folderexperiments = "C:/Users/canall/new-results-time-transformation-project/$(experiment)/"
        # check_data = latent_tf_data  #umapout' #umapout' # latent_tf_data #umapout_pcs_10'


        df_latent_tf_data = DataFrames.DataFrame(z1 = latent_tf_data[:,1], 
                                                z2 = latent_tf_data[:,2], 
                                                timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
                                                Celltype  = repeat(cell_annotation, outer=ntimepoints), 
        )

        plot_original_latent = false
        if plot_original_latent
                plot_transformed_data = df_latent_tf_data |>  
                @vlplot(width=500,
                        height=400,
                        :circle, 
                        x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
                        y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
                        column = {:"timepoint:n",  header={title = "Timepoint", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
                        # color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
                        color={"Celltype:n", scale={range=tableau20_extended}, legend={disable=true, title="Clusters",orient="right"}},
                        #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
                        size={value=15},
                        labelFontSize=50,
                        config={legend={titleFontSize=20, labelFontSize=20}}
                )
                # VegaLite.save(folderexperiments*"$(dataset)_TIME_$(check).png", plot_transformed_data)
        end

        # # ---------- Explained variance -----------
        # experiment = "explained_variance"
        # folderexperiments = "C:/Users/canall/new-results-time-transformation-project/$(experiment)/"
        # pcs, explained_variance_ratio = prcomps2(translognormcounts) #translognormcounts
        # zoom = true
        # pcs_var_plot = plot_explained_variance(explained_variance_ratio, zoom; n_pcs=50)
        # Plots.savefig(pcs_var_plot, folderexperiments * dataset * "timedataset_variance_explained_zoom_$(zoom)_translognormcounts.png")
        # # pcs_2, explained_variance_ratio_2 = prcomps(tf_data) #translognormcounts


        # ---------- Parameters for the models -----------


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


        # ---------- RUN MODELS ---------------------

        # to make some space in the memory
        pcs = nothing
        tSNE = nothing
        umapout = nothing
        df_reduced_dimensions_tf_data = nothing
        tf_latent_data = nothing


        translognormcounts = log.(normalization(tf_data) .+ 1)  #log
        # pcs = scVI.prcomps(log.(normalization(tf_adata.countmatrix) .+ 1))
        elapsed_time_seconds = @elapsed begin
                pcs = scVI.prcomps(translognormcounts)
        end

        # firstpcs = pcs[:,1:100]

        elapsed_time_seconds = @elapsed begin
                Random.seed!(big_TSNE_args.seed) 
                tSNE = tsne(rescale(tf_data, dims=1), 2, big_TSNE_args.pca_dims, 1000, big_TSNE_args.perplexity, progress=true); 
        end
        minutes, seconds = divrem(elapsed_time_seconds, 60)
        println("Execution time tSNE with $(big_TSNE_args.pca_dims) pcs: ", minutes, " minutes, ", seconds, " seconds")



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

        # umapout = umap(Transpose(log.(normalization(tf_adata.countmatrix) .+ 1)))
        # umapout = umap(Transpose(translognormcounts))

        # tf_m_scvi = get_scvi_model(tf_adata, 3841)


        Random.seed!(3841)
        library_log_means, library_log_vars = init_library_size(tf_adata) 

        tf_m_scvi = scVAE(size(tf_adata.countmatrix,2);
                n_layers=1,
                n_latent=2,
                library_log_means=library_log_means,
                library_log_vars=library_log_vars
                )


        training_args = TrainingArgs(
                train_test_split=false, 
                lr = 1e-3, #1e-3 for 10 dim
                batchsize=128, 
                max_epochs=200, # 50 for 10-dim 
                weight_decay=Float32(1e-6),
                register_losses = true,
                verbose=false
                )

        train_model!(tf_m_scvi, tf_adata, training_args)

        plot_loesses(tf_m_scvi.loss_registry)
        # loss_registry = tf_m_scvi.loss_registry

        # # function plot_loesses(loss_registry::Dict)
        #         epochs = length(loss_registry["total_loss"])
        #         df_loss = DataFrame(
        #             epoch = repeat(1:epochs, 3),
        #             loss_value = vcat(loss_registry["total_loss"], loss_registry["kl_z"], loss_registry["reconstruction"]),
        #             loss_type = repeat(["Total Loss", "KL Divergence", "Reconstruction"], inner=epochs)
        #         )
        
        #         plot_loss = df_loss |>
        #         @vlplot(
        #             :line,
        #             width=400,
        #             height=300,
        #             x={field=:epoch, scale = {zero = false}, type="quantitative", title="Epoch", axis={titleFontSize=12, labelFontSize=12}}, # 
        #             y={field=:loss_value, scale = {zero = false}, type="quantitative", title="Loss Value", axis={titleFontSize=12, labelFontSize=12}}, #
        #             column = {"loss_type:n", header={title = nothing, labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}},
        #             resolve={scale={x="independent",y="independent"}}
        #         )
        #         return plot_loss
        # #     end

        tf_latent_data = get_latent_representation(tf_m_scvi, tf_adata.countmatrix)


        df_reduced_dimensions_tf_data = DataFrame(z1 = vcat(latent_tf_data[:,1], pcs[:,1], tSNE[:,1], umapout[1,:], tf_latent_data[1,:]), 
                                z2 = vcat(latent_tf_data[:,2], pcs[:,2], tSNE[:,2], umapout[2,:], tf_latent_data[2,:]),
                                Celltype = repeat(get_celltypes(tf_adata), outer = 5),
                                timepoint = repeat(tf_adata.obs.timepoint,  outer = 5),
                                Method = vcat(fill("Orig. Transf. $(picked_model)", tf_ncells), fill("PCS", tf_ncells), fill("tSNE", tf_ncells), fill("UMAP", tf_ncells), fill("scVI", tf_ncells)), 
        )


        plot_reduced_dimensions_tf_data = df_reduced_dimensions_tf_data |>
        @vlplot(width=400,
                height=300,
                :circle, 
                x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
                y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
                column = {"timepoint:n", header={title = "Timepoint", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}},
                row = {"Method:n", header={title = nothing, labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, #sort=["Orig. Transf. $(picked_model)", "PCS, tSNE, UMAP, scVI"]
                color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
                color={"Celltype:n", scale={range=tableau20_extended}, legend={disable=true, title="Clusters",orient="right"}},
                #color={:timepoint, legend={disable=true, title="Clusters",orient="right"}},
                size={value=15},
                config={legend={titleFontSize=20, labelFontSize=15}},
                resolve={scale={x="independent",y="independent"}} 
        )


        # VegaLite.save(folderresults*"$(transnum).pdf", plot_reduced_dimensions_tf_data)
        VegaLite.save(folderresults*"$(transnum).png", plot_reduced_dimensions_tf_data)

        CSV.write(folderresults*"$(transnum)_df_reduced_dimensions_tf_data.csv", df_reduced_dimensions_tf_data)
end


# -------- GET ALL BASE VARIABLES READY --------------
picked_dataset = "Lymph"
picked_model = "umap_sae"

folderdata = "./Data/Synthetic_data/$(picked_dataset)/$(picked_model)/"
folderresults = "./Data/Synthetic_data/$(picked_dataset)/$(picked_model)/Results/"

if !isdir(folderdata)
        println("Folder '$(folderdata)' does not exist.")
else
        println("Folder '$(folderdata)' exists.")
end

if !isdir(folderresults)
        println("Folder '$(folderresults)' does not exist. Creating it now.")
        mkpath(folderresults)
else
        println("Folder '$(folderresults)' already exists.")
end

cell_annotation = CSV.read("./Data/Datasets/$(picked_dataset)_cell_annotation.csv", DataFrames.DataFrame)
cell_annotation = cell_annotation[:,:x]



# -------- READ TRANSFORMED DATA -----------------

vec_transnum = ["t21"] 

# transnum = "t21"

for transnum in vec_transnum

        data_dict = load(folderdata*"matrices_$(transnum).jld2")["$(transnum)"]
        println(data_dict["transformation_info"])
        tf_data = data_dict["generated"]
        tf_ncells, tf_ngenes = size(tf_data)
        latent_tf_data = data_dict["latent"]

        ntimepoints = 4

        tf_adata = AnnData(countmatrix=tf_data, 
                celltypes = repeat(cell_annotation, outer=ntimepoints),
                obs = DataFrame(timepoint =  cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1))
        )

        translognormcounts = log.(normalization(tf_data) .+ 1)  #logNormCounts(tf_data, 2)

        # # --------- Look at the loaded transformation (if needed) ---------------
        # experiment = "UMAP_parameters"
        # check = "latent_data"
        # folderexperiments = "C:/Users/canall/new-results-time-transformation-project/$(experiment)/"
        # check_data = latent_tf_data  #umapout' #umapout' # latent_tf_data #umapout_pcs_10'


        df_latent_tf_data = DataFrames.DataFrame(z1 = latent_tf_data[:,1], 
                                                z2 = latent_tf_data[:,2], 
                                                timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1),
                                                Celltype  = repeat(cell_annotation, outer=ntimepoints), 
        )

        plot_original_latent = false
        if plot_original_latent
                plot_transformed_data = df_latent_tf_data |>  
                @vlplot(width=500,
                        height=400,
                        :circle, 
                        x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
                        y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
                        column = {:"timepoint:n",  header={title = "Timepoint", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
                        # color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
                        color={"Celltype:n", scale={range=tableau20_extended}, legend={disable=true, title="Clusters",orient="right"}},
                        #color={:timepoint, legend={disable=false, title="Time",orient="right"}},
                        size={value=15},
                        labelFontSize=50,
                        config={legend={titleFontSize=20, labelFontSize=20}}
                )
                # VegaLite.save(folderexperiments*"$(dataset)_TIME_$(check).png", plot_transformed_data)
        end

        # # ---------- Explained variance -----------
        # experiment = "explained_variance"
        # folderexperiments = "C:/Users/canall/new-results-time-transformation-project/$(experiment)/"
        # pcs, explained_variance_ratio = prcomps2(translognormcounts) #translognormcounts
        # zoom = true
        # pcs_var_plot = plot_explained_variance(explained_variance_ratio, zoom; n_pcs=50)
        # Plots.savefig(pcs_var_plot, folderexperiments * dataset * "timedataset_variance_explained_zoom_$(zoom)_translognormcounts.png")
        # # pcs_2, explained_variance_ratio_2 = prcomps(tf_data) #translognormcounts


        # ---------- Parameters for the models -----------


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


        # ---------- RUN MODELS ---------------------

        # to make some space in the memory
        pcs = nothing
        tSNE = nothing
        umapout = nothing
        df_reduced_dimensions_tf_data = nothing
        tf_latent_data = nothing


        translognormcounts = log.(normalization(tf_data) .+ 1)  #log
        pcs, explained_variance2 = prcomps2(translognormcounts)
        # pcs = scVI.prcomps(log.(normalization(tf_adata.countmatrix) .+ 1))
        elapsed_time_seconds = @elapsed begin
                pcs = scVI.prcomps(translognormcounts)
        end

        # firstpcs = pcs[:,1:100]

        elapsed_time_seconds = @elapsed begin
                Random.seed!(big_TSNE_args.seed) 
                tSNE = tsne(rescale(tf_data, dims=1), 2, big_TSNE_args.pca_dims, 1000, big_TSNE_args.perplexity, progress=true); 
        end
        minutes, seconds = divrem(elapsed_time_seconds, 60)
        println("Execution time tSNE with $(big_TSNE_args.pca_dims) pcs: ", minutes, " minutes, ", seconds, " seconds")



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

        # umapout = umap(Transpose(log.(normalization(tf_adata.countmatrix) .+ 1)))
        # umapout = umap(Transpose(translognormcounts))

        # tf_m_scvi = get_scvi_model(tf_adata, 3841)


        Random.seed!(3841)
        library_log_means, library_log_vars = init_library_size(tf_adata) 

        tf_m_scvi = scVAE(size(tf_adata.countmatrix,2);
                n_layers=1,
                n_latent=2,
                library_log_means=library_log_means,
                library_log_vars=library_log_vars
                )


        training_args = TrainingArgs(
                train_test_split=false, 
                lr = 1e-3, #1e-3 for 10 dim
                batchsize=128, 
                max_epochs=200, # 50 for 10-dim 
                weight_decay=Float32(1e-6),
                register_losses = true,
                verbose=false
                )

        train_model!(tf_m_scvi, tf_adata, training_args)

        plot_loesses(tf_m_scvi.loss_registry)
        # loss_registry = tf_m_scvi.loss_registry

        # # function plot_loesses(loss_registry::Dict)
        #         epochs = length(loss_registry["total_loss"])
        #         df_loss = DataFrame(
        #             epoch = repeat(1:epochs, 3),
        #             loss_value = vcat(loss_registry["total_loss"], loss_registry["kl_z"], loss_registry["reconstruction"]),
        #             loss_type = repeat(["Total Loss", "KL Divergence", "Reconstruction"], inner=epochs)
        #         )
        
        #         plot_loss = df_loss |>
        #         @vlplot(
        #             :line,
        #             width=400,
        #             height=300,
        #             x={field=:epoch, scale = {zero = false}, type="quantitative", title="Epoch", axis={titleFontSize=12, labelFontSize=12}}, # 
        #             y={field=:loss_value, scale = {zero = false}, type="quantitative", title="Loss Value", axis={titleFontSize=12, labelFontSize=12}}, #
        #             column = {"loss_type:n", header={title = nothing, labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}},
        #             resolve={scale={x="independent",y="independent"}}
        #         )
        #         return plot_loss
        # #     end

        tf_latent_data = get_latent_representation(tf_m_scvi, tf_adata.countmatrix)


        df_reduced_dimensions_tf_data = DataFrame(z1 = vcat(latent_tf_data[:,1], pcs[:,1], tSNE[:,1], umapout[1,:], tf_latent_data[1,:]), 
                                z2 = vcat(latent_tf_data[:,2], pcs[:,2], tSNE[:,2], umapout[2,:], tf_latent_data[2,:]),
                                Celltype = repeat(get_celltypes(tf_adata), outer = 5),
                                timepoint = repeat(tf_adata.obs.timepoint,  outer = 5),
                                Method = vcat(fill("Orig. Transf. $(picked_model)", tf_ncells), fill("PCS", tf_ncells), fill("tSNE", tf_ncells), fill("UMAP", tf_ncells), fill("scVI", tf_ncells)), 
        )


        plot_reduced_dimensions_tf_data = df_reduced_dimensions_tf_data |>
        @vlplot(width=400,
                height=300,
                :circle, 
                x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
                y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
                column = {"timepoint:n", header={title = "Timepoint", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}},
                row = {"Method:n", header={title = nothing, labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, #sort=["Orig. Transf. $(picked_model)", "PCS, tSNE, UMAP, scVI"]
                # color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
                color={"Celltype:n", scale={range=tableau20_extended}, legend={disable=true, title="Clusters",orient="right"}},
                #color={:timepoint, legend={disable=true, title="Clusters",orient="right"}},
                size={value=15},
                config={legend={titleFontSize=20, labelFontSize=15}},
                resolve={scale={x="independent",y="independent"}} 
        )


        # VegaLite.save(folderresults*"$(transnum).pdf", plot_reduced_dimensions_tf_data)
        VegaLite.save(folderresults*"$(transnum).png", plot_reduced_dimensions_tf_data)

        CSV.write(folderresults*"$(transnum)_df_reduced_dimensions_tf_data.csv", df_reduced_dimensions_tf_data)
end




# ---------- Quantitative measure study ------------------

latent_tf_data = min_max_scale(latent_tf_data) # It re-arranges the data in the 1-square for computational stability


df_t1 = df_latent_tf_data[df_latent_tf_data.timepoint .== 1, :]  # all cells from timepoint 1
df_t2 = df_latent_tf_data[df_latent_tf_data.timepoint .== 2, :]
df_t3 = df_latent_tf_data[df_latent_tf_data.timepoint .== 3, :]
df_t4 = df_latent_tf_data[df_latent_tf_data.timepoint .== 4, :]
        

# ------ KERNEL ESTIMATION ------
timepoint1_data = [df_t1.z1 df_t1.z2]
timepoint2_data = [df_t2.z1 df_t2.z2]
timepoint3_data = [df_t3.z1 df_t3.z2]
timepoint4_data = [df_t4.z1 df_t4.z2]

kde_npoints = 200 #64
kde_t1 = kde(timepoint1_data; npoints=(kde_npoints,kde_npoints))  # or kde((timepoint1_data[:,1], timepoint1_data[:,2]))
kde_t2 = kde(timepoint2_data; npoints=(kde_npoints,kde_npoints))  
kde_t3 = kde(timepoint3_data; npoints=(kde_npoints,kde_npoints))  
kde_t4 = kde(timepoint4_data; npoints=(kde_npoints,kde_npoints))  

plot_kde_t1 = plot(kde_t1,  linewidth=2, legend = false);
plot_kde_t2 = plot(kde_t2,  linewidth=2, legend = false);
plot_kde_t3 = plot(kde_t3, linewidth=2, legend = false);
plot_kde_t4 = plot(kde_t4, linewidth=2, legend = false);

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




# -------- Try with different number of points at each timepoint  -----------------


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

timepoint = cat(collect(fill(i, length(cell_annotation)) for i in 1:ntimepoints)..., dims=1)
ncells = length(cell_annotation)
ncells_per_timepoint = [Int(0.8*ncells), Int(0.4*ncells), Int(0.6*ncells), Int(0.3*ncells)]
selected_indices_dtp = select_cells_per_timepoint(timepoint, ncells_per_timepoint)

tf_data_dtp = tf_data[selected_indices_dtp,:]

translognormcounts = log.(normalization(tf_data) .+ 1)  
translognormcounts_dtp = log.(normalization(tf_data_dtp) .+ 1) 
timepoint_dtp = timepoint[selected_indices_dtp]

# pcs = scVI.prcomps(log.(normalization(tf_adata.countmatrix) .+ 1))
elapsed_time_seconds = @elapsed begin
        pcs_dtp = scVI.prcomps(translognormcounts_dtp)
end

elapsed_time_seconds = @elapsed begin
        Random.seed!(big_UMAP_args.seed)
        # umapout_pcs_10 = umap(pcs[:,1:big_UMAP_args.pca_dims]', min_dist=0.5) 
        umapout_dtp = umap(Transpose(translognormcounts_dtp), min_dist=0.5)
        # umapout = umap(Transpose(translognormcounts), min_dist=0.5)
        # umapout_origdata = umap(tf_data')
end
minutes, seconds = divrem(elapsed_time_seconds, 60)
println("Execution time: ", minutes, " minutes, ", seconds, " seconds")


df_scatter_transformed = DataFrames.DataFrame(z1 = vcat(umapout[1,:], umapout_dtp[1,:]), 
        z2 = vcat(umapout[2,:], umapout_dtp[2,:]), 
        timepoint = vcat(timepoint, timepoint_dtp),
        Celltype  = vcat(repeat(cell_annotation, outer=ntimepoints), repeat(cell_annotation, outer=ntimepoints)[selected_indices_dtp]), 
        Method = vcat(fill("All points UMAP", length(timepoint)), fill("Rnadomly selected points UMAP", length(timepoint_dtp)))
)


plot_transformed_data = df_scatter_transformed |>  
        @vlplot(width=500,
        height=400,
        :circle, 
        x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}}, 
        # column = {:"timepoint:n",  header={title = "Timepoint", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}}, 
        row = {:"Method:n"},
        # color={"Celltype:n", scale={scheme="tableau10"}, legend={disable=true, title="Clusters",orient="right"}},
        color={:timepoint, legend={disable=false, title="Time",orient="right"}},
        size={value=30},
        labelFontSize=50,
        config={legend={titleFontSize=20, labelFontSize=20}}
)

# =========================================
#       Modify style of the plots
# =========================================

df_reduced_dimensions_tf_data = CSV.read(folderresults * "$(transnum)_df_reduced_dimensions_tf_data.csv", DataFrame)


plot_reduced_dimensions_tf_data = df_reduced_dimensions_tf_data |>
@vlplot(
    width=600,
    height=450,
    :circle,
    x={:z1, title="z1", axis={titleFontSize=15, labelFontSize=15, tickCount=5}},
    y={:z2, title="z2", axis={titleFontSize=15, labelFontSize=15, tickCount=5}},
    column={"timepoint:n", header={title="Timepoint", labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}},
    row={"Method:n", header={title=nothing, labelFontSize=15, titleFontSize=15, labelFontWeight="bold"}},
    color={"Celltype:n", scale={range=tableau20_extended}, legend={disable=true, title="Clusters", orient="right"}},
    size={value=10},  # Smaller point size here!
    config={legend={titleFontSize=20, labelFontSize=15}},
    resolve={scale={x="independent", y="independent"}}
)

VegaLite.save(folderresults*"$(transnum).png", plot_reduced_dimensions_tf_data)
VegaLite.save(folderresults*"$(transnum).pdf", plot_reduced_dimensions_tf_data)


CSV.write(folderresults * "$(transnum)_df_reduced_dimensions_tf_data.csv", df_reduced_dimensions_tf_data)
