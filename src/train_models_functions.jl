# ===========================================
# Functions used for model training 
# ===========================================



# -------------------------------------------
#  single-cell variational autoencoder (scVI)
# -------------------------------------------


"""
    get_scvi_model(adata, seed=3841; kwargs...)

Creates and trains a standard scVI model.

# Arguments
- `adata`: An AnnData object containing the dataset.
- `seed::Int`: Random seed for reproducibility (default: `3841`).
- `kwargs...`: Additional keyword arguments for model configuration.

# Returns
- `m`: A trained scVI model.
"""
function get_scvi_model(adata, seed::Int=3841; kwargs...)
    # seed = 3841 #8271 
    Random.seed!(seed)
    library_log_means, library_log_vars = init_library_size(adata) 

    m = scVAE(size(adata.countmatrix,2);
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

    train_model!(m, adata, training_args)
    return m
end


# -------------------------------------------
#  single-cell supervised variational autoencoder (sae)
# -------------------------------------------


"""
    get_sae_model(adata, labels; seed=3841)

Creates and trains a supervised autoencoder (SAE) model.

# Arguments
- `adata`: An AnnData object containing the dataset.
- `labels::AbstractVector{<:Number}`: Labels for supervised training.
- `seed::Int`: Random seed for reproducibility (default: `3841`).

# Returns
- `m`: A trained SAE model.
"""
function get_sae_model(adata, labels::AbstractMatrix{<:Number}; seed::Int=3841)    # library_log_means, library_log_vars = init_library_size(adata, batch_key=:age)
    library_log_means, library_log_vars = init_library_size(adata)

    function encoder_loss(m, x, y)
        z, qz_m, qz_v, ql_m, ql_v, library = scVI.inference(m,x)
        enc_loss = Flux.mse(qz_m, y)
    end

    Random.seed!(seed)
    m = scVAE(size(adata.countmatrix,2);
            n_latent=2)

    training_args = TrainingArgs(
        max_epochs=200, # 50 for 10-dim 
        weight_decay=Float32(1e-6),
        register_losses = true
    )

    if training_args.register_losses
        m.loss_registry["MSE_encoder"] = []
    end

    dataloader = Flux.DataLoader((adata.countmatrix', labels'), batchsize=training_args.batchsize, shuffle=true)
    encoder_opt = Flux.Optimiser(Flux.Optimise.WeightDecay(training_args.weight_decay), ADAM(training_args.lr))
    encoder_ps = Flux.params(m.z_encoder.encoder, m.z_encoder.mean_encoder)
    progress = Progress(length(dataloader)*training_args.max_epochs);

    for epoch in 1:training_args.max_epochs
        for d in dataloader
            curloss, back = Flux.pullback(encoder_ps) do 
                encoder_loss(m, d...)    
            end
            grad = back(1f0)
            Flux.Optimise.update!(encoder_opt, encoder_ps, grad)
            if training_args.progress
                next!(progress; showvalues=[(:loss, curloss)]) # @showprogress (in front of the loop)
            end
        end
        push!(m.loss_registry["MSE_encoder"], encoder_loss(m, adata.countmatrix', labels'))
    end

    latent = get_latent_representation(m, adata.countmatrix)
    # @vlplot(:circle, x=latent[1,:], y=latent[2,:], color = adata.celltypes)

    # use this as the encoder and only train the rest 
    ncells, ngenes = size(adata.countmatrix)

    training_args = TrainingArgs(
        max_epochs=200, # 50 for 10-dim 
        weight_decay=Float32(1e-6),
        register_losses = true
        )

    if training_args.train_test_split
        trainsize = training_args.trainsize
        validationsize = 1 - trainsize # nothing 
        train_inds = shuffle!(collect(1:ncells)[1:Int(ceil(trainsize*ncells))])
    else
        train_inds = collect(1:ncells);
    end

    if training_args.register_losses
        m.loss_registry["kl_l"] = []
        m.loss_registry["kl_z"] = []
        m.loss_registry["reconstruction"] = []
        m.loss_registry["total_loss"] = []
    end
    
    ps = Flux.params(m.z_encoder.var_encoder, m.decoder) 
    opt = Flux.Optimiser(Flux.Optimise.WeightDecay(training_args.weight_decay), ADAM(training_args.lr))
    dataloader = Flux.DataLoader(adata.countmatrix', batchsize=training_args.batchsize, shuffle=true)
    progress = Progress(length(dataloader)*training_args.max_epochs);

    for epoch in 1:training_args.max_epochs
        kl_weight = scVI.get_kl_weight(training_args.n_epochs_kl_warmup, training_args.n_steps_kl_warmup, epoch, 0)
        for d in dataloader
            curloss, back = Flux.pullback(ps) do 
                scVI.loss(m, d; kl_weight=kl_weight)    
            end
            grad = back(1f0)
            Flux.Optimise.update!(opt, ps, grad)
            next!(progress; showvalues=[(:loss, curloss)])
        end
        training_args.register_losses && scVI.register_losses!(m, Float32.(adata.countmatrix[train_inds,:]'); kl_weight=kl_weight)
        push!(m.loss_registry["MSE_encoder"], encoder_loss(m, adata.countmatrix', labels'))
    end
    return m
end

# plot_total_loss = Plots.plot(m.loss_registry["total_loss"], title="Total Loss part scVI", titlefontsize=10, legend=false)
# plot_kl_z = Plots.plot(m.loss_registry["kl_z"], title="KL z", titlefontsize=10, legend=false)
# plot_reconstruction = Plots.plot(m.loss_registry["reconstruction"], title="Reconstruction", titlefontsize=10, legend=false)
# plot_kl_l = Plots.plot(m.loss_registry["kl_l"], title="KL l", titlefontsize=10, legend=false)
# plot_MSE_encoder = Plots.plot(m.loss_registry["MSE_encoder"], title="MSE Encoder", titlefontsize=10, legend=false)


# plot_loss = plot(
#     plot_total_loss,
#     plot_kl_z,
#     plot_reconstruction,
#     plot_kl_l,
#     plot_MSE_encoder,
#     layout=(2, 3),
#     size=(1000, 400)
# )


# generated_data = decodefromlatent(m, latent, adata.countmatrix)
# latent = get_latent_representation(m, adata.countmatrix)
# latent_from_generated = get_latent_representation(m, generated_data; give_mean = true)

# df_latent = DataFrame(z1 = latent_from_generated[1, :], 
#                     z2 = latent_from_generated[2, :],
#                     Celltype = cell_annotation)

# df_latent|>
#     @vlplot(width=350,
#         height=250,
#         :circle, 
#         x={:z1, title=" ", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
#         y={:z2, title=" ", axis={titleFontSize=10, labelFontSize=10, tickCount=5}}, 
#         color={
#             "Celltype:n", 
#             # scale={scheme="tableau10"}, 
#             scale={range=tableau20_extended}, 
#             legend={disable = true, title="Cell annotation", orient="right", direction="vertical"}
#         },
#         size={value=5},
#         config={legend={titleFontSize=20, labelFontSize=15}},
#         resolve={scale={x="independent", y="independent"}})

# save_trained_model("$(picked_dataset)", picked_model, m, latent_data, generated_data; original_labels = original_labels)   # uncomment this line if you want to save the trained model

"""
    decodefromlatent(m, latent, data)

Decodes reconstructed counts from latent space.

# Arguments
- `m`: A trained scVI or SAE model.
- `latent::AbstractMatrix{<:Number}`: Latent space representation of the data.
- `data::AbstractMatrix{<:Number}`: Original data to reconstruct.

# Returns
- A reconstructed dataset as a matrix of sampled counts.
"""
function decodefromlatent(m, latent::AbstractMatrix{<:Number}, data::AbstractMatrix{<:Number})
    # re-encode reconstructed counts 
    Random.seed!(32)
    z, qz_m, qz_v, ql_m, ql_v, library = scVI.inference(m,data')
    px_scale, theta, mu, zi_logits = m.decoder(latent, library)
    theta = exp.(theta)
    if m.gene_likelihood == :nb
        samp = rand.(NegativeBinomial.(theta, theta ./ (theta .+ (mu .+ eps(Float16)))))
    elseif m.gene_likelihood == :zinb
        samp = rand.(NegativeBinomial.(theta, theta ./ (theta .+ (mu .+ eps(Float16)))))
        zi_probs = scVI.logits_to_probs(zi_logits)
        is_zero = rand(Float32, size(mu)) .<= zi_probs
        samp[is_zero] .= 0.0
    else
        error("Not implemented")
    end
    return Float32.(samp')
end