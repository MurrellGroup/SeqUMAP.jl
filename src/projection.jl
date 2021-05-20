#define custom Distances Metric...
#see discussion of req methods here: https://github.com/JuliaStats/Distances.jl/issues/95

struct CorrectedKmer <: Metric 
    k::Int
end 

correctedkmer(a::AbstractArray, b::AbstractArray, k::Int) = CorrectedKmer(k)(a, b)
correctedkmer(a::Number, b::Number, k::Int) = CorrectedKmer(k)(a, b)

#looking back at RAD paper not sure about this... check!
function evaluate(dist::CorrectedKmer,a,b)
    D = sqeuclidean(a,b)
    #dist.N*log(1 -  minimum([D / (dist.N * 2),1])) / - dist.k
    D / (dist.k * (sum(a) + sum(b)))
end

@eval @inline (dist::CorrectedKmer)(a::AbstractArray, b::AbstractArray) = evaluate(dist, a, b)
@eval @inline (dist::CorrectedKmer)(a::Number, b::Number) = evaluate(dist, a, b)

#top level function 
function sequmap(seqs::Array{String,1}, ndim::Int;
    k = 5,
    lookup_dic = NT_DICT,
    pca = true,
    pca_maxoutdim = 5,
    n_neighbors = 12, 
    min_dist = 0.7,
    repulsion_strength = 0.1,
    metric = SqEuclidean(),
    umap_kwargs = Pair{Symbol,Any}[] #can also put all umap args together
    )
    vecs = []
    missing_chars = Set{Char}()
    for seq in seqs
        vec, missing_chars = kmer_embed(seq, k, kmer_count!; lookup_dic = lookup_dic, missing_chars = missing_chars);
        push!(vecs, vec)
    end
    if length(missing_chars) > 0
        @warn "Ignored the following characters missing from lookup: $(unique(missing_chars)). Check your sequence type!"
    end
    X = hcat(vecs...);
    X = convert(Array{Float64,2}, X)
    if pca 
        M = fit(PCA, X; maxoutdim = pca_maxoutdim)
        PCAembedding = transform(M, X);
        proj = umap(PCAembedding, ndim; n_neighbors = n_neighbors, min_dist = min_dist, metric = metric, umap_kwargs...)
    else
        proj = umap(X, ndim; n_neighbors = n_neighbors, min_dist = min_dist, metric = metric, umap_kwargs...)
    end
    return proj
end

