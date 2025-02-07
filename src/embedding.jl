#break into kmers for integers
@views collect_kmers(seq::AbstractVector{Int}, k::Int) = [seq[1+i:k+i] for i in 0:length(seq)-k]

#weights vector length is the number of unique kmers
#each kmer is counted once, contributes to one pos
function kmer_count!(kmer_counts::AbstractVector{Int}, kmer::AbstractVector{Int}, n_chars::Int)
    kmer_counts[get_kmer_index(kmer,n_chars)] += unsigned(1)
end

#this function gets the index in lexicographic order...
#remember 1-based indexing in Julia... 
function get_kmer_index(kmer::AbstractVector{Int}, n_chars::Int)
    kmer_ix = 0
    for c in kmer
        kmer_ix = kmer_ix * n_chars + (c - 1)
    end
    return kmer_ix + 1
end 

#with a `kmer_contribution()` function... 
function kmer_embed(seq, k::Int, kmer_contribution::Function; lookup_dic=AA_DICT, missing_chars=Set{Char}())
    n = length(union(values(lookup_dic)))
    kmer_counts = zeros(Int, n^k) 
    
    #break into kmers
    seqvec, missing_chars = string2encoding(seq, lookup_dic; missing_chars)
    kmers = collect_kmers(seqvec, k)
    for kmer in kmers
        if !(0 in kmer)
            kmer_contribution(kmer_counts, kmer, n)
        end
    end
    return kmer_counts, missing_chars
end
