module SeqUMAP

using Distances
using FASTX
import MultivariateStats: fit, PCA, transform
import UMAP: umap


include("encoding.jl")
include("embedding.jl")
include("projection.jl")
include("io.jl")

#encoding.jl...
export AA_DICT,
NT_DICT,
IUPACbool,
resolve_base,
resolve_seq,
string2encoding,

#embedding.jl...
collect_kmers,
kmer_count!,
get_kmer_index,
kmer_embed,

#projection.jl...
CorrectedKmer,
sequmap,
seqpca,

#io.jl
read_fasta,
read_fastq,
clean_nt, 
clean_aa

end