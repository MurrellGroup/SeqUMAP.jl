# SeqUMAP.jl

Alignment-free embedding of sequences via kmer-counting and UMAP.jl.

### Setup

From the Julia REPL:

```julia-repl
using Pkg
Pkg.add(PackageSpec(;name="SeqUMAP",url="https://github.com/MurrellGroup/SeqUMAP.jl.git"))
```

### Usage

Import sequences from a file by running 

```julia-repl
seqnames, seqs = read_fasta("my_sequences.fa");
```

Optionally remove gaps and convert 'U' => 'T' with

```julia-repl
seqs_clean = clean_nt.(seqs; to_strip = ['-']);
```

Then, obtain a SeqUMAP projection, using the default settings, with 

```
proj = sequmap(seqs_clean, 2)
```

The full list of options are:
```
proj = sequmap(seqs_clean, ndim; 
    k = 5, lookup_dic = NT_DICT, pca = true, pca_maxoutdim = 5, 
    n_neighbors = 12, min_dist = 0.7, repulsion_strength = 0.1, 
    metric = SqEuclidean(), umap_kwargs = Pair{Symbol,Any}[]
    )
```

where `seqs` is an array of strings, `ndim` is the number of projected dimensions and `k` is the kmer size. If using nucleotide sequences with an {'A', 'T', 'C', 'G'} alphabet, pass `SeqUMAP.NT_DICT` as the character lookup. If `pca` is specified, PCA will be run on kmer vectors prior to UMAP embedding with `pca_maxoutdim` components. 

UMAP parameters `n_neighbors`, `min_dist`, `repulsion_strength` and `metric` can also be specified, with the default options shown here. If you want to pass UMAP additional keyword arguments, you can do so by providing them in `umap_kwargs`.

Both AA and NT dictionaries are included for sequence encoding. For AA, k=2 is usually sufficient due to the larger alphabet. For NT, start with k=5 or k=6. 

```
AA_DICT = Dict(
    'A' => 1, 'C' => 2, 'D' => 3, 'E' => 4, 'F' => 5, 
    'G' => 6, 'H' => 7, 'I' => 8, 'K' => 9, 'L' => 10, 
    'M' => 11, 'N' => 12, 'P' => 13, 'Q' => 14, 'R' => 15, 
    'S' => 16, 'T' => 17, 'V' => 18, 'W' => 19, 'Y' => 20, 
    '*' => 21, 'X' => 22
    )

NT_DICT = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)
```

Embeddings can be easily visualized with PyPlot or other common Julia plotting libraries. 

```julia-repl
using PyPlot

fig, ax = subplots()
ax.scatter(proj[1, :], proj[2, :]; s = 2.0, linewidth = 0.0)
```
