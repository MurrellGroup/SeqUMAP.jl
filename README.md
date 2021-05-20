# SeqUMAP.jl

[In progress.] Encode and embed sequences as kmer vectors (both NT and AA) and generate projections using UMAP.jl.

### Setup

```
using Pkg
Pkg.add(PackageSpec(;name="SeqUMAP",url="https://github.com/murrellb/SeqUMAP.jl.git"))
```

### Usage

Embed and project sequences using

```
proj = sequmap(seqs, ndim; k = 2, lookup_dic = NT_DICT, n_neighbors = 30, min_dist = 1e-3)
```

where `seqs` is an array of strings, `ndim` is the number of projected dimensions and `k` is the kmer size. By default, sequences are embedded by counting all unique kmers and the distance between seuqences is estimated using corrected kmer distance, provided as the `CorrectedKmer(k)` metric. The following lookup dictionaries are included for sequence encoding:

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