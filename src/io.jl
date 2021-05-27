"""
    read_fasta(filename; seqtype=String)

Read a fasta file, returning sequences as the specified type. 

# Examples
```julia-repl
seqnames, seqs = read_fasta(filename)
````
"""
function read_fasta(filename; seqtype=String)
    stream = open(FASTA.Reader, filename)
    seqnames, seqs = String[], seqtype[]
    for entry in stream
        push!(seqs, FASTA.sequence(seqtype, entry))
        push!(seqnames, FASTA.identifier(entry))
    end
    return seqnames, seqs
end

"""
    read_fastq_records(filename; seqtype=String)
Read .fastq file contents.

# Examples
```julia-repl
seqs, phreds, seqnames = read_fastq(filename)
```
"""
function read_fastq(filename; seqtype=String)
    stream = open(FASTQ.Reader, filename)
    seqs, phreds, seqnames = seqtype[], Array{Int8, 1}[], String[]
    for entry in stream
        push!(seqs, FASTQ.sequence(seqtype, entry))
        push!(phreds, FASTQ.quality(entry, :sanger))
        push!(seqnames, FASTQ.identifier(entry))
    end
    return seqs, phreds, seqnames
end

"""
    clean_nt(seq::String; to_strip = ['.', '-'])
Enforce case, remove gaps and convert RNA to DNA in a single step. Alternatively use Biosequences types for more control.

# Examples
```julia-repl
degapped_dna = clean_nt("AUG---GUGUAG")
"ATGGTGTAG"
```
"""
function clean_nt(seq::String; to_strip = ['.', '-'])
    dna_seq = replace(uppercase(seq), 'U' => 'T')
    for c in to_strip
        dna_seq = replace(dna_seq, c => "")
    end
    return dna_seq
end

"""
    clean_nt(seq::String; to_strip = ['.', '-'])
Enforce case and remove gaps in a single step. Alternatively use Biosequences types for more control.

# Examples
```julia-repl
degapped_aa = clean_aa("mvk--lvss*")
"MVKLVSS*"
```
"""
function clean_aa(seq::String; to_strip = ['.', '-'])
    aa_seq = uppercase(seq)
    for c in to_strip
        aa_seq = replace(aa_seq, c => "")
    end
    return aa_seq
end