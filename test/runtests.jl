using SeqUMAP
using Test

@testset "SeqUMAP.jl" begin
    seqs = [join(rand(('A','C','T','G'), 2000)) for _ in 1:100]
    seqs_clean = clean_nt.(seqs, to_strip=['-'])
    proj = sequmap(seqs_clean, 2)
    @test size(proj) == (2, 100)
end
