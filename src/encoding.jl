#Amino acid...
const AA_DICT = Dict(
    'A' => 1, 'C' => 2, 'D' => 3, 'E' => 4, 'F' => 5, 
    'G' => 6, 'H' => 7, 'I' => 8, 'K' => 9, 'L' => 10, 
    'M' => 11, 'N' => 12, 'P' => 13, 'Q' => 14, 'R' => 15, 
    'S' => 16, 'T' => 17, 'V' => 18, 'W' => 19, 'Y' => 20, 
    '*' => 21, 'X' => 22
    )

#possibly add ambig handling for aas...

#Nucleotide...
const NT_DICT = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)

#disambiguate! 
const IUPACbool = Dict{Char,Array{Bool,1}}(Dict())
IUPACbool['A']=[true,false,false,false]
IUPACbool['C']=[false,true,false,false]
IUPACbool['G']=[false,false,true,false]
IUPACbool['T']=[false,false,false,true]
IUPACbool['U']=[false,false,false,true]
IUPACbool['R']=[true,false,true,false]
IUPACbool['Y']=[false,true,false,true]
IUPACbool['S']=[false,true,true,false]
IUPACbool['W']=[true,false,false,true]
IUPACbool['K']=[false,false,true,true]
IUPACbool['M']=[true,true,false,false]
IUPACbool['B']=[false,true,true,true]
IUPACbool['D']=[true,false,true,true]
IUPACbool['H']=[true,true,false,true]
IUPACbool['V']=[true,true,true,false]
IUPACbool['N']=[true,true,true,true];

function resolve_base(c)
    uc = uppercase(c)
    if uc in keys(IUPACbool)
        return ['A','C','G','T'][rand((1:4)[IUPACbool[uc]])]
    else
        return uc
    end
end

function resolve_seq(s)
    return join(resolve_base.(collect(s)))
end

#encode string
function string2encoding(seq::String, lookup_dic::Dict{Char,Int}; missing_chars = Set{Char}())
    vec = Int64[]
    for c in seq
        if !haskey(lookup_dic, c)
            push!(missing_chars, c)
            push!(vec, 0)
        else
            push!(vec, lookup_dic[c])
        end
    end
    return vec, missing_chars
end