Base.@kwdef mutable struct SYSTEM
    theory::String = ""
    sysname::String = ""
    nelem::Int64 = 0
    natom::Int64 = 0
    nelec::Int64 = 0
    nstate::Int64 = 0
    alat::Vector{Float64} = zeros(3)
    avec::Vector{Vector{Float64}} =[zeros(3) for i=1:3]
    bvec::Vector{Vector{Float64}} =[zeros(3) for i=1:3]
    num_kgrid::Vector{Int64} = zeros(Int64,3)
    kvec::Vector{Vector{Float64}} = [zeros(3) for i=1:3]
    num_rgrid::Vector{Int64} = zeros(Int64,3)
    rvec::Vector{Vector{Float64}} = [zeros(3) for i=1:3]
    atomic_red_coor::Vector{Vector{Float64}} = [zeros(3) for i=1:3]
end

Base.@kwdef mutable struct WFN
    nk::Vector{Int64} =zeros(Int64,3)
    nb::Int64 = 0
    nr::Vector{Int64} =zeros(Int64,3)
    orbital::Array{ComplexF64,5} =zeros(ComplexF64,1,1,1,1,1) 
end
Base.@kwdef mutable struct WFN_sort
    nk::Vector{Int64} =zeros(Int64,3)
    nb::Int64 = 0
    ng::Int64 = 0
    orbital::Array{ComplexF64,3} =zeros(ComplexF64,1,1,1) 
end

# Base.@kwdef mutable struct DM
#     nR::Vector{Int64} =zeros(Int64,3)
#     densitymatrix::Array{ComplexF64,6} = zeros(ComplexF64,1,1,1,1,1,1)
# end

Base.@kwdef mutable struct WignerFunction
    nr::Vector{Int64} =zeros(Int64,3)
    wignerfunction::Array{ComplexF64,6} = zeros(ComplexF64,1,1,1,1,1,1)
end

Base.@kwdef mutable struct Gvector
    vector::Vector{Vector{Int64}}=[zeros(Int64,3),]
    norm::Vector{Float64}=zeros(1)
end
