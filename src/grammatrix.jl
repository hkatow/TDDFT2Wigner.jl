function grammatrix(wfn::WFN,salmon::SYSTEM,index1::Vector{Int64},index2::Vector{Int64})::ComplexF64
    ib1=index1[1]; ik1=index1[2]
    ib2=index2[1]; ik2=index2[2]
    nr=salmon.num_rgrid
    # alat=salmon.alat
    # dr = alat ./ nr
    # dv=dr[1]*dr[2]*dr[3]
    v=ComplexF64(0.0)    
    for ir3=1:nr[3],ir2=1:nr[2],ir1=1:nr[1]
        v+=conj(wfn.orbital[ir1,ir2,ir3,ib1,ik1])*wfn.orbital[ir1,ir2,ir3,ib2,ik2]
    end
    # v = v * dv
    return v
end