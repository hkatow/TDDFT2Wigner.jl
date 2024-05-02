"""
Read wave function generated by real time evolution simulation of SALMON.
"""
function read_wfn_rt(restart_dir::String,salmon::SYSTEM,rank::Vector{Int64})::WFN
    # note: file index begins from 0
    nfile=rank[1]+1
    
    # system parameters
    nk=salmon.num_kgrid
    nr=salmon.num_rgrid
    nb=salmon.nstate
    nk_all=prod(nk)
    nr_all=prod(nr)
    
    # volume element
    a=salmon.avec
    V=abs(a[1]'*cross(a[2],a[3]))
    dV=V/nr_all

    wfn=WFN()
    wfn.nk=nk
    wfn.nr=nr
    wfn.nb=nb
    wfn.orbital=zeros(ComplexF64,nr[1],nr[2],nr[3],nb,nk_all)
    
    # Suppose that nk states are stored in each file
    for irank=rank[1]:rank[2]
        if irank < 10
            file="$restart_dir/rank_00000$irank/wfn.bin"
        else
            file="$restart_dir/rank_0000$irank/wfn.bin"
        end 
        println("Read orbital function from $file")
        ib=irank+1
        io=open(file,"r")
        for ik=1:nk_all
            # println("ik",ik)
            for is=1:1
                v=zeros(ComplexF64,nr_all)            
                read!(io,v)
                ir=1
                for ir3=1:nr[3], ir2=1:nr[2], ir1=1:nr[1] 
                    wfn.orbital[ir1,ir2,ir3,ib,ik]=v[ir]*sqrt(dV)
                    ir+=1
                end
            end
        end
        close(io)
    end

    return wfn
end
