function sort_Gvector_wfn(salmon::SYSTEM,gvec::Gvector,wfn::WFN;cutoff::Float64=-1.0)
    nk=wfn.nk
    nk_all=prod(nk)
    nb=wfn.nb
    nr=wfn.nr
    if cutoff < 0
        println("cutoff energy :",cutoff," hartree ",27.2114*cutoff," eV)")
        cutoff = maximum(gvec.norm)
    end
    index=sortperm(gvec.norm)
    nG=length(index)    
    ng=1 # size of reduced g vector set
    for iG=1:nG
        if gvec.norm[index[iG]] > cutoff
            ng = iG-1
            break
        end
    end
    println("Reduced g vector size ng = ",ng)
    orbital_small=zeros(ComplexF64,ng,nb,nk_all)
    gvec_small=[zeros(Int64,3) for ig=1:ng]
    norm_small=zeros(ng)
    for ik=1:nk_all
        for ib=1:nb
            for iG=1:ng
                ig=index[iG]
                # orbital_small[ig,ib,ik]=wfn.orbital[]
                gvec_small[iG]=gvec.vector[ig]
                norm_small[iG]=gvec.norm[ig]
                ir=copy(gvec.vector[ig])
                for id=1:3
                    (ir[id] < 0)&&(ir[id]+=nr[id])
                    ir[id]+=1
                end
                # (ib==1 && ik==1)&&(println("ir ",ir))
                orbital_small[iG,ib,ik]=wfn.orbital[ir[1],ir[2],ir[3],ib,ik]
            end

        end
    end
    gvec_small=Gvector(gvec_small,norm_small)
    wfn_small=WFN_sort(nk,nb,ng,orbital_small)
    return gvec_small, wfn_small 
end

function real2Gspace(wfn::WFN)
    nk=wfn.nk
    nk_all=nk[1]*nk[2]*nk[3]
    nb=wfn.nb
    nr=wfn.nr
    # temp=zeros(ComplexF64,nr[1],nr[2],nr[3])
    orbital=zeros(ComplexF64,nr[1],nr[2],nr[3],nb,nk_all)
    wfg=WFN()
    wfg.nk=copy(nk)
    wfg.nb=nb
    wfg.nr=copy(nr)

    for ik=1:nk_all
        for ib=1:nb
            # for ir3=1:nr[3],ir2=1:nr[2],ir1=1:nr[1]
            #     temp[ir1,ir2,ir3]=wfn.orbital[ir1,ir2,ir3,ib,ik]
            # end
            temp=copy(wfn.orbital[:,:,:,ib,ik])
            v=fft(temp)
            orbital[:,:,:,ib,ik]=v/sqrt(prod(nr))
        end
    end
    wfg.orbital=orbital

    return wfg
end

function get_Gvector(salmon::SYSTEM)
    G=Gvector()
    nr=salmon.num_rgrid
    ng=prod(nr)
    bv=salmon.bvec
    temp1=[zeros(Int64,3) for ig=1:ng]
    temp2=zeros(ng)
    # nhalf=div.(nr,2)
    ig=1
    for ir3=1:nr[3]
        # folding ig3
        if mod(nr[3],2) == 0
            if ir3 > div(nr[3],2)+1
                ig3=ir3-nr[3]
            else
                ig3=ir3
            end
        else
            if ir3 > div(nr[3]-1,2)+1
                ig3=ir3-nr[3]
            else
                ig3=ir3
            end
        end
        for ir2=1:nr[2]
            # folding ig2
            if mod(nr[2],2) == 0
                if ir2 > div(nr[2],2)+1
                    ig2=ir2-nr[2]
                else
                    ig2=ir2
                end
            else
                if ir2 > div(nr[2]-1,2)+1
                    ig2=ir2-nr[2]
                else
                    ig2=ir2
                end
            end
            for ir1=1:nr[1]
                # folding ig1
                if mod(nr[1],2) == 0
                    if ir1 > div(nr[1],2)+1
                        ig1=ir1-nr[1]
                    else
                        ig1=ir1
                    end
                else
                    if ir1 > div(nr[1]-1,2)+1
                        ig1=ir1-nr[1]
                    else
                        ig1=ir1
                    end
                end

                # v=[ir1-1,ir2-1,ir3-1]
                v=[ig1-1,ig2-1,ig3-1]
                
                temp1[ig]=v


                u=zeros(3)
                for id=1:3
                    u+=v[id] .* bv[id]
                end
                temp2[ig]=u'*u/2.0 
                ig+=1
            end
        end
    end
    G.vector=temp1
    G.norm=temp2
    gmax=maximum(G.norm)
    println("maximum kinetic energy of G vector ",gmax," hartree (",27.2114*gmax," eV)")
    return G
end