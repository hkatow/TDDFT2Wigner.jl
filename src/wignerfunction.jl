"""
wignerfunction() computes Wigner distribution function on a given point (x,p).

For Bloch wave vector k, the momentum p must satisfy a condition that G'=p-2k-G 
is placed on the G vector grid. 

"""
function wignerfunction(salmon::SYSTEM,wfg::WFN,gvec::Gvector,x::Vector{Float64},p::Vector{Float64};verbose=0)
    nG=length(gvec.norm)
    nk_all=prod(wfg.nk)
    nr=wfg.nr
    # nb=wfg.nb
    nb=16
    gmax=div.(nr,2)
    tol=10^(-5)
    wigner=ComplexF64(0.0)
    flag=false
    for ik=1:nk_all
        k=salmon.kvec[ik]
        for ib=1:nb
            for iG=1:nG
                G1=copy(gvec.vector[iG])
                # Extract wfn amplitude C_{km,p-2k-G}
                # reciplocal coordinate
                G2=-2(p+k)-G1
                res=mod.(G2,1)
                if all(i -> i < tol, res) 
                    IG1=copy(gvec.vector[iG])
                    IG2=Int64.(round.(G2))
                    # mapping IG1 and IG2 on positive region 
                    for id=1:3
                        (IG1[id] < 0 )&&(IG1[id]+=nr[id])
                        IG1[id]+=1
                        (IG2[id] < 0 )&&(IG2[id]+=nr[id])
                        IG2[id]+=1
                    end
                    # Check if G2 is in given G space
                    condition=true
                    for id=1:3
                        condition *= (0 < IG2[id])&&(IG2[id] <=nr[id])
                    end
                    # if condition is true, there is a corrsponding G vector
                    if condition
                        flag=true
                        # Index of G1 and G2
                        # iG2=IG[1]+nr[1]*(IG[2]-1)+nr[1]*nr[2]*(IG[3]-1)
                        # phase factor exp{-i(2G1+2k-p)x}
                        v=2(G1 + k +p)
                        v_cart=zeros(3)
                        for id=1:3
                            v_cart+=v[id] * salmon.bvec[id]
                        end
                        u_cart=zeros(3)
                        for id=1:3
                            u_cart+=x[id] * salmon.avec[id]
                        end
                        phase=exp(-im*(v_cart'*u_cart))
                        # Compute Wigner Distribution function
                        C1=wfg.orbital[IG1[1],IG1[2],IG1[3],ib,ik]
                        C2=wfg.orbital[IG2[1],IG2[2],IG2[3],ib,ik]
                        wigner+=phase* conj(C1)*C2
                    end
                end
            end
        end
    end
    (verbose==1)&&(println("flag ",flag))
    return wigner
end


function get_index(vec::Vector{Float64},salmon::SYSTEM)::Vector{Int64}
    # find corresponding index to a given vector in periodic boundary condition
    index=zeros(Int64,3)
    # 
    v=zeros(3)
    for id=1:3
        v[id]=mod(vec[id],salmon.alat[id])
    end   

    return index
end

function gen_qvec(nq::Vector{Int64},salmon::SYSTEM)::Vector{Vector{Float64}}
    dq=[zeros(3) for id=1:3]
    for id=1:3
        dq[id]=salmon.avec[id] ./ nq
    end
    nq_all=nq[1]*nq[2]*nq[3]
    # println(nr_all)
    qvec=[zeros(3) for iq=1:nq_all]
    iq=1
    for iq3=0:nq[3]-1
        for iq2=0:nq[2]-1
            for iq1=0:nq[1]-1
                qvec[iq]=iq1*dq[1]+iq2*dq[2]+iq3*dq[3]
                iq+=1
            end
        end
    end
    return qvec
end
