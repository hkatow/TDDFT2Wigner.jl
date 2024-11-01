"""
wignerfunction() computes Wigner distribution function on a given point (x,p).

For Bloch wave vector k, the momentum p must satisfy a condition that G'=p-2k-G 
is placed on the G vector grid. 

"""
function wignerfunction(salmon::SYSTEM,wfg::WFN,gvec::Gvector,x::Vector{Float64},p::Vector{Float64},nb::Int64;occ=nothing, isocc::Bool=false,verbose=0)
    nG=length(gvec.norm)
    nk_all=prod(wfg.nk)
    nr=wfg.nr
    # nb=wfg.nb
    # nb=16
    gmax=div.(nr,2)
    tol=10^(-5)
    wigner=ComplexF64(0.0)
    flag=false

    # G range depends on whether nr is odd or even
    Gmin=zeros(Int64,3)
    Gmax=zeros(Int64,3)
    for id=1:3
        if mod(nr[id],2) == 0
            Gmin[id] = -div(nr[id],2)+1
            Gmax[id] =  div(nr[id],2)
        else
            Gmin[id] = -div(nr[id]-1,2)
            Gmax[id] =  div(nr[id]-1,2)
        end
    end

    for ik=1:nk_all
        k=salmon.kvec[ik]
        # check 
        if isinteger(2(p+k))
            (verbose==1)&&(println("ik=$ik 2(p+k) is integer"))
            for ib=1:nb
                for iG=1:nG
                    G1=copy(gvec.vector[iG])
                    # Extract wfn amplitude C_{km,p-2k-G}
                    # reciplocal coordinate
                    G2=-2(p+k)-G1
                    # res : remainder
                    # res=mod.(round.(G2,digits=5),1)
                    res=abs.(round.(G2)-round.(G2,digits=5))

                    if all(i -> i < tol, res) 
                        (verbose==1)&&(println("ib=$ib iG=$iG G2 is integer"))
                        IG1=copy(gvec.vector[iG])
                        IG2=Int64.(round.(G2))
                        # Check if G2 is in [ Gmin, Gmax ]
                        condition=true
                        for id=1:3
                            condition *= (Gmin[id] <= IG2[id])&&(IG2[id] <= Gmax[id])
                        end
                        # if condition is true, there is a corrsponding G vector
                        if condition
                            (verbose==1)&&(println("              G2 is in [Gmin, Gmax]"))
                            flag=true
                            # Index of G1 and G2
                            for id=1:3
                                ( IG1[id] < 0 )&&(IG1[id]+=nr[id])
                                IG1[id]+=1
                                ( IG2[id] < 0 )&&(IG2[id]+=nr[id])
                                IG2[id]+=1
                            end
                            # iG2=IG[1]+nr[1]*(IG[2]-1)+nr[1]*nr[2]*(IG[3]-1)
                            # phase factor exp{-i(2G1+2k-p)x}

                            v=2(G1 + k +p)
                            v_cart=zeros(3)
                            for id=1:3
                                v_cart+=v[id] * salmon.bvec[id]
                            end
                            u_cart=zeros(3)
                            u_cart=x
                            # for id=1:3
                            #     u_cart+=x[id] #* salmon.avec[id]
                            # end
                            phase=(1/2pi)^2*exp(-im*(v_cart'*u_cart))
                            # Compute Wigner Distribution function
                            C1=wfg.orbital[IG1[1],IG1[2],IG1[3],ib,ik]
                            C2=wfg.orbital[IG2[1],IG2[2],IG2[3],ib,ik]
                            if isocc
                                wigner+=occ[ib,ik]*phase* conj(C1)*C2
                            else
                                wigner+=phase* conj(C1)*C2
                            end
                        end
                    end
                end
            end
        end
    end
    (verbose==1)&&(println("G point match ",flag))
    return wigner
end

function isinteger(vec::Vector{Float64})::Bool
    tol=10^(-5)
    res=abs.(round.(vec)-round.(vec,digits=5))
    val=false
    if all(i -> i < tol, res) 
        val=true
    end
    return val
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
