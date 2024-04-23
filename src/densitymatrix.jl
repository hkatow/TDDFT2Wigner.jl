function densitymatrix(wfn::WFN,salmon::SYSTEM)::DM
    t_init=now()
    println("Construct density matrix on supercell. ",t_init)

    nr=salmon.num_rgrid
    nk=salmon.num_kgrid
    kv=salmon.kvec
    nr_all=nr[1]*nr[2]*nr[3]
    nk_all=nk[1]*nk[2]*nk[3]
    nrho=nk .* nr
    nrho_all=nrho[1]*nrho[2]*nrho[3]
    println("Matrix size : ",nrho[1],"*",nrho[2],"*",nrho[3],"=",nrho_all)

    # dv : volume element
    a=salmon.avec
    volume=abs(a[1]'*cross(a[2],a[3]))
    dv=volume/nr_all

    # k point weight
    k_weight=1/nk_all

    nb=salmon.nstate
    rho=DM()
    # density matrix is constructed on Nk*V supercell
    rho_temp=zeros(ComplexF64,nrho[1],nrho[2],nrho[3],nrho[1],nrho[2],nrho[3])
    # generate R grid on super cell
    Rvec=get_Rvec(salmon)
    
    # t=now()
    # compute rho(x,y)
    iY=1    
    @inbounds for iY3=1:nrho[3],iY2=1:nrho[2],iY1=1:nrho[1]
        # index : supercell -> unit cell
        iy3=mod(iY3,nk[3])+1
        iy2=mod(iY2,nk[2])+1
        iy1=mod(iY1,nk[1])+1

        t_now=now()
        h=hour(t_now)
        m=minute(t_now)
        s=second(t_now)

        (mod(iY,100) == 0)&&(println("$iY/$(nrho_all) $h:$m:$s"))
        Y=Rvec[iY]
        iX=1
        for iX3=1:nrho[3],iX2=1:nrho[2],iX1=1:nrho[1]
            # index : supercell -> unit cell
            ix3=mod(iX3,nk[3])+1
            ix2=mod(iX2,nk[2])+1
            ix1=mod(iX1,nk[1])+1

            # x=salmon.rvec[ix]
            X=Rvec[iX]
            v=ComplexF64(0.0)
            for ik=1:nk_all
                for ib=1:nb
                    # println("ik",ik)
                    ekx=exp(im*(kv[ik]'*X))
                    eky=exp(im*(kv[ik]'*Y))
                    c1=ekx*wfn.orbital[ix1,ix2,ix3,ib,ik]
                    c2=eky*wfn.orbital[iy1,iy2,iy3,ib,ik]
                    v+=conj(c1)*c2
                end
            end
            rho_temp[iX1,iX2,iX3,iY1,iY2,iY3]=v*k_weight*dv
            iX+=1
        end    
        iY+=1
    end
    t_fin=now()
    dt=t_fin-t_init
    println("density matrix calculation finished within $(dt.value/60000) minutes")
    rho.nR=nrho
    rho.densitymatrix=rho_temp    
    return rho 
end

function get_Rvec(salmon::SYSTEM)
    # dr::r point interval
    dr=[zeros(3) for id=1:3]
    for id=1:3
        dr[id]=salmon.avec[id] ./ salmon.num_rgrid[id]
    end
    # println(dr)
    nr=salmon.num_rgrid
    nk=salmon.num_kgrid
    nR=nk .* nr
    nR_all=nR[1]*nR[2]*nR[3]

    # println(nr_all)
    Rvec=[zeros(3) for iR=1:nR_all]
    iR=1
    for iR3=0:nR[3]-1
        for iR2=0:nR[2]-1
            for iR1=0:nR[1]-1
                Rvec[iR]=iR1*dr[1]+iR2*dr[2]+iR3*dr[3]
                iR+=1
            end
        end
    end
    return Rvec
end

function dm_diagonal_sum(dm::DM,salmon::SYSTEM)::ComplexF64
    nr=dm.nr
    v=ComplexF64(0)
    for ix3=1:nr[3],ix2=1:nr[2],ix1=1:nr[1]
        v+=dm.densitymatrix[ix1,ix2,ix3,ix1,ix2,ix3]
    end
    return v
end