function read_salmon_input(file::String,file_kpt)::SYSTEM
    io=open(file,"r")
    doc=readlines(io)    
    close(io)
    nline=length(doc)
    salmon=SYSTEM()
    for il=1:nline
        line=doc[il]
        words=split(line,"=")
        # println(words)
        nw=length(words)
        if nw == 2
            option=strip(words[1])
            if option == "theory"
                v = split(words[2],"'") .|> strip           
                salmon.theory = v[2]
            elseif option == "sysname"
                v = split(words[2],"'") .|> strip
                salmon.sysname = v[2]
            elseif option == "al(1:3)"
                v=replace(words[2],"d0" => "")
                println(v)
                v=split(v,",")
                salmon.alat=parse.(Float64,strip.(v))
            elseif option == "nelem"
                salmon.nelem=parse(Int64,words[2])
            elseif option == "natom"
                salmon.natom=parse(Int64,words[2])
                salmon.atomic_red_coor = [zeros(3) for ia=1:salmon.natom]
            elseif option == "nelec"
                salmon.nelec=parse(Int64,words[2])
            elseif option == "nstate"
                salmon.nstate=parse(Int64,words[2])
            elseif option == "num_kgrid(1:3)"
                v=split(words[2],",")
                salmon.num_kgrid=parse.(Int64,v)
            elseif option == "num_rgrid(1:3)"
                v=split(words[2],",")
                salmon.num_rgrid=parse.(Int64,v)
            end
        end
    end
    iflag=0
    for il=1:nline
        line=doc[il]
        word=strip(line)
        if word == "&atomic_red_coor"
            iflag = il
        end
    end
    for il= iflag+1:nline
        line=doc[il]
        v=strip(line)
        # if the line is not a comment, start to read atomic coordinate
        # println(typeof(v[1]))
        if !occursin("!","$(v[1])") 
            # println(il)
            for ia = 1:salmon.natom
                line=doc[il+ia-1]
                v=split(line)
                # println(v)
                for id=1:3
                    salmon.atomic_red_coor[ia][id]=parse(Float64,v[1+id])
                end
            end
            break
        end
    end

    # Generate lattice vector
    for id=1:3
        v=zeros(3)
        v[id]=salmon.alat[id]
        salmon.avec[id]=v        
    end
    salmon.bvec=get_bvec(salmon.avec)

    # Generate r grid
    # dr::r point interval
    dr=[zeros(3) for id=1:3]
    for id=1:3
        dr[id]=salmon.avec[id] ./ salmon.num_rgrid[id]
    end
    # println(dr)
    nr=salmon.num_rgrid
    nr_all=nr[1]*nr[2]*nr[3]
    # println(nr_all)
    salmon.rvec=[zeros(3) for ir=1:nr_all]
    ir=1
    for ir3=0:nr[3]-1
        for ir2=0:nr[2]-1
            for ir1=0:nr[1]-1
                salmon.rvec[ir]=ir1*dr[1]+ir2*dr[2]+ir3*dr[3]
                ir+=1
            end
        end
    end

    # read k grid
    io=open(file_kpt,"r")
    doc=readlines(io)
    close(io)
    nk=salmon.num_kgrid
    nk_all=nk[1]*nk[2]*nk[3]
    salmon.kvec=[zeros(3) for ik=1:nk_all]
    nline=length(doc)
    ik=1
    for il=6:nk_all+5
        v=split(doc[il])
        salmon.kvec[ik]=parse.(Float64,[v[2],v[3],v[4]])
        ik+=1        
    end

    return salmon
end

function get_bvec(avec::Vector{Vector{Float64}})
    bvec=[zeros(3) for id=1:3]
    volume=cross(avec[1],avec[2])'*avec[3]
    bvec[1]=2pi .* ( cross(avec[2],avec[3]) ./volume)
    bvec[2]=2pi .* ( cross(avec[3],avec[1]) ./volume)
    bvec[3]=2pi .* ( cross(avec[1],avec[2]) ./volume)
    return bvec
end

function read_wfn(file::String,salmon::SYSTEM)::WFN
    # doc=FortranFiles.FortranFile(file,"r",access="direct",recl=56623104;rec=1)
    io=open(file,"r")
    # FortranFiles.read(doc,Float64)
    # println(doc)


    nk=salmon.num_kgrid
    nr=salmon.num_rgrid
    nb=salmon.nstate
    nk_all=nk[1]*nk[2]*nk[3]
    nr_all=nr[1]*nr[2]*nr[3]
    
    # volume element
    a=salmon.avec
    V=abs(a[1]'*cross(a[2],a[3]))
    dV=V/prod(nr)

    wfn=WFN()
    wfn.nk=nk
    wfn.nr=nr
    wfn.nb=nb
    wfn.orbital=zeros(ComplexF64,nr[1],nr[2],nr[3],nb,nk_all)
    println("Read orbital function from $file")
    for im=1:1
    for ik=1:nk_all
        t=now()
        if mod(ik,10)==0
            println("kpt $ik/$nk_all  ",hour(t),":",minute(t),":",second(t))
        end
        for ib=1:nb
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
    end
    end 
    # wfn.orbital=dV*wfn.orbital
    close(io)
    t=now()
    println("Finished reading orbitals.  ",t)
    return wfn
end