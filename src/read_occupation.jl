function read_occupation(file::String,salmon::SYSTEM)
    # in SALMON code
    # dir_file_out = trim(odir)//"occupation.bin"
    # open(iu1_w,file=dir_file_out,form='unformatted')
    # write(iu1_w) system%rocc(1:system%no,1:system%nk,1:system%nspin)
    # close(iu1_w)
    nk=salmon.num_kgrid
    nb=salmon.nstate
    nk_all=nk[1]*nk[2]*nk[3]

    fio=FortranFile(file)
    occupation=zeros(nb,nk_all)
    # io=open(file,"r")
    # read!(io,occupation)
    # for ik=1:nk_all
    #     for ib=1:nb
    occupation=read(fio,(Float64,nb,nk_all))
        # end
    # end
    # close(io)

    return occupation
end