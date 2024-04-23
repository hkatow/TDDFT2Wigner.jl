function get_px_rx_plane(salmon::SYSTEM,prange::Vector{Float64},rrange::Vector{Float64};ngrid::Int64)
    k=salmon.kvec
    nk_all=prod(salmon.num_kgrid)
    # find ekx
    dkx=0.0
    for ik=1:nk_all
        if k[1] > 0 && k[2] == 0 && k[3] == 0
            (dkx)
        end
    end
    return
end

# function plot_density(dm::DM,salmon::SYSTEM,dir::String,interpolation::Int64)
#     nr=dm.nR
#     rho=[real(dm.densitymatrix[ix1,ix2,ix3,ix1,ix2,ix3]) for ix1=1:nr[1],ix2=1:nr[2],ix3=1:nr[3]]
#     rho_1=zeros(Float64,nr[2],nr[3])
#     rho_2=zeros(Float64,nr[1],nr[3])
#     rho_3=zeros(Float64,nr[1],nr[2])

    
#     for ix3=1:nr[3], ix2=1:nr[2]
#         v1=0.0;v2=0.0;v3=0.0
#         for ix1=1:nr[1]
#             v1+=rho[ix1,ix2,ix3]
#             v2+=rho[ix2,ix1,ix3] 
#             v3+=rho[ix2,ix3,ix1]             
#         end
#         rho_1[ix2,ix3]=v1
#         rho_2[ix2,ix3]=v2
#         rho_3[ix2,ix3]=v3
#     end
    
#     # plot 
#     ni=interpolation
#     mkpath(dir)
#     ENV["GKSwstyle"]="null"
#     gr()
#     plt=plot()
#     rho_dense=resample(rho_1,(ni*nr[2],ni*nr[3]))
#     contour!(plt,rho_dense,fill=true,title="rho_1",legend=false)
#     savefig(plt,"$dir/rho_1.pdf")

#     plt=plot()
#     rho_dense=resample(rho_2,(ni*nr[1],ni*nr[3]))
#     contour!(plt,rho_dense,fill=true,title="rho_2",legend=false)
#     savefig(plt,"$dir/rho_2.pdf")

#     plt=plot()
#     rho_dense=resample(rho_2,(ni*nr[1],ni*nr[2]))
#     contour!(plt,rho_dense,fill=true,title="rho_3",legend=false)
#     savefig(plt,"$dir/rho_3.pdf")

#     return
# end