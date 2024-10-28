function get_rxpx_grid(salmon::SYSTEM;prange::Vector{Float64},rrange::Vector{Float64},nr::Int64)
    return get_grid(salmon;prange=prange,rrange=rrange,nr=nr,rdim=1,pdim=1)
end
function get_rypy_grid(salmon::SYSTEM;prange::Vector{Float64},rrange::Vector{Float64},nr::Int64)
    return get_grid(salmon;prange=prange,rrange=rrange,nr=nr,rdim=2,pdim=2)
end
function get_rzpz_grid(salmon::SYSTEM;prange::Vector{Float64},rrange::Vector{Float64},nr::Int64)
    return get_grid(salmon;prange=prange,rrange=rrange,nr=nr,rdim=3,pdim=3)
end

function get_grid(salmon::SYSTEM;prange::Vector{Float64},rrange::Vector{Float64},nr::Int64,rdim::Int64=1,pdim::Int64=1)
    # Possible p point is included in {k+G} grid.
    # k is a Bloch wave vector
    # We assume that kx = n*Gx/m where n,m are positive integers.
    k=salmon.kvec
    nk_all=prod(salmon.num_kgrid)
    # find ekx
    klist=[k[ik][pdim] for ik=1:nk_all]
    # println("kx")
    # println(klist)
    dk=maximum(klist)
    tol=10^(-5)
    for ik=1:nk_all
        if abs(klist[ik]) < dk && abs(klist[ik]) > tol
            dk = abs(klist[ik])
        end
    end
    
    # dkx is now smallest non-zero value along kx axis
    fp1=mod(prange[1],dk)
    fp2=mod(prange[2],dk)
    
    p1=prange[1]+fp1
    p2=prange[2]-fp2
    np=Int64(round((p2-p1)/dk))
    
    # drx=abs(salmon.avec[1])/salmon.num_rgrid[1]
    dr=(rrange[2]-rrange[1])/nr
    println("Generate grid on (x,p) space.")
    println("dk ",dk)
    println("np ",np)
    println("p margin : fp1,fp2 ",fp1," ",fp2)
    println("dr ",dr)
    println("nr ",nr)
    xpgrid=[[zeros(3),zeros(3)] for ir=1:nr,ip=1:np]
    for ip=1:np
        p=p1+(ip-1)*dk
        for ir=1:nr
            r=rrange[1]+(ir-1)*dr

            xpgrid[ir,ip][1][rdim]=r
            xpgrid[ir,ip][2][pdim]=p
        end
    end
    return (grid=xpgrid, np=np)
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