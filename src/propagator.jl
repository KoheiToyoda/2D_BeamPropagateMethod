

#ビームの集光をテーパー型の屈折率分布として表現する。
function setNfocus(N::Array{ComplexF64,2},NA,n,separatingpoint,)
    
end

function setNwaveguide!(N::Array{ComplexF64,2},xstep, zstep, diameter, angle,baseN , propN, starting_separate)
    Nx = size(N,1)
    Nz = size(N,2)
    xrange_max = xstep * Nx / 2
    zrange_max = zstep * Nz
    xrange = range(-xrange_max,xrange_max,length = Nx)
    zrange = range(0,zrange_max,length = Nz)
    fill!(N,baseN)
    for (i,x_pos) in enumerate(xrange)
        for (k,z_pos) in enumerate(zrange)
            #セパレート前
#                if k < Nz * starting_separate

                if x_pos^2 < (diameter/2)^2 
                    N[i, k] = propN
                else 
                    N[i, k] = baseN
                end
                #斜め方向にして正射影を取ればいいと思うんだけど、
                #いまはまっすぐの導波路を作ろう。
                #if (x_pos-angle*x_pos)
            #end
            end
    end
end

function setNwaveguide2!(N, crange_x, crange_z, diameter, angle, baseN , propN, starting_separate)
    Nx = size(N,1)
    Nz = size(N,2)
    xrange = range(-crange_x/2, crange_x/2 , length = Nx)
    zrange = range(0, crange_z , length = Nz)
    @show xrange[1], xrange[end], xrange[Int64(floor(Nx / 2)+1)]
    fill!(N,baseN)
    for (i,x_pos) in enumerate(xrange)
        for (k,z_pos) in enumerate(zrange)
            #セパレート前
#                if k < Nz * starting_separate
                if z_pos < starting_separate
                    if (x_pos)^2 <= (diameter/2)^2 
                        N[i, k] = propN
                    else 
                        N[i, k] = baseN
                    end
                    #斜め方向にして正射影を取ればいいと思うんだけど、
                    #いまはまっすぐの導波路を作ろう。
                    #if (x_pos-angle*x_pos)
                #end
                end
                if z_pos >= starting_separate
                    if (x_pos - tan(deg2rad(angle))*(z_pos-starting_separate))^2 <= (diameter/2)^2
                        N[i, k] = propN
                    elseif (x_pos + tan(deg2rad(angle))*(z_pos-starting_separate))^2 <= (diameter/2)^2
                        N[i, k] = propN
                    else
                        N[i, k] = baseN
                    end
                end

            end
    end
end
