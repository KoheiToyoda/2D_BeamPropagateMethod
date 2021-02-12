### このコードは２D用！！！

#計算に関係ないパッケージ
using Pkg
include("propagator.jl")
include("Params.jl")
#計算に使うパッケージ
using LinearAlgebra
using Plots
using KTOptical 

# 自作のパッケージを使うときにはGithubにアップロードして、
# GithubのHTTPSをadd HTTPS_URLする。

using CSV
using DataFrames
using BenchmarkTools
using JLD2
using Polynomials
um = Params.um

#using FileIO
# 計算条件###################
#計算レンジ
λ = 1.06um
crange = Params.crange(x = 30λ, z = 2500λ, t = 0.1)

#計算ステップ
#steps = Params.steps(x = 0.2λ, z = 1λ, t = 0.1)
Nx = 501#Int(floor(crange.x / steps.x))
Nz = 2501#Int(floor(crange.z / steps.z))
Nt = 1#Int(floor(crange.t / steps.t))
N = Params.N(Nx,Nz,Nt)

steps = Params.steps(x = crange.x/Nx, z = crange.z/Nz, t = 0.1)

#材料情報
nb = 1.47
Δn0 = -0.16
Nref = nb + Δn0/3
mode = Params.HG_mode(0,0)
core_diameter = 3λ
beam_diameter = 1λ

material = Params.materiarl(nb = nb, Δn0 = Δn0, τ = 0.1,α = 0)
#ビーム情報
beam = Params.beam(w = beam_diameter, U0 = 1, wavelength = 1.064um)

#println("計算環境")
#versioninfo()
#Pkg.status()

println("計算条件")
@show crange
@show steps
@show N
#@show material
@show mode
@show beam

#println("//////////////////////////////")
function returnPades(T::Integer,N::Integer)
    # Julia に型注釈は必要ないが、Padeの数は整数型に限定したい。という説明
    # 未実装、広角(>10°)の場合はエラーが出る
end

# 各Zに対して、calcstepを回していく。
#
# Ax = B をつくってトーマスアルゴリズムで解く
# juliaの場合はB\Aでよい。(バックスラッシュ！）
# Aが N.x＊N.x(z=k+1の係数)
# すなわち、(x,z) = (i,k)において、
#
# BがN.xサイズのベクトル(z=kの係数)
# 
###### 各A,Bに入る屈折率項は座標情報が必要である。
# 
#F_k_before 既知のビーム伝搬
#F_k_after 未知のビーム伝搬の格納先(F_k+1)
#k 今のZカウント
# 藪哲郎光導波路132Pより
#

function calcStep1!(F_k_before, k, matN, Nref)
    F_k_after = similar(F_k_before)
    k0 = 2π / beam.wavelength
    # a,b,cは左辺用
    a = -1/(steps.x)^2
    c = -1/(steps.x)^2

    # d は 右辺用
    B = zeros(ComplexF64,N.x,1)

    # 透明境界条件 藪哲郎光導波路
    alpha = F_k_before[1] / F_k_before[2]
    beta = F_k_before[N.x] / F_k_before[N.x-1]

    # 各回 iを固定して、jxj　の行列を作って計算する。
    # よって、最初のループは固定方向のx(i)で、
    # 内部で、diagmで行列方向のN.y x N.y = jxjを作る。

    # Ax=BのA (z = k+1 における係数行列)を作る。
    # mapでベクトルを作っておいて、diagmで対角行列にする。
    # bが位置によって値が異なるのでmapで対応。
    # diagmは配列を正方行列の対角に配置する。
    # pair（=>) でズレを表現できる。

    b(i,k)     = 1im*4*k0*((matN[i, k]+matN[i+1, k])/2)/steps.z + 2/steps.x^2 - k0^2(((matN[i, k]+matN[i+1, k])/2)^2 - Nref^2)
    b_end(i,k) = 1im*4*k0*((matN[i, k]+matN[i-1, k])/2)/steps.z + 2/steps.x^2 - k0^2(((matN[i, k]+matN[i-1, k])/2)^2 - Nref^2)
    
    Apre = map( i -> begin
        if i == N.x
            return b_end(i, k)
        else
            return b(i, k)
        end
    end, 
    1:N.x)

    A = diagm(Apre)
    A += diagm(1 => fill(a,N.x-1))
    A += diagm(-1 => fill(a,N.x-1))


        #反転処理
    if real(-(1/(1im * steps.x)) * log(F_k_before[1] / F_k_before[2])) < 0
        #realパートの符号反転処理がないと反射する。
        tr_kl =  imag(-(1/(1im * steps.x)) * log(F_k_before[1] / F_k_before[2]))
        tr_kl -= real(-(1/(1im * steps.x)) * log(F_k_before[1] / F_k_before[2]))
        alpha = exp(1im * tr_kl * -steps.x)
    end

    if real(-(1/(1im * steps.x)) * log(F_k_before[N.x] / F_k_before[N.x-1])) < 0
        tr_kr =  imag(-(1/(1im * steps.x)) * log(F_k_before[N.x] / F_k_before[N.x-1]))
        tr_kr -= real(-(1/(1im * steps.x)) * log(F_k_before[N.x] / F_k_before[N.x-1]))
        beta = exp(-1im * tr_kr * steps.x)
    end

    A[1,1] += a * alpha
    A[N.x,N.x] += c* beta 

    e = 1/steps.x^2
    F_left = alpha * F_k_before[1]
    F_right = beta * F_k_before[N.x]

    d(i, k)      = e * F_k_before[i-1]  + ((4im*Nref*k0)/steps.z - 2e + (((matN[i, k]+matN[i+1, k]) / 2)^2 - Nref^2)*k0^2)* F_k_before[i] + e * F_k_before[i+1]
    d_left(i,k)  = e * F_left  +          ((4im*Nref*k0)/steps.z - 2e + (((matN[i, k]+matN[i+1, k]) / 2)^2 - Nref^2)*k0^2)* F_k_before[i] + e * F_k_before[i+1]
    d_right(i,k) = e * F_k_before[i-1]  + ((4im*Nref*k0)/steps.z - 2e + (((matN[i, k]+matN[i-1, k]) / 2)^2 - Nref^2)*k0^2)* F_k_before[i] + e * F_right

    # https://docs.julialang.org/en/v1/manual/functions/#Do-Block-Syntax-for-Function-Arguments
    B = map( i -> begin
                        if i == 1
                            return d_left(i,k)
                        elseif i == N.x
                            return d_right(i,k)
                        else
                            return d(i,k)
                        end
                end, 
                1:N.x)
    #@show F_k_before

    #########################

    F_k_after = B\A
    #F_k_after[]を使って、新しい屈折率マップを作る。
    return F_k_after
end

function retN()
    return 1
end

# 試験的。BPMのモード測定を使う
function showMode()
end

# 電界の初期条件
function initial_set(Mode, Beamparam)
    KTOptical.setParam(Beamparam.w, 0, Beamparam.wavelength)

    x = range(-crange.x/2, crange.x/2, length = N.x)
    #y = range(-crange.y/2, crange.y/2 - steps.x ,length = N.y)

    #ここをLG_E HG_Eののなんかラッパ的な奴にできない？
    E = HG_E.(Mode.m, Mode.n, x)
    #plot()
    return E
end

function renewN!(matN, E3d)
end

function main()
    #セル個数
    F_result = zeros(ComplexF64, N.x, N.z)
    #@show size(F0)
    
    E = initial_set(mode, beam)
    #E2 = initial_set
    gr()
    #!!!!!!!!!!!!!!!!!!!!!!!!
    x = range(-crange.x/2, crange.x/2 ,length = N.x)
    z = range(0, crange.z ,length = N.z)

    #F_k_1stは現在のF_k
    #F_k_2ndは更新されたF_k
    @show N.x, N.x
    F_k_1st = zeros(ComplexF64, N.x) #Zerosに型指定忘れないこと！！！
    matN = zeros(ComplexF64,N.x, N.z)
    #@show size(F0)
    @show size(E)
    #@show size(F_k_1st[:,1])
    F_k_1st = E
    F_result[:,1] .= E
#    function setNwaveguide2!(N, crange_x, crange_z, diameter, angle, baseN , propN, starting_separate)

    setNwaveguide2!(matN, crange.x, crange.z, core_diameter, 0.2, material.nb, material.nb + material.Δn0, 600λ)
    #@show matN
    
    #左辺は更新されたF_k_1stが入る。
    # S
    #初期条件を F_K_1st に入れる。
    
    for t in 1:N.t
        @show "Zmax:", N.z
        for k in 1:N.z
            # x 固定、　y方向移動
            if k%100 == 0
                @show t,k
            end
            F_k_1st = F_result[:,k] = calcStep1!(F_k_1st, k, matN, Nref)
        end
    #    @save "/savefile/F_"* string(t)*".jld2" F_k_2nd
    end
    # @show F_result
    println("calc done")
    println("plotting.....") 
    Ezx = abs.(F_result[: , :])
#  @show Ezx
    @show typeof(Ezx)
    angle = [30,70]
    p1 = plot(x, z, Ezx',
      st=:surface, camera = (angle[1] ,80))
    display(p1)
    p3 = plot(x,z,Ezx',
    st=:surface, camera = (0,90))
    display(p3)
    p4 = plot(x, abs.(E),zlim=(0,0.75),cmax = (0,0.75))
    display(p4)
    print()
    p5 = plot(x,z, abs.(matN)',
    st=:surface, camera = (0,90))
    display(p5)
    p6 = plot(x, abs.(matN[:,1]))
    display(p6)
    #//////////////////////////////////////////////////////////////////////

    tEzx = DataFrame(Ezx)
    tMatN = DataFrame(matN)
    tEzx |> CSV.write("test2.csv" , delim = "," , writeheader = false)
    tMatN |> CSV.write("test_MatN.csv" , delim = "," , writeheader = false)
end

@time main()

 