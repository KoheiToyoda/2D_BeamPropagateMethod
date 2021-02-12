
#コード分割のため、
#パラメーターはすべてここで定義

# struct の初期化のため、Parameters をインクルードする。

#negative value の累乗は.0を付ける

# 実験条件####################
module Params
   export  m
   export  mm
   export  um
   export  nm

    m = 1
    mm = 10^(-3.0m)
    um = 10^(-6.0m)
    nm = 10^(-9.0m)

    using Parameters

    #〇構造体の型注釈はフィールドにつけると速い。
    #その際、抽象型(Abstract type)ではなく具象型(Concrete type)を指定する。
    #〇Parameters.jlの構造体初期値機能は@with_kw
    @with_kw struct crange
        # セルの距離(um)
        x::Float64 = 20um
        z::Float64 = 20um
        t::Float64 = 10.0
    end
    @with_kw struct steps
        # セルのステップ（um
        x::Float64 = 0.1um
        z::Float64 = 0.1um
        t::Float64 = 0.1
    end

    # Int16 * Int16 = Int16 により、xとyの値によっては、
    # 配列作成時にあふれだしてしまう。
    @with_kw struct N
        x::Int64    
        z::Int64
        t::Int64
    end
    @with_kw struct materiarl
        nb::Float64 =1.5
        Δn0::Float64 = 0.03
        τ::Float64 = 0.1
        α::Float64 = 0.0
    end
    @with_kw struct beam
        w::Float64 = 3um
        U0::Float64 = 100.0
        wavelength::Float64 = 1.06um
    end

    @with_kw struct HG_mode
        m::Int8 = 0
        n::Int8 = 0
    end
    @with_kw struct LG_mode
        # userはl,pで代入するが、
        # m,nのフィールドも持っておき、コードライティング中では、m,nでやらせる。
        l::Int8 = 1
        p::Int8 = 0
        m = l
        n = p
    end
end