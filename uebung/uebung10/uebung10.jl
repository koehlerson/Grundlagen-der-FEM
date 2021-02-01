### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 18227574-5675-11eb-1296-4bce242be471
using JuAFEM, SparseArrays

# ╔═╡ 34d635de-5675-11eb-334c-097e04f6abe7
begin
	using PlutoUI
	import Plots
	function doassemble(cellvalues::CellVectorValues, facevalues::FaceValues, dh::DofHandler,E,ν;t=Vec{2}((0,5))) 
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

    K = create_sparsity_pattern(dh)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    
    b = Vec{2}((0,0))
    ℂ = elasticity(E,ν)

    @inbounds for (cellcount,cell) in enumerate(CellIterator(dh)) #für Element in Elemente
        fill!(Ke, 0)
        fill!(fe, 0)
        reinit!(cellvalues, cell)

        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)

            for i in 1:n_basefuncs
                δu  = shape_value(cellvalues, q_point, i)
                δε = shape_symmetric_gradient(cellvalues, q_point, i)
                fe[i] += (δu ⋅ b) * dΩ
                for j in 1:n_basefuncs
                    ε = shape_symmetric_gradient(cellvalues,q_point,j)
                    Ke[i, j] += (δɛ ⊡ ℂ ⊡ ɛ) * dΩ
                end
            end
        end
        for face in 1:nfaces(cell)
            if onboundary(cell, face) &&
                   ((cellcount, face) ∈ getfaceset(dh.grid, "traction"))
                reinit!(facevalues, cell, face)
                for q_point in 1:getnquadpoints(facevalues)
                    dΓ = getdetJdV(facevalues, q_point)
                    for i in 1:n_basefuncs
                        δu = shape_value(facevalues, q_point, i)
                        fe[i] += (δu ⋅ t) * dΓ
                    end
                end
            end
        end
        assemble!(assembler, celldofs(cell), fe, Ke)
    end
    return K,f
end
		function create_cook_grid(nx, ny, t::Union{Type{QuadraticTriangle},Type{Triangle},Type{Quadrilateral}})
		corners = [Vec{2}((0.0,   0.0)),
				   Vec{2}((48.0, 44.0)),
				   Vec{2}((48.0, 60.0)),
				   Vec{2}((0.0,  44.0))]
		grid = generate_grid(t, (nx, ny), corners);
		# facesets sind Markierungen denen wir Randbedingungen zuweisen können!
		addfaceset!(grid, "clamped", x -> norm(x[1]) ≈ 0.0);
		addfaceset!(grid, "traction", x -> norm(x[1]) ≈ 48.0);
		return grid
	end
		function elasticity(E, ν; dim=2)
		λ = E*ν / ((1 + ν) * (1 - 2ν))
		μ = E / (2(1 + ν))
		δ(i,j) = i == j ? 1.0 : 0.0
		f = (i,j,k,l) -> λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
		return SymmetricTensor{4,dim}(f)
	end
	nothing
end

# ╔═╡ 8aa3fb16-622a-11eb-2938-ef508e6d76f2
include("../../definitions/plotting.jl")

# ╔═╡ 30471402-5675-11eb-176c-5330058c812b
include("../../definitions/fem.jl")

# ╔═╡ 3dfad08e-5675-11eb-2f13-591a06338472
md"""
# Lineares Viereckselement

Willkommen in der zehnten Übung von Grundlagen der FEM. In dieser Übung beschäftigen wir uns mit dem linearen Viereckselement. In der letzten Übung haben wir die Ansatzfunktion des quadratischen Dreiecks sowie `reinit!` zusammen implementiert. Wir haben uns Automatisches Ableiten zu Nutze gemacht, was wir in dieser Übung nicht verwenden werden, sondern wir werden die klassische Methode verwenden. Mit klassischer Methode ist hier gemeint, dass wir Ansatzfunktionen sowie deren Ableitung selber programmieren. Anschließend werden wir Viereckselement und Dreieckselement miteinander vergleichen.
"""

# ╔═╡ 93afe74e-5675-11eb-2b88-ad6ef3c67a10
HTML(open(f->read(f, String), "assets/quadrilateral.svg"))

# ╔═╡ b72c3244-5676-11eb-101f-a92e00229685
integration = QuadratureRule{2, RefCube}(2)

# ╔═╡ baa8a842-5678-11eb-1bcc-5d31818d4a2e
integration.weights

# ╔═╡ bfd21d4e-5678-11eb-00b1-c530b298c084
integration.points

# ╔═╡ 17a0f648-573f-11eb-2632-c70ddf47ed59
struct LinearesViereck <: Interpolation{2,RefCube,1} end

# ╔═╡ 4469ac5c-5752-11eb-09d3-cb65c9b15076
md"""
### Aufgabe
Implementieren sie `value` und `gradient` für das lineare Viereckselement
"""

# ╔═╡ b2854bcc-5743-11eb-053e-a7e81fd4c948
md"""
Welche Ansatzfunktion soll graphisch geprüft werden? $(@bind ifunc Slider(1:4,default=2))
"""

# ╔═╡ 00b273ee-5744-11eb-248c-2d75533976fe
md"""
Ansatzfunktion $(Int(ifunc)) wird angezeigt
"""

# ╔═╡ 1560cfda-573f-11eb-1113-8338e62772e5
function JuAFEM.value(ip::LinearesViereck, i::Int, ξ::Vec{2})
    0
end

# ╔═╡ 0631202a-5744-11eb-01f6-53f629bd00bd
begin
	ξ₁ = collect(-1:0.1:1)
	ξ₂ = collect(-1:0.1:1)
	Plots.surface(ξ₁,ξ₂, (x1,x2) -> JuAFEM.value(LinearesViereck(), ifunc, Vec{2}((x1,x2))), xlabel="ξ₁", ylabel="ξ₂", zlabel="value") 
end

# ╔═╡ f186b6da-5743-11eb-1ba0-d79ed268d133
function JuAFEM.gradient(ip::LinearesViereck, i::Int, ξ::Vec{2})
	[0,0]
end

# ╔═╡ fe6f3e4a-5751-11eb-0c51-2f65e72047ae
md"""
### Verständnisfrage: Welche Punkte müssen nun mit `value` und `gradient` ausgewertet werden?

$(solution(md"Alle Gaußpunkte werden ausgewertet und in `CellVectorValues` gespeichert",blur=true))

"""

# ╔═╡ 91c71c7a-573f-11eb-3f00-3d380bd2288f
begin
	JuAFEM.getnbasefunctions(::LinearesViereck) = 4
	JuAFEM.nvertexdofs(::LinearesViereck) = 1
	JuAFEM.faces(::LinearesViereck) = ((1,2), (2,3), (3,4), (4,1))
	function JuAFEM.reference_coordinates(::LinearesViereck)
    return [Vec{2, Float64}((-1.0, -1.0)),
            Vec{2, Float64}(( 1.0, -1.0)),
            Vec{2, Float64}(( 1.0,  1.0,)),
            Vec{2, Float64}((-1.0,  1.0,))]
	end
	
	#JuAFEM.getngeobasefunctions(cv::CellVectorValues) = size(cv.M, 1)
	#JuAFEM.getn_scalarbasefunctions(cv::CellScalarValues) = size(cv.N, 1)
	JuAFEM.getdim(ip::LinearesViereck) = 2
	JuAFEM.getrefshape(ip::LinearesViereck) = RefCube
	nothing
end

# ╔═╡ ec71ec04-5749-11eb-2959-c373ff53eded
cv_juafem = CellVectorValues(QuadratureRule{2,RefCube}(2), Lagrange{2,RefCube,1}())

# ╔═╡ 724920ae-574a-11eb-0cc5-fbf634c07728
cv_juafem.N

# ╔═╡ 78e94b28-574a-11eb-1c7a-2fa5b0548997
cv_implementation = CellVectorValues(QuadratureRule{2,RefCube}(2), LinearesViereck())

# ╔═╡ f2756446-574e-11eb-0a9d-356169379bf5
begin
	if cv_juafem.N ≈ cv_implementation.N && cv_juafem.dNdξ ≈ cv_implementation.dNdξ && cv_juafem.M ≈ cv_implementation.M && cv_juafem.dMdξ ≈ cv_implementation.dMdξ
		correct()
	else
		keep_working()
	end
end

# ╔═╡ 74953d70-574c-11eb-2900-8955e4d656de
cv_implementation.N

# ╔═╡ 7f341230-5752-11eb-1ead-17968440015e
md"""
Wie wird nun also `CellVectorValues` initialisiert? Dafür habe ich ihnen den Konstruktur unten angegeben. 

### Verständnisfrage
Welche Schleifen sind nötig um `CellVectorValues` zu füllen?

"""

# ╔═╡ 07901bd6-5746-11eb-09e8-4fea99c0b68b
function JuAFEM.CellVectorValues(::Type{T}, quad_rule::QuadratureRule{dim,shape}, func_interpol::LinearesViereck,
        geom_interpol::Interpolation=LinearesViereck()) where {dim,T,shape<:JuAFEM.AbstractRefShape}
	
    n_qpoints = length(JuAFEM.getweights(quad_rule))

    # Function interpolation
    n_func_basefuncs = getnbasefunctions(func_interpol) * dim
    N    = fill(zero(Vec{dim,T})      * T(NaN), n_func_basefuncs, n_qpoints)
    dNdx = fill(zero(Tensor{2,dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)
    dNdξ = fill(zero(Tensor{2,dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)

    # Geometry interpolation
    n_geom_basefuncs = getnbasefunctions(geom_interpol)
    M    = fill(zero(T)          * T(NaN), n_geom_basefuncs, n_qpoints)
    dMdξ = fill(zero(Vec{dim,T}) * T(NaN), n_geom_basefuncs, n_qpoints)

    for (qp, ξ) in enumerate(quad_rule.points)
        basefunc_count = 1
        for basefunc in 1:getnbasefunctions(func_interpol)
            N_temp = JuAFEM.value(func_interpol, basefunc, ξ)
			dNdξ_temp = JuAFEM.gradient(func_interpol, basefunc, ξ)
            for comp in 1:dim
                N_comp = zeros(T, dim)
                N_comp[comp] = N_temp
                N[basefunc_count, qp] = Vec{dim,T}((N_comp...,))

                dN_comp = zeros(T, dim, dim)
                dN_comp[comp, :] .= dNdξ_temp
                dNdξ[basefunc_count, qp] = Tensor{2,dim,T}((dN_comp...,))
                basefunc_count += 1
            end
        end
        for basefunc in 1:n_geom_basefuncs
			M[basefunc, qp] = JuAFEM.value(geom_interpol, basefunc, ξ)
			dMdξ[basefunc, qp] = Vec{2}(JuAFEM.gradient(geom_interpol, basefunc, ξ))
        end
    end

    detJdV = fill(T(NaN), n_qpoints)
    MM = Tensors.n_components(Tensors.get_base(eltype(dNdx)))

    JuAFEM.CellVectorValues{dim,T,shape,MM}(N, dNdx, dNdξ, detJdV, M, dMdξ, quad_rule.weights)
end

# ╔═╡ 2b896730-5756-11eb-381b-8f76acd3347b
meshsizes = [10 5
			 20 10
			 40 20
			 80 40
		     160 80
			 320 160]

# ╔═╡ 2ae066a8-5756-11eb-03a1-0b98318f3ca4
function convergence_study(sizes, interpolation, E, ν)
	u_study = []
	if isa(interpolation, Lagrange{2,RefTetrahedron,1})
		meshcelltype = Triangle
		integration = QuadratureRule{2, RefTetrahedron}(1)
		elementvalues = CellVectorValues(integration, interpolation)
		integration_rand = QuadratureRule{1, RefTetrahedron}(2)
		randvalues = FaceVectorValues(integration_rand, interpolation)
	elseif isa(interpolation, Lagrange{2,RefTetrahedron,2})
		meshcelltype = QuadraticTriangle
		integration = QuadratureRule{2, RefTetrahedron}(3)
		elementvalues = CellVectorValues(integration, interpolation)
		integration_rand = QuadratureRule{1, RefTetrahedron}(2)
		randvalues = FaceVectorValues(integration_rand, interpolation)
	elseif isa(interpolation, LinearesViereck) || isa(interpolation, Lagrange{2,RefCube,1})
		meshcelltype = Quadrilateral
		integration = QuadratureRule{2, RefCube}(2)
		elementvalues = CellVectorValues(integration, interpolation)
		integration_rand = QuadratureRule{1, RefCube}(2)
		randvalues = FaceVectorValues(integration_rand, Lagrange{2,RefCube,1}())
	end

	point = Vec(47.9,59.9)
	for meshsize in eachrow(meshsizes)
		grid = create_cook_grid(meshsize[1], meshsize[2], meshcelltype)
		dh = DofHandler(grid)
		push!(dh, :u, 2)
		close!(dh)
		ch = ConstraintHandler(dh)
		∂Ω = getfaceset(grid, "clamped")
		dbc = Dirichlet(:u, ∂Ω, (x, t) -> [0,0],[1,2])
		add!(ch, dbc)
		close!(ch)
		update!(ch, 0.0)
		K, f = doassemble(elementvalues, randvalues, dh, E, ν)
		apply!(K, f, ch)
		u = K \ f
		u_measured = measure_function(point,u,dh,interpolation)
		push!(u_study, [ndofs(dh),u_measured[2]])
	end
	return u_study
end

# ╔═╡ 5ffda758-575c-11eb-17f8-bf9ec208cc9c
E = 210

# ╔═╡ 62268e84-575c-11eb-1762-7f9984695b68
ν = 0.3

# ╔═╡ 6bb363ba-5756-11eb-2ab2-61d6d262ccf4
lineartrianglestudy = convergence_study(meshsizes, Lagrange{2,RefTetrahedron,1}(), E, ν)

# ╔═╡ 7556574c-5756-11eb-1572-2730d14e0984
quadratictrianglestudy = convergence_study(meshsizes,Lagrange{2,RefTetrahedron,2}(),E, ν)

# ╔═╡ 9ca12ab4-575d-11eb-2b70-0f9b1b91090f
linearquad = convergence_study(meshsizes, LinearesViereck(), E, ν)

# ╔═╡ 6701644c-575c-11eb-2587-13c64a2e9bb9
begin
	Plots.plot(getindex.(linearquad,1), getindex.(linearquad,2),xaxis=:log, label="lineares Viereck", marker="hex")
	Plots.plot!(getindex.(lineartrianglestudy,1), getindex.(lineartrianglestudy,2),xaxis=:log, label="lineares Dreieck", marker="hex")
	Plots.plot!(getindex.(quadratictrianglestudy,1), getindex.(quadratictrianglestudy,2),xaxis=:log, label="quadratisches Dreieck",legend=:bottomright,xlabel="Anzahl Freiheitsgrade", ylabel="Vertikale Verschiebung an Punkt A", marker="hex")
end

# ╔═╡ Cell order:
# ╠═18227574-5675-11eb-1296-4bce242be471
# ╟─8aa3fb16-622a-11eb-2938-ef508e6d76f2
# ╟─30471402-5675-11eb-176c-5330058c812b
# ╟─34d635de-5675-11eb-334c-097e04f6abe7
# ╟─3dfad08e-5675-11eb-2f13-591a06338472
# ╟─93afe74e-5675-11eb-2b88-ad6ef3c67a10
# ╠═b72c3244-5676-11eb-101f-a92e00229685
# ╠═baa8a842-5678-11eb-1bcc-5d31818d4a2e
# ╠═bfd21d4e-5678-11eb-00b1-c530b298c084
# ╠═17a0f648-573f-11eb-2632-c70ddf47ed59
# ╟─4469ac5c-5752-11eb-09d3-cb65c9b15076
# ╟─0631202a-5744-11eb-01f6-53f629bd00bd
# ╟─b2854bcc-5743-11eb-053e-a7e81fd4c948
# ╟─00b273ee-5744-11eb-248c-2d75533976fe
# ╠═1560cfda-573f-11eb-1113-8338e62772e5
# ╠═f186b6da-5743-11eb-1ba0-d79ed268d133
# ╟─f2756446-574e-11eb-0a9d-356169379bf5
# ╟─fe6f3e4a-5751-11eb-0c51-2f65e72047ae
# ╟─91c71c7a-573f-11eb-3f00-3d380bd2288f
# ╠═ec71ec04-5749-11eb-2959-c373ff53eded
# ╠═724920ae-574a-11eb-0cc5-fbf634c07728
# ╠═78e94b28-574a-11eb-1c7a-2fa5b0548997
# ╠═74953d70-574c-11eb-2900-8955e4d656de
# ╟─7f341230-5752-11eb-1ead-17968440015e
# ╠═07901bd6-5746-11eb-09e8-4fea99c0b68b
# ╠═2ae066a8-5756-11eb-03a1-0b98318f3ca4
# ╠═2b896730-5756-11eb-381b-8f76acd3347b
# ╠═5ffda758-575c-11eb-17f8-bf9ec208cc9c
# ╠═62268e84-575c-11eb-1762-7f9984695b68
# ╠═6bb363ba-5756-11eb-2ab2-61d6d262ccf4
# ╠═7556574c-5756-11eb-1572-2730d14e0984
# ╠═9ca12ab4-575d-11eb-2b70-0f9b1b91090f
# ╟─6701644c-575c-11eb-2587-13c64a2e9bb9
