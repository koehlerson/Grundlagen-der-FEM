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

# ╔═╡ b481d88c-55b0-11eb-20c0-b5a183424ced
using JuAFEM, SparseArrays

# ╔═╡ d97d9764-55dd-11eb-23aa-25d63bf69cdb
begin
	using PlutoUI
	import Plots
	function doassemble(cellvalues::CellVectorValues, facevalues::FaceValues, dh::DofHandler,E,ν) 
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
                t = Vec{2}((0,5))
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
		function create_cook_grid(nx, ny, t::Union{Type{QuadraticTriangle},Type{Triangle}})
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

# ╔═╡ c4fbe160-55f1-11eb-2bf1-076572f7c8ee
include("../../definitions/fem.jl")

# ╔═╡ d4d866f8-55b0-11eb-29d6-fbf5c8ed7012
md"""
# Quadratisches Dreieckselement

Willkommen in der neunten Übung von Grundlagen der FEM. In dieser Übung beschäftigen wir uns mit dem quadratischen Dreieckselement. In Übung 7 haben sie kennen gelernt, wie man anhand der $\mathbf{B}$ Matrix sein eigenes Element implementieren kann. Heute beschäftigen wir uns mit der Implementierung des quadratischen Dreieckselements basierend auf der Ansatzfunktion Schleife $i$ und $j$. Wir folgen also in dieser Übung eher den Weg von JuAFEM. Da sie in der letzten Übung viele Ansatzfunktionen schon erfolgreich aufstellen konnten, müssen wir heute nur wiederholen und an der richtigen Stelle eingreifen.
"""

# ╔═╡ 78a2a2de-55b7-11eb-3d71-8901db10496c
begin
	quadratictriangle = HTML(open(f->read(f, String), "assets/quadratictriangle.svg"))
	nothing
end

# ╔═╡ 0719861c-55bc-11eb-3cef-a33cb2fff055
begin
datatable = md"""
| Data       |      | values    |
| :------------- | :----------: | -----------: |
|  E |    | 100 GPa    |
|  ν |    | 0.3     |
| node 1:    | [x y] | [0.0 0.0] mm |
| node 2:    | [x y] | [1.0 0.3] mm |
| node 3:    | [x y] | [0.5 0.8] mm |
| node 4:    | [x y] | [0.5 0.15] mm |
| node 5:    | [x y] | [0.75 0.55] mm |
| node 6:    | [x y] | [0.25 0.4] mm |
"""
	nothing
end

# ╔═╡ 223f35e6-55bd-11eb-333f-dfd8ec7b54e4
md"""
$(quadratictriangle)
$(datatable)"""

# ╔═╡ 4770c46a-55bd-11eb-1d91-37681d0435fe
begin 
	n1 = Node((0.0,0.0))
	n2 = Node((1.0, 0.3))
	n3 = Node((0.5,0.8))
	n4 = Node((0.5,0.15))
	n5 = Node((0.75, 0.55))
	n6 = Node((0.25, 0.4))
	element = QuadraticTriangle((1,2,3,4,5,6))
	testgrid = Grid([element],[n1,n2,n3,n4,n5,n6])
end

# ╔═╡ 80d75f84-55bd-11eb-0c96-6f563a2d7c4f
struct QuadratischesDreieck  <: Interpolation{2,RefTetrahedron,2} end

# ╔═╡ 0e0478fe-55cd-11eb-31dc-174be39419c5
md"""
Die Ansatzfunktionen des quadratischen Dreiecks lauten

$$
\mathbf{N} =
\begin{bmatrix}
N_{1} \\
N_{2} \\
N_{3} \\
N_{4} \\
N_{5} \\
N_{6}
\end{bmatrix} = 
\begin{bmatrix}
\lambda_1(2\lambda_1 - 1) \\
\lambda_2(2\lambda_2 - 1) \\
\lambda_3(2\lambda_3 - 1) \\
4 \lambda_1 \lambda_2 \\
4 \lambda_2 \lambda_3\\
4 \lambda_3 \lambda_1
\end{bmatrix}, \quad \text{\small wobei sich die Flächenkoordinaten} \quad 
\begin{bmatrix}
\lambda_1 \\
\lambda_2 \\
\lambda_3 
\end{bmatrix}=
\begin{bmatrix}
\xi_1 \\
\xi_2 \\
1 - \xi_1 - \xi_2 
\end{bmatrix}$$
durch die Parameterkoordinaten $\xi_1$ und $\xi_2$ ausdrücken lassen. 
"""

# ╔═╡ e13f84aa-55d3-11eb-2d5a-9fe94f40d221
begin
	physikalischekoordinaten = html"""
<span style="color: blue">x,y,z</span>
"""
	nothing
end

# ╔═╡ 404a2f1e-55d6-11eb-2956-5797d3a45630
md"""
## Aufgabe 
Füllen sie die Funktion `value`. Ziel ist es auf unser Struct `QuadratischesDreieck` die Funktion zu dispatchen. Die Funktion wird später bspw. von `shape_value` aufgerufen. Eingabewerte der Funktion sind zum einen unsere Interpolation `ip::QuadratischesDreieck`, `i::Int` - welche Basisfunktion? und die Parameterkoordinaten als zwei dimensionaler Vektor `ξ::Vec{2}`
"""

# ╔═╡ b32f3e1e-55dd-11eb-2215-65c5daa82be4
md"""
Welche Ansatzfunktion soll graphisch geprüft werden? $(@bind ifunc Slider(1:6,default=3))
"""

# ╔═╡ 170693ce-55de-11eb-17cc-ab39215510b9
md"""
Ansatzfunktion $(Int(ifunc)) wird angezeigt
"""

# ╔═╡ b2165910-55c7-11eb-24d6-b3c5625a09a4
function JuAFEM.value(ip::QuadratischesDreieck, i::Int, ξ::Vec{2})
	0
end

# ╔═╡ cd61ffda-55db-11eb-295c-298e9be0f34c
begin
	ξ₁ = collect(0:0.01:1)
	ξ₂ = collect(0:0.01:1)
	Plots.surface(ξ₁,ξ₂, (x1,x2) -> (1 - x1 - x2) > 0.0 ? JuAFEM.value(QuadratischesDreieck(), ifunc, Vec{2}((x1,x2))) : 0.0, xlabel="ξ₁", ylabel="ξ₂", zlabel="value") 
end

# ╔═╡ 7a83d192-565f-11eb-0c7d-97feeb3d78c1
begin
	test1 = JuAFEM.value(QuadratischesDreieck(), 1, Vec{2}((0.5,0))) == 0 && JuAFEM.value(QuadratischesDreieck(), 1, Vec{2}((1,0))) == 1 && JuAFEM.value(QuadratischesDreieck(), 1, Vec{2}((0,0))) == 0
	test2 = JuAFEM.value(QuadratischesDreieck(), 2, Vec{2}((0,0.5))) == 0 && JuAFEM.value(QuadratischesDreieck(), 2, Vec{2}((0,1))) == 1 && JuAFEM.value(QuadratischesDreieck(), 2, Vec{2}((0,0))) == 0
	test3= JuAFEM.value(QuadratischesDreieck(), 3,Vec{2}((0,0.5))) == 0 && JuAFEM.value(QuadratischesDreieck(), 3,Vec{2}((0,1))) == 0 && JuAFEM.value(QuadratischesDreieck(), 3,Vec{2}((0,0))) == 1
	test4= JuAFEM.value(QuadratischesDreieck(), 4, Vec{2}((0.5,0.5))) == 1 && JuAFEM.value(QuadratischesDreieck(), 4, Vec{2}((0,1))) == 0 && JuAFEM.value(QuadratischesDreieck(), 4, Vec{2}((1,0))) == 0
	test5= JuAFEM.value(QuadratischesDreieck(), 5, Vec{2}((0.0,0.5))) == 1 && JuAFEM.value(QuadratischesDreieck(), 5, Vec{2}((0,1))) == 0 && JuAFEM.value(QuadratischesDreieck(), 5, Vec{2}((1,0))) == 0
	test6= JuAFEM.value(QuadratischesDreieck(), 6, Vec{2}((0.5,0.0))) == 1 && JuAFEM.value(QuadratischesDreieck(), 6, Vec{2}((0,1))) == 0 && JuAFEM.value(QuadratischesDreieck(), 6, Vec{2}((1,0))) == 0
	if test1 && test2 && test3 && test4 && test5 && test6
		correct()
	else
		keep_working()
	end
end

# ╔═╡ c85736ae-55d3-11eb-3bec-d5430291117c
md"""
Mithilfe der

1. ersten Ableitungen der Ansatzfunktionen nach den Parameterkoordinaten
2. den physikalischen Knotenkoordinaten und 
3. dem isoparametrischen Konzept 
lässt sich die Jacobimatrix aufstellen.

Die Jacobimatrix stellt den Zusammenhang zwischen **physikalischen** ($physikalischekoordinaten) und **parametrischen** ($$\xi_1,\xi_2,\xi_3$$) Raum her.
"""

# ╔═╡ 288d6f10-55bf-11eb-07ba-cd5e0be8e4bd
begin
	JuAFEM.getnbasefunctions(::QuadratischesDreieck) = 6
	JuAFEM.nfacedofs(::QuadratischesDreieck) = 1
	JuAFEM.faces(::QuadratischesDreieck) = ((1,2,4), (2,3,5), (3,1,6))
	JuAFEM.nvertexdofs(::QuadratischesDreieck) = 1
	JuAFEM.vertices(::QuadratischesDreieck) = (1,2,3)
	function JuAFEM.reference_coordinates(::QuadratischesDreieck)
		return [Vec{2, Float64}((1.0, 0.0)),
				Vec{2, Float64}((0.0, 1.0)),
				Vec{2, Float64}((0.0, 0.0)),
				Vec{2, Float64}((0.5, 0.5)),
				Vec{2, Float64}((0.0, 0.5)),
				Vec{2, Float64}((0.5, 0.0))]
	end
	
	JuAFEM.getngeobasefunctions(cv::CellVectorValues) = size(cv.M, 1)
	JuAFEM.getn_scalarbasefunctions(cv::CellScalarValues) = size(cv.N, 1);
end

# ╔═╡ e0c18ac8-55d6-11eb-030b-b3683e485a25
md"""
Wir haben aber nur die Ansatzfunktion implementiert. Wie sollen wir jetzt also die Ableitung davon bekommen?

Dafür gibt es **zwei** Wege. Eine klassische und eine moderne Methode.

### Die klassische Methode
Sie implementieren eine Funktion die ihnen die Ableitungen zurück gibt. Sie leiten also per Hand ab und definieren sich im Code eine Funktion die zur Ableitung assoziiert wird. Das haben wir mehr oder weniger in Übung 7 bei der Implementierung des linearen Dreiecks gemacht. Die B-Matrix ist nämlich genau diese Ableitung nur in besonderer Weise als Matrix aufgebaut.

### Die moderne Methode
Sie lassen ihren PC die Arbeit tun. Schlaue Mathematiker und Informatiker haben sich überlegt, wie man Funktionen (im Informatik bzw Coding Sinne) ableitet ohne ihre Ableitung anzugeben (bzw. zu kennen). Die Idee ist relativ nahe liegend: Wir wissen wie man per Hand ableitet sofern eine Funktion gegeben ist, wieso also auch nicht das den PC übernehmen lassen. Beispiel:
"""

# ╔═╡ a6b5566e-55c9-11eb-3cab-7be38a79f69c
dNdξ, N = gradient(ξ -> JuAFEM.value(QuadratischesDreieck(), 4, ξ), Vec{2}((0.5,0.5)), :all)

# ╔═╡ 810d5ca0-55d7-11eb-20a8-83614d7c5f90
md"""
Hier habe ich die 4. Ansatzfunktion nach Koordinatenvektor $\xi$ abgeleitet und die Ableitung an der Stelle (0.5,0.5) ausgewertet. Zur Erinnerung: Die 4. Ansatzfunktion war gegeben durch:
$$N_4 = 4\xi_1 \ \xi_2$$
Diese Funktionalität nennt man Automatic Differentiation (abgekürzt AD, automatisches ableiten)
"""

# ╔═╡ dc558e7a-55d7-11eb-3342-e9397159dd82
md"""
In dieser Übung werden wir die moderne Methode benutzen, da sie uns einige Arbeit spart und wir dafür unser Augenmerk auf die Jacobimatrix richten können. Die Jacobimatrix ist wie folgt definiert:

$$\mathbf{J} = \sum_{I=1}^{n}
\begin{bmatrix}
 x_I M_{I,\xi_1}   &  x_I M_{I,\xi_2}   \\
 y_I M_{I,\xi_1}   &  y_I M_{I,\xi_2}
\end{bmatrix}=
\begin{bmatrix}
 x_1 & x_2 & ... & x_{n}  \\
 y_1 & y_2 & ... & y_{n}
\end{bmatrix}
\begin{bmatrix}
 M_{1,\xi_1} & M_{1,\xi_2} \\
 M_{2,\xi_1} & M_{2,\xi_2} \\
 \vdots    & \vdots     \\
 M_{n,\xi_1} & M_{n,\xi_2}
\end{bmatrix}$$

!!! note "Zusatzinformation"
    In Tensorprodukt Schreibweise entspricht das 
    $$\mathbf{x}_i \otimes \mathbf{M}_{i,\mathbf{\xi}}$$

Hier notiert $M$ die geometrische Interpolation. Diese ist für isoparametrische Elemente die gleiche wie für unsere Verschiebung (die wir mit N notieren).

Also ist die Jacobimatrix als **äußeres** (dyadisches, $\otimes$) Produkt aus Knotenkoordinaten und Ableitung der Ansatzfunktion nach den parametrischen Koordinaten zu verstehen. Sie ist eine extrem wichtige Entität, nicht nur im Mechanik Kontext.

Ziel ist es die Inverse der Matrix zu haben, denn dadurch können wir folgende Relation berechnen:

$$
 \mathbf{N}_{I,x} = \begin{bmatrix}
N_{Ix,x} & N_{Iy,x}\\ 
N_{Ix,y} & N_{Iy,y}
\end{bmatrix}= 
\begin{bmatrix}
 N_{Ix,\xi_1} & N_{Iy,\xi_1}\\ 
N_{Ix,\xi_2} & N_{Iy,\xi_2}
\end{bmatrix}\ \mathbf{J}^{-1}$$

!!! note "Zusatzinformation"
    In JuAFEM ist die Ansatzfunktion eines Freiheitsgrad (hier examplisch von Knoten 	1) wie folgt definiert:

    $$\mathbf{N}_1 = \begin{bmatrix}
    N^1(x, y)\\
    0
    \end{bmatrix}
    \quad
    \mathbf{N}_2 = \begin{bmatrix}
    0\\
    N^1(x, y)
    \end{bmatrix}$$
    
    Wobei der Index von $\mathbf{N}$ bis 12 läuft (12 Freiheitsgrade) und der Index von $N$ bis 6 (6 Ansatzfunktionen für jeweils eine Richtung).

Also die Ableitung der Ansatzfunktion nach der physikalischen Koordinate (x,y). Genau diese Ableitung brauchen wir in der Elementsteifigkeitsmatrix (siehe schwache Form).

Neben der Inverse der Jacobi, möchten wir auch noch die Determinante bestimmen, denn dieser Wert kommt bei der Auflösung des Integrals vor.

$$\int\limits_\Omega f(\mathbf{x}) d \Omega \approx \sum\limits_{q = 1}^{n_q} f(\mathbf{x}_q) \det \mathbf{J} \ w_q$$

Die Anteile die durch das Integral auftauchen sind also zum einen die Determinante von $\mathbf{J}$ und zum anderen die Wichtung des Integrationspunkt $w_q$. Daher der Funktionsname `getdetJdV` in JuAFEM.
"""

# ╔═╡ 8e967584-55e9-11eb-3b33-1b2efe208352
md"""
### Aufgabe
Vervollständigen sie `reinit!`. Ziel ist es die Werte in `cv::CellVectorValues` zu überschreiben bzw. zu aktualisieren. Dazu eine Übersicht, der relevanten Variablen von `cv`

**CellVectorValues**
- `cv.qr_weights` Vektor aller Wichtungen; Eintrag $i$ entspricht $w_i$
- `cv.detJdV` Vektor der pro Gaußpunkt $\det \mathbf{J} \ w_q$ hält; Eintrag $i$ entspricht $\det \mathbf{J} \ w_i$
- `cv.dMdξ[j,i]` Ableitung der geometrischen Ansatzfunktionen ($\mathbf{M}$). Zeile entspricht der $j$-ten Ansatzfunktion, Spalte entspricht $i$-ten Gausspunkt
- `cv.dNdξ[j,i]` Ableitung der Verschiebung Ansatzfunktionen nach der **parametrischen** Koordinate. Zeile entspricht der $j$-ten Ansatzfunktion, Spalte entspricht $i$-ten Gausspunkt
- `cv.dNdx[j,i]` Ableitung der Verschiebung Ansatzfunktionen nach der **physikalischen Koordinate**. Zeile entspricht der $j$-ten Ansatzfunktion, Spalte entspricht $i$-ten Gausspunkt
- `x::AbstractVector{Vec{dim,T}}` selbe Struktur wie `xe`
"""

# ╔═╡ c98cb286-5efc-11eb-1e7e-5f4d68433267
cellvalues_implementation = CellVectorValues(QuadratureRule{2, RefTetrahedron}(3), QuadratischesDreieck())

# ╔═╡ 4b7add68-5665-11eb-23c0-7fc310fcbea9
xe = getcoordinates(testgrid, 1)

# ╔═╡ 191c38d2-5663-11eb-0a5e-67907ff0294b
cellvalues_implementation.qr_weights

# ╔═╡ db0894dc-5662-11eb-3b21-a144fb7ba097
cellvalues_implementation.detJdV

# ╔═╡ d292ca02-5662-11eb-2b77-cd5f5722b543
cellvalues_implementation.dNdξ[1,1]

# ╔═╡ ccf5fbae-5668-11eb-0285-65f8ceed4970
cellvalues_implementation.dMdξ[1,1]

# ╔═╡ 722363ca-55e5-11eb-3cc0-2befb87163f2
function JuAFEM.reinit!(cv::CellVectorValues{dim}, x::AbstractVector{Vec{dim,T}}) where {dim,T}
    n_geom_basefuncs = JuAFEM.getngeobasefunctions(cv)
    n_func_basefuncs = JuAFEM.getn_scalarbasefunctions(cv)
    @assert length(x) == n_geom_basefuncs
    isa(cv, CellVectorValues) && (n_func_basefuncs *= dim)


    @inbounds for i in 1:length(cv.qr_weights) #Loop über alle Integrationspunkte

    end
end

# ╔═╡ dd2fce0a-55eb-11eb-09c9-1f68f18bdb69
begin
	dh_implementation_test = DofHandler(testgrid)
	push!(dh_implementation_test, :u, 2)
	close!(dh_implementation_test)
	quad_order = 3
	cellvalues_juafem = CellVectorValues(QuadratureRule{2, RefTetrahedron}(quad_order), Lagrange{2,RefTetrahedron,2}())
	try
	for cell in CellIterator(dh_implementation_test)
		reinit!(cellvalues_juafem, cell)
		reinit!(cellvalues_implementation, cell)
	end
	catch
		keep_working()
	end
	if getdetJdV(cellvalues_implementation,1) ≈ getdetJdV(cellvalues_juafem,1) && cellvalues_implementation.dNdx ≈ cellvalues_juafem.dNdx
		correct()
	else
		keep_working()
	end
end

# ╔═╡ 9b9148e6-55f2-11eb-34e8-994c039d472b
meshsizes = [10 5
			 20 10
			 40 20
			 80 40
		     160 80
			 320 160]

# ╔═╡ f033c83e-55f1-11eb-0837-59ee3ba3f4d4
function convergence_study(sizes, interpolation, E, ν)
	u_study = []
	integration = QuadratureRule{2, RefTetrahedron}(3)
	elementvalues = CellVectorValues(integration, interpolation)
	integration_rand = QuadratureRule{1, RefTetrahedron}(2)
	randvalues = FaceVectorValues(integration_rand, interpolation)
	point = Vec(47.9,59.9)
	for meshsize in eachrow(meshsizes)
		meshcelltype = Triangle
		isa(interpolation, Lagrange{2,RefTetrahedron,1}) ? (meshcelltype = Triangle) : (meshcelltype = QuadraticTriangle)
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

# ╔═╡ c06aacce-5756-11eb-2991-dfcd2c36cf3e
E = 210

# ╔═╡ c4004c72-5756-11eb-33b9-a55883aab05e
ν = 0.3

# ╔═╡ a03e7f44-55f2-11eb-2237-cf19909b0ea8
linearstudy = convergence_study(meshsizes, Lagrange{2,RefTetrahedron,1}(), E, ν)

# ╔═╡ 477a5ffa-564b-11eb-0d31-6118dc420d9c
quadraticstudy = convergence_study(meshsizes,QuadratischesDreieck(), E, ν)

# ╔═╡ 429e6ffc-564c-11eb-185c-7556eb50d474
begin
	Plots.plot(getindex.(linearstudy,1), getindex.(linearstudy,2),xaxis=:log, label="lineares Dreieck", marker="hex")
	Plots.plot!(getindex.(quadraticstudy,1), getindex.(quadraticstudy,2),xaxis=:log, label="quadratisches Dreieck",legend=:bottomright,xlabel="Anzahl Freiheitsgrade", ylabel="Vertikale Verschiebung an Punkt A", marker="hex")
end

# ╔═╡ Cell order:
# ╠═b481d88c-55b0-11eb-20c0-b5a183424ced
# ╟─c4fbe160-55f1-11eb-2bf1-076572f7c8ee
# ╟─d97d9764-55dd-11eb-23aa-25d63bf69cdb
# ╟─d4d866f8-55b0-11eb-29d6-fbf5c8ed7012
# ╟─78a2a2de-55b7-11eb-3d71-8901db10496c
# ╟─0719861c-55bc-11eb-3cef-a33cb2fff055
# ╟─223f35e6-55bd-11eb-333f-dfd8ec7b54e4
# ╠═4770c46a-55bd-11eb-1d91-37681d0435fe
# ╠═80d75f84-55bd-11eb-0c96-6f563a2d7c4f
# ╟─0e0478fe-55cd-11eb-31dc-174be39419c5
# ╟─e13f84aa-55d3-11eb-2d5a-9fe94f40d221
# ╟─404a2f1e-55d6-11eb-2956-5797d3a45630
# ╟─cd61ffda-55db-11eb-295c-298e9be0f34c
# ╟─7a83d192-565f-11eb-0c7d-97feeb3d78c1
# ╟─b32f3e1e-55dd-11eb-2215-65c5daa82be4
# ╟─170693ce-55de-11eb-17cc-ab39215510b9
# ╠═b2165910-55c7-11eb-24d6-b3c5625a09a4
# ╟─c85736ae-55d3-11eb-3bec-d5430291117c
# ╟─288d6f10-55bf-11eb-07ba-cd5e0be8e4bd
# ╟─e0c18ac8-55d6-11eb-030b-b3683e485a25
# ╠═a6b5566e-55c9-11eb-3cab-7be38a79f69c
# ╟─810d5ca0-55d7-11eb-20a8-83614d7c5f90
# ╟─dc558e7a-55d7-11eb-3342-e9397159dd82
# ╟─8e967584-55e9-11eb-3b33-1b2efe208352
# ╠═c98cb286-5efc-11eb-1e7e-5f4d68433267
# ╟─4b7add68-5665-11eb-23c0-7fc310fcbea9
# ╠═191c38d2-5663-11eb-0a5e-67907ff0294b
# ╠═db0894dc-5662-11eb-3b21-a144fb7ba097
# ╠═d292ca02-5662-11eb-2b77-cd5f5722b543
# ╠═ccf5fbae-5668-11eb-0285-65f8ceed4970
# ╠═722363ca-55e5-11eb-3cc0-2befb87163f2
# ╟─dd2fce0a-55eb-11eb-09c9-1f68f18bdb69
# ╠═f033c83e-55f1-11eb-0837-59ee3ba3f4d4
# ╠═9b9148e6-55f2-11eb-34e8-994c039d472b
# ╠═c06aacce-5756-11eb-2991-dfcd2c36cf3e
# ╠═c4004c72-5756-11eb-33b9-a55883aab05e
# ╠═a03e7f44-55f2-11eb-2237-cf19909b0ea8
# ╠═477a5ffa-564b-11eb-0d31-6118dc420d9c
# ╟─429e6ffc-564c-11eb-185c-7556eb50d474
