### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# ╔═╡ 4392a236-50f0-11eb-012f-f1cf6276de5a
using JuAFEM, SparseArrays

# ╔═╡ 6c6968e8-51f9-11eb-3862-8b0c81e071d5
include("../../definitions/fem.jl")

# ╔═╡ 4fd8b53c-50f0-11eb-30e7-035bff3609b0
begin
	"""
	    doassemble(cellvalues::CellValues, facevalues::FaceValues, dh::DofHandler) -> K, f
	
	Referenz JuAFEM Implementierung der Assemblierung
	"""
	function doassemble(cellvalues::CellValues, facevalues::FaceValues, dh::DofHandler) 
		n_basefuncs = getnbasefunctions(cellvalues)
		Ke = zeros(n_basefuncs, n_basefuncs)
		fe = zeros(n_basefuncs)

		K = create_sparsity_pattern(dh)
		f = zeros(ndofs(dh))
		assembler = start_assemble(K, f)

		b = Vec{2}((0,0))
		ℂ = elasticity(210,0.3)

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
					   ((cellcount, face) ∈ getfaceset(grid, "traction"))
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
	function elasticity(E, ν; dim=2)
		λ = E*ν / ((1 + ν) * (1 - 2ν))
		μ = E / (2(1 + ν))
		δ(i,j) = i == j ? 1.0 : 0.0
		f = (i,j,k,l) -> λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
		return SymmetricTensor{4,dim}(f)
	end
	
	"""
	create_cook_grid(nx, ny) -> JuAFEM.Grid

	# Argumente
	- `nx` Anzahl der Elemente in x-Richtung
	- `ny` Anzahl der Elemente in y-Richtung
	"""
	function create_cook_grid(nx, ny)
		corners = [Vec{2}((0.0,   0.0)),
				   Vec{2}((48.0, 44.0)),
				   Vec{2}((48.0, 60.0)),
				   Vec{2}((0.0,  44.0))]
		grid = generate_grid(Triangle, (nx, ny), corners);
		# facesets sind Markierungen denen wir Randbedingungen zuweisen können!
		addfaceset!(grid, "clamped", x -> norm(x[1]) ≈ 0.0);
		addfaceset!(grid, "traction", x -> norm(x[1]) ≈ 48.0);
		return grid
	end
	
	nothing
	
end

# ╔═╡ f2572c6a-50f2-11eb-3f00-e5c828b12b1f
md""" 
# Das eigene, erste Element
Willkommen in der siebten Übung im Modul Grundlagen der FEM. In dieser Übung werden sie ihr erstes Element selber implementieren. Dazu müssen wir zunächst wiederholen, was ein Element beschreibt und inwiefern es in der Assemblierung benötigt wird.
"""

# ╔═╡ dc9af4a8-5102-11eb-2aa1-9bf7df11d6ab
md"""
In der letzten Übung habe ich ihnen folgende Funktion vorgestellt
```julia
function doassemble(cellvalues::CellValues, facevalues::FaceValues, dh::DofHandler) 
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

	K = create_sparsity_pattern(dh)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
	
	b = Vec{2}((0,0))
	ℂ = elasticity(210,0.3)

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
                   ((cellcount, face) ∈ getfaceset(grid, "traction"))
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
```
### Verständnisfrage: Welche Funktionsaufrufe sind Element spezifisch?
Gehen sie die Assemblierungsfunktion genau durch und überlegen sie sich welche vier Funktionen Element spezifisch sind
"""

# ╔═╡ 105afeae-51ea-11eb-13ea-fbc35ac0a3e0
solution(md"""
1. `reinit!(cellvalues, cell)`
2. `getdetJdV(cellvalues, q_point)`
3. `shape_value(cellvalues, q_point, i)`
4. `shape_symmetric_gradient(cellvalues, q_point, i)`
	""")

# ╔═╡ 078a9c26-51ea-11eb-2d31-f7096973f6a9
md"""
JuAFEM ist modular aufgebaut und **nicht** Element zentrisch. Alte FE Software ist oftmals rund um die Element Entität aufgebaut und dadurch wird oftmals für User die direkte Berechnung der Elementsteifigkeitsmatrix als Schnitstelle angeboten.

!!! note "Hinweis"
    Im JuAFEM `doassemble` bauen wir die Steifigkeitsmatrix sukzessive über zwei Schleifen `i` bzw `j` auf
"""

# ╔═╡ d611296c-51c6-11eb-092c-2b9bd68e2115
md""" ## Lineares Dreieckselement"""

# ╔═╡ d7be5ef8-51c4-11eb-076b-1bf9c779fc89
html"""
<iframe width="700" height="400" src="https://www.youtube.com/embed/ASkg-oQOk8U" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ e90070e6-51c6-11eb-0dea-9b96c08dc498
md"""
Die Herleitung des Elements können sie im Video über diesen Text finden.

Da wir mit linearen Ansätzen arbeiten wird die $\mathbf{B}$ Matrix konstant sein, da sie die ersten Ableitungen der Ansatzfunktionen speichert. Wird nun weiterhin angenommen, dass das Material konstant linear elastisch für jedes Element ist, können wir also den Ausdruck

$$\mathbf{k}^e = \int_{E} (\mathbf{B}^e)^T\mathbb{C}\mathbf{B}^e \ \text{d} V$$

vereinfachen zu

$$\mathbf{k}^e = (\mathbf{B}^e)^T\mathbb{C}\mathbf{B}^e \int_{E} \text{d} V = A^e (\mathbf{B}^e)^T\mathbb{C}\mathbf{B}^e$$
$$[6\times 6] = [6 \times 3] [3 \times 3] [3 \times 6].$$

Hierbei ist $\mathbf{B}^e$ definiert durch

$$\mathbf{B}^e =
\begin{bmatrix}
N_{1,x} & 0 & N_{2,x} & 0 & N_{3,x} & 0 \\
0 & N_{1,y} & 0 & N_{2,y} & 0 & N_{3,y} \\
N_{1,y} & N_{1,x} & N_{2,y} & N_{2,x} & N_{3,y} & N_{3,x}
\end{bmatrix}$$

Die einzelnen Einträge ergeben sich hierbei aus den Knotenkoordinaten des Elements, nämlich

$$
\mathbf{N}_{,x} =
\begin{bmatrix}
N_{1,x} \\
N_{2,x} \\
N_{3,x} \\
\end{bmatrix} = 
\frac{1}{2A^e}
\begin{bmatrix}
y_2 - y_3 \\
y_3 - y_1 \\
y_1 - y_2 \\
\end{bmatrix}

\qquad 

\mathbf{N}_{,y} =
\begin{bmatrix}
N_{1,y} \\
N_{2,y} \\
N_{3,y} \\
\end{bmatrix} = 
\frac{1}{2A^e}
\begin{bmatrix}
x_3 - x_2 \\
x_1 - x_3 \\
x_2 - x_1 \\
\end{bmatrix}$$

Die Fläche des Elements kann über die Determinante der Matrix $\mathbf{A}^e$ bestimmt werden

$$
2A^e = \det \mathbf{A}^e = (x_2 - x_1)(y_3 - y_1) + (x_3 - x_1)(y_1 - y_2) \qquad \text{mit } \mathbf{A}^e = 
\begin{bmatrix}
1 & 1 & 1 \\
x_1 & x_2 & x_3\\
y_1 & y_2 & y_3
\end{bmatrix}$$
"""

# ╔═╡ 2e443d14-51e1-11eb-3a22-798a7427a384
begin
datatable=md"""
| Data       |      | values    |
| :------------- | :----------: | -----------: |
|  E |    | 210 GPa    |
|  ν |    | 0.3     |
| node no.1:    | [x y] | [0.0 0.0] mm |
| node no.2:    | [x y] | [1.0 0.3] mm |
| node no.3:    | [x y] | [0.5 0.8] mm |
"""
	nothing
end

# ╔═╡ 612589b8-51e1-11eb-1f0c-dd16569a600d
begin
	trianglesvg = HTML(open(f->read(f, String), "assets/triangleelement.svg"))
	nothing
end

# ╔═╡ 1487b332-51cc-11eb-26bf-db14dc9a1fdb
md"""
|     Fig.7.1   |            |
|---------------|------------|
|$(trianglesvg) |$(datatable)|
"""

# ╔═╡ ca281dee-50fa-11eb-31fb-a3a87f5ab886
begin 
	n1 = Node((0.0,0.0))
	n2 = Node((1.0, 0.3))
	n3 = Node((0.5,0.8))
	element = Triangle((1,2,3))
	testgrid = Grid([element],[n1,n2,n3])
end

# ╔═╡ 7d153014-510e-11eb-016e-cfa348898cf6
xe = getcoordinates(testgrid,1)

# ╔═╡ a1a60248-51ca-11eb-12f0-1926619a0a8e
md""" ### Verständnisfrage: Was berechnen sie in welcher Reihenfolge mit den Element Koordinaten?"""

# ╔═╡ acdc4702-51e9-11eb-16cb-6bb735781120
solution(md"""
1. Koordinaten extrahieren
2. Berechnung von $\det \mathbf{A}^e$
3. Gradienten der Ansatzfunktion $i$ bestimmten (also $\mathbf{N_{,x}}$,$\mathbf{N_{,y}}$)
4. Matrix $\mathbf{B}$ befüllen""")

# ╔═╡ ce9e717e-510f-11eb-24a5-e99ee402c102
function Bᵉ(coordinates)
	
	detAₑ = missing
	B = missing
	return B, detAₑ
end

# ╔═╡ b31f4640-511d-11eb-1125-fd87dba171cb
B , detA = Bᵉ(xe)

# ╔═╡ 39691180-51ce-11eb-35b7-c90c4f99f8a5
md"""
Hoffentlich stellen sie sich jetzt die Frage: Aber was hat das mit dem linearen Dreieck der letzten Übung zu tun?
**Sehr viel!**

Die große Frage ist, was die Spalten der $\mathbf{B}$ Matrix darstellen. Dafür initialisieren wir uns ein JuAFEM Element mit einem Gauss Punkt und linearen Ansatzfunktionen
"""

# ╔═╡ a20c82ae-50fb-11eb-248a-2f78ae6271cf
cellvalues = CellVectorValues(QuadratureRule{2, RefTetrahedron}(1), Lagrange{2, RefTetrahedron, 1}())

# ╔═╡ 6e430212-51ce-11eb-065e-1901668c52ca
md"""
Nun initalisiere ich die $\mathbf{B}$ Matrix und befülle sie.

### Verständnisfrage was sind die Spalten der B-Matrix?
"""

# ╔═╡ 574200aa-51e9-11eb-13b0-bf026ca3acd3
solution(md"""
Die Spalten der $\mathbf{B}$ Matrix sind die Voigt Notation des Ansatzfunktionsgradienten.

Oder anders gesagt: Unser Element hat 6 Ansatzfunktionen. Der Gradient dieser Ansatzfunktionen wird jeweils ein symmetrischer Tensor 2. Stufe sein. Da wir uns in 2D aufhalten (Dreieckselement) wird der Tensor über eine $4\times 4$ Koeffizientenmatrix beschreibbar sein, wobei Eintrag [1,2] identisch mit [2,1] ist (Symmetrie). Daher ist die Voigtnotation eines Symmetrischen Tensors 2. Stufe in 2 physikalischen Dimension ein Vektor mit **drei** Einträgen. Genau diese **drei** Einträge sind die $i$-te Spalte der B-Matrix.""")

# ╔═╡ 88e13b1a-510c-11eb-0efd-876526d8a2e8
begin
	dh = DofHandler(testgrid)
	push!(dh, :u, 2)
	close!(dh)
end

# ╔═╡ c6f7b842-510d-11eb-26db-63761118c696
for cell in CellIterator(dh)
	reinit!(cellvalues, cell)
end

# ╔═╡ d2a08698-5120-11eb-30dd-b1e6d2895721
getdetJdV(cellvalues,1)

# ╔═╡ 59f4ffa0-5105-11eb-36d0-f35c1b294d19
n_base_funcs = getnbasefunctions(cellvalues)

# ╔═╡ e31715ac-50fb-11eb-2a94-01de7b0801a4
begin
	Bu = zeros(3, n_base_funcs)
	for i in 1:getnbasefunctions(cellvalues) #
		Bu[:, i] += tovoigt(symmetric(shape_gradient(cellvalues, 1, i));offdiagscale=2.)
	end
end

# ╔═╡ 7aae873e-5105-11eb-0679-b35f1195afb0
Bu

# ╔═╡ afc850e0-5104-11eb-3f26-eba3a66750c4
tovoigt(elasticity(210,0.3))

# ╔═╡ b07af112-51ce-11eb-10af-7d221d36bccf
md"""
### Aufgabe 
Stellen sie die Elementsteifigkeitsmatrix des Dreiecks auf sowohl für ihre selbst geschriebene Funktion `Bᵉ` als auch für meine vorgegebene JuAFEM $\mathbf{B}$ Matrix. Überprüfen sie anschließend, ob die Matrizen übereinstimmen.
"""

# ╔═╡ 23e78adc-510e-11eb-25b7-318d2d556f39
missing * tovoigt(elasticity(210,0.3)) * missing * getdetJdV(cellvalues,1)

# ╔═╡ 011b869c-511e-11eb-120d-c1ff428aeb48
missing * tovoigt(elasticity(210,0.3)) * missing * missing

# ╔═╡ 4b9e22dc-5124-11eb-2d38-0daee470156a
try
		if B' * tovoigt(elasticity(210,0.3)) * B * (detA/2) ≈ Bu' * tovoigt(elasticity(210,0.3)) * Bu * getdetJdV(cellvalues,1)
		correct()
	else
		keep_working()
	end
catch
	still_missing()
end

# ╔═╡ ed687ef0-51e2-11eb-1b85-eb516ec769d7
md"""
### Aufgabe
Verifizieren sie ihr Element mithilfe der selben Konvergenzstudie der letzten Übung.

Vervollständigen sie die unten stehende `doassemble` Funktion und rufen sie an geeigneter Stelle ihre Funktion `Bᵉ` auf. Bilden sie anschließend die Elementsteifigkeitsmatrix. Als Referenz schauen sie sich die `doassemble` der vorherigen Übung an. Wir werden nun nicht mehr die Elementsteifigkeit sukzessive aufbauen, sondern durch ein Matrix-Matrix-Matrix Produkt
"""

# ╔═╡ 4303c606-51e5-11eb-21a5-7f21ea4e7560
function doassemble(facevalues::FaceValues, dh::DofHandler) 
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

	K = create_sparsity_pattern(dh)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
	
	b = Vec{2}((0,0))
	ℂ = elasticity(210,0.3)

    @inbounds for (cellcount,cell) in enumerate(CellIterator(dh)) #für Element in Elemente
        fill!(Ke, 0)
        fill!(fe, 0)
		
		#wow so leer hier
		
		
		#Ke .= 

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


# ╔═╡ d4122386-51e5-11eb-11c7-777a78f3f95a
meshsizes = [10 5
			 20 10
			 40 20
			 80 40
		     160 80
			 320 160]

# ╔═╡ b9ca16be-51e5-11eb-32c2-61787b20c1bb
function convergence_study(sizes)
	u_study = []
	ansatzfunktionen = Lagrange{2, RefTetrahedron, 1}()
	integration_rand = QuadratureRule{1, RefTetrahedron}(2)
	randvalues = FaceVectorValues(integration_rand, ansatzfunktionen)
	point = Vec(47.9,59.9)
	for meshsize in eachrow(meshsizes)
		grid = create_cook_grid(meshsize[1], meshsize[2])
		dh = DofHandler(grid)
		push!(dh, :u, 2)
		close!(dh)
		ch = ConstraintHandler(dh)
		∂Ω = getfaceset(grid, "clamped")
		dbc = Dirichlet(:u, ∂Ω, (x, t) -> [0,0],[1,2])
		add!(ch, dbc)
		close!(ch)
		update!(ch, 0.0)
		K, f = doassemble(randvalues, dh)
		apply!(K, f, ch)
		u = K \ f
		u_measured = measure_function(point,u,dh,ansatzfunktionen)
		push!(u_study, u_measured[2])
	end
	return u_study
end

# ╔═╡ cfae31c2-51e5-11eb-24eb-bd5b7b4b3907
study = convergence_study(meshsizes)

# ╔═╡ 3db9e0a6-51e8-11eb-2334-7f31e5e3d796
if isapprox(study,[7.81095,8.38392,8.595,8.67185,8.69987,8.70929], atol=1e-4)
	correct()
else
	keep_working()
end

# ╔═╡ 909e0ad2-51e7-11eb-1c7f-950597c1cf30
begin
	import Plots
	Plots.plot(study, label = "u am Punkt A", legend=:bottomright)
end

# ╔═╡ Cell order:
# ╠═4392a236-50f0-11eb-012f-f1cf6276de5a
# ╟─4fd8b53c-50f0-11eb-30e7-035bff3609b0
# ╟─6c6968e8-51f9-11eb-3862-8b0c81e071d5
# ╟─f2572c6a-50f2-11eb-3f00-e5c828b12b1f
# ╟─dc9af4a8-5102-11eb-2aa1-9bf7df11d6ab
# ╟─105afeae-51ea-11eb-13ea-fbc35ac0a3e0
# ╟─078a9c26-51ea-11eb-2d31-f7096973f6a9
# ╟─d611296c-51c6-11eb-092c-2b9bd68e2115
# ╟─d7be5ef8-51c4-11eb-076b-1bf9c779fc89
# ╟─e90070e6-51c6-11eb-0dea-9b96c08dc498
# ╟─2e443d14-51e1-11eb-3a22-798a7427a384
# ╟─1487b332-51cc-11eb-26bf-db14dc9a1fdb
# ╟─612589b8-51e1-11eb-1f0c-dd16569a600d
# ╠═ca281dee-50fa-11eb-31fb-a3a87f5ab886
# ╠═7d153014-510e-11eb-016e-cfa348898cf6
# ╟─a1a60248-51ca-11eb-12f0-1926619a0a8e
# ╟─acdc4702-51e9-11eb-16cb-6bb735781120
# ╠═ce9e717e-510f-11eb-24a5-e99ee402c102
# ╠═b31f4640-511d-11eb-1125-fd87dba171cb
# ╟─39691180-51ce-11eb-35b7-c90c4f99f8a5
# ╠═a20c82ae-50fb-11eb-248a-2f78ae6271cf
# ╟─6e430212-51ce-11eb-065e-1901668c52ca
# ╠═e31715ac-50fb-11eb-2a94-01de7b0801a4
# ╠═7aae873e-5105-11eb-0679-b35f1195afb0
# ╟─574200aa-51e9-11eb-13b0-bf026ca3acd3
# ╠═88e13b1a-510c-11eb-0efd-876526d8a2e8
# ╠═c6f7b842-510d-11eb-26db-63761118c696
# ╠═d2a08698-5120-11eb-30dd-b1e6d2895721
# ╠═59f4ffa0-5105-11eb-36d0-f35c1b294d19
# ╠═afc850e0-5104-11eb-3f26-eba3a66750c4
# ╟─b07af112-51ce-11eb-10af-7d221d36bccf
# ╠═23e78adc-510e-11eb-25b7-318d2d556f39
# ╠═011b869c-511e-11eb-120d-c1ff428aeb48
# ╟─4b9e22dc-5124-11eb-2d38-0daee470156a
# ╟─ed687ef0-51e2-11eb-1b85-eb516ec769d7
# ╠═4303c606-51e5-11eb-21a5-7f21ea4e7560
# ╠═b9ca16be-51e5-11eb-32c2-61787b20c1bb
# ╠═d4122386-51e5-11eb-11c7-777a78f3f95a
# ╠═cfae31c2-51e5-11eb-24eb-bd5b7b4b3907
# ╟─3db9e0a6-51e8-11eb-2334-7f31e5e3d796
# ╠═909e0ad2-51e7-11eb-1c7f-950597c1cf30
