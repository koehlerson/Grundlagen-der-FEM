### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# ╔═╡ fa20e096-48f1-11eb-1c40-f7c48ff8b2a4
using JuAFEM, SparseArrays

# ╔═╡ fce316a4-50ed-11eb-335d-7d24c4032a43
using Plots

# ╔═╡ 3940b9e4-4e75-11eb-29ec-4d1da14f089c
begin
	include("../../definitions/fem.jl")
	include("../../definitions/plotting.jl")
end

# ╔═╡ ee0b0c68-48ed-11eb-35e0-d513bd952297
md""" 
# FEM mit einem Software Framework - JuAFEM
Willkommen in der sechsten Übung im Modul Grundlagen der FEM. In dieser Übung werden wir uns mit der FE Toolbox JuAFEM beschäftigen. Im speziellen mit einem typischen Problem: Cook's membrane. Das Randwertproblem ist im Bild unter diesem Text gegeben.
"""

# ╔═╡ 23b1189e-48ee-11eb-3652-37e022de561e
HTML(open(f->read(f, String), "assets/cooksmembrane.svg"))

# ╔═╡ 1dc8fe14-48f1-11eb-37fb-ab620a99c16c
md"
Grundsätzlich werde ich sie interaktiv durch den Ablauf von JuAFEM führen und sie werden im Anschluss das Gelernte nutzen, um eine Konvergenzstudie mit der Verschiebung im Punkt A durchzuführen.
"

# ╔═╡ 6efe578a-4e98-11eb-33e3-d7edbb880dd5
md"""
## Vorbereitungen
"""

# ╔═╡ 5ded3c3e-48f2-11eb-0883-5ddf347145ac
md"
Zunächst gebe ich ihnen eine Funktion an die Hand, die das Mesh parametrisch anhand von `nx` und `ny` für sie erstellt. 
"

# ╔═╡ 2dfa3f64-48f2-11eb-1c16-0b6477748f43
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

# ╔═╡ baf3f6b6-4e86-11eb-03a1-8d8d66bb944e
md"
Einer der wichtigsten Sachen in der FEM sind die gewählten Ansatzfunktionen.
In JuAFEM wählt man sie über

```julia
Lagrange{dim, RefShape, order}()
```

wobei `dim` für die physikalische Dimension, `RefShape` für die Topologie (tetrahedrische `RefTetrahedron` oder quaderartige `RefCube`) und `order` für den Ansatzgrad stehen. 
"

# ╔═╡ 262a7aee-48f7-11eb-0db3-3be29e5a9d6b
ansatzfunktionen = Lagrange{2, RefTetrahedron, 1}()

# ╔═╡ a6be058e-4e88-11eb-0dfa-5110b45ac242
md"""
Eine nächste wichtige Auswahl ist, wie wir die Integrale auflösen. Rufen sie sich noch einmal die Gleichungen für die Elementsteifigkeitsmatrix in das Gedächtnis und achten sie genau auf die Stelle, an der wir Integrale haben. In Übung 5 haben wir diese über analytische Integration gelöst.

$\boldsymbol{\tilde{k}}^e = \int_0^l {\boldsymbol{B}^e}^T EA \boldsymbol{B}^e d\tilde{x}$
	
$\boldsymbol{B}^e = [-\frac{1}{l} \quad \frac{1}{l}] \quad \text{const in } \tilde{x}$
	
$\Rightarrow \boldsymbol{\tilde{k}}^e = EA {\boldsymbol{B}^e}^T\boldsymbol{B}^e\int_0^l d\tilde{x}$
	
$\boldsymbol{\tilde{k}}^e = \frac{EA}{L}\begin{bmatrix}
	1 & -1 \\
	-1 & 1 \end{bmatrix}$

(Übung 5 für mehr Details.)

Das ist jedoch **nicht generell möglich**. In der vorherigen Übung nutzten wir den Fakt, dass $E,A$ und $\mathbf{B}$ konstant sind über die Integrationsvariable. 

Dieser Fakt kann schon durch eine andere Wahl von Ansatzgrad bzw. veränderlichen Querschnitt verletzt werden. Daher sollte man etwas flexibleres nutzen und das ist numerische Integration (auch numerische Quadratur). Wir werden im Rahmen dieser Veranstaltung die Gauss Quadratur nutzen.

Mithilfe numerischer Quadratur kann mein Integral wie folgt umschreiben

$$\int_{\Omega} f(\mathbf{x}) \ \text{d}\Omega \approx \ \sum_{q=1}^{n_q} f(\mathbf{x}_q)w_q$$

wobei für Polynome das $\approx$ zur Gleichheit wird, sofern genügend $n_q$ Punkte verwendet werden.

Generell spricht man bei den Punkten $\mathbf{x}_q$ von den Integrationspunkten (Materialpunkten) und bei $w_q$ von den Wichtungen (engl. weights).
"""

# ╔═╡ c9abb7a8-4e8a-11eb-1ae5-11d947d575f4
md"""
JuAFEM stellt diese Quadratur für uns bereit unter folgender Syntax

```julia
QuadratureRule{dim, RefShape}(order)
```

Mehr theoretische Details werden sie in noch zu folgenden Übungen erfahren.

"""

# ╔═╡ 2e611812-48f7-11eb-29b9-c5a3e1651bfa
integration = QuadratureRule{2, RefTetrahedron}(1)

# ╔═╡ 62061634-4e8f-11eb-187c-177d497addd4
md"""
`CellValues`, im speziellen für uns `CellVectorValues`, sind der Dreh und Angelpunkt in JuAFEM. Sie stellen die **Elementwerte** eines **Elements** dar. Sie halten die Ansatzfunktionen (`ansatzfunktionen`) sowie (`integration`), also die ausgewerteten Funktionen ausgewertet im Integrationspunkt.
"""

# ╔═╡ 39d110c4-48f7-11eb-3da9-37d44c9c73f3
elementvalues = CellVectorValues(integration, ansatzfunktionen)

# ╔═╡ ca537870-4e90-11eb-0487-c7159de71eeb
md"""
Da wir eine Last am Rand haben, müssen wir auch über den Rand integrieren und brauchen hierfür **Randwerte** die für die jweiligen Elemente gelten, die am Rand des Körpers liegen.
"""

# ╔═╡ d9a7d7c0-490f-11eb-3c59-470d0a89297a
integration_rand = QuadratureRule{1, RefTetrahedron}(2)

# ╔═╡ 04a20cfc-4910-11eb-15cd-31526230c2f6
randvalues = FaceVectorValues(integration_rand, ansatzfunktionen)

# ╔═╡ f3f77000-4e90-11eb-38d2-ff24eeabd75d
md" Nun können wir ein **Mesh** (von Mathematikern oft `Grid` genannt) initialisieren"

# ╔═╡ 9d7940d6-4905-11eb-3d1d-b36e18ec60d6
grid = create_cook_grid(320,160)

# ╔═╡ 9582f632-4e92-11eb-204e-71d11cc15f37
md"""
Damit wir später über alle **Elemente iterieren** (Elemente Schleife) können, initialisieren wir uns den `DofHandler`. Er ist dafür verantwortlich einen Elemente Iterator zur Verfügung zu stellen und degress of freedom (dof) zu händeln.
"""

# ╔═╡ 91a2e654-4905-11eb-2fa1-430af651e6de
begin
	dh = DofHandler(grid)
	#wir pushen die gesuchte Lösung mit Namen u und Dimension 2
	push!(dh, :u, 2)
	close!(dh)
end

# ╔═╡ cd0fcf3a-4e92-11eb-17ec-157c6c4da712
md"""
Nun können wir mithilfe eines kleinen Tricks uns veranschaulichen wie der Aufbau der Steifigkeitsmatrix ist. Die Funktion wird von JuAFEM bereit gestellt.

`create_sparsity_pattern` allokiert allen nötigen Speicher für die Steifigkeitsmatrix $\mathbf{K}$ und stellt schon einmal die grundlegende Struktur der Matrix bereit.
"""

# ╔═╡ b72722a8-4907-11eb-079f-119a41d0d91f
K_structure = create_sparsity_pattern(dh)

# ╔═╡ c805bdfa-4907-11eb-32a8-f5c93f962213
begin
	import UnicodePlots
	fill!(K_structure.nzval, 1.0)
	UnicodePlots.spy(K_structure; height = 15)
end

# ╔═╡ 5d2ae784-50ef-11eb-316f-256a96f3bea8
md"""
!!! note "Achtung!"
    Beobachten sie die Struktur der Matrix, wenn wir das Mesh feiner wählen
"""

# ╔═╡ 5f8f3e0e-4e93-11eb-216c-cfe9f4dd20cd
md"""
Neben dem `DofHandler` gibt es den `ConstraintHandler`, der sich um das Management der Randbedingungen kümmert.

Er wird über `ConstraintHandler(dh)` erstellt, wobei `dh` der `DofHandler` des Problems ist. Zu Randbedingungen gehören essentiell zwei Informationen

- der Rand (also Ortsinformation)
- die Werte am Rand

Über `getfaceset` können vom `grid` markierte Ränder rangezogen werden. Für die Dirichlet Randbedingung ist `"clamped"` die relevante Markierung.

Eine Dirichlet Randbedingung wird über folgende Syntax erstellt in JuAFEM

```julia
Dirichlet(:u, ∂Ω, valuefunction)
```
wobei `:u` der Name der gesuchten Lösung ist (wurde zuvor in `DofHandler` festgelegt), ∂Ω stellt den zu beschreibenden Rand dar und `valuefunction` ist irgendeine Funktion die den Input `x,t` auf einen Dirichlet Randwert abbildet. Im unten stehenden Beispiel wurde dafür eine anonyme Funktion verwendet `(x,t) -> 0` (gesprochen: x,t bilden auf die konstante Funktion 0 ab).

Danach muss die Randbedingung über `add!` noch hinzugefügt werden und über `update!` die Zeitvariable `t` gesetzt werden. Für statische Probleme wird diese immer 0.0 sein.
"""

# ╔═╡ 1821efe8-4908-11eb-14f6-3f24ca70f1d1
begin
	ch = ConstraintHandler(dh)
	∂Ω = getfaceset(grid, "clamped")
	dbc = Dirichlet(:u, ∂Ω, (x, t) -> [0,0],[1,2])
	add!(ch, dbc)
	close!(ch)
	update!(ch, 0.0)
end

# ╔═╡ 534ed35c-4e98-11eb-2631-d38c5677ae24
md"""
## Assemblierung

### Elastizitätstensor

Der Elastizitätstensor wird wie in den voran gegangenen Übungen erstellt
"""

# ╔═╡ 4ea13f36-490c-11eb-3c24-c378588dde20
function elasticity(E, ν; dim=2)
	λ = E*ν / ((1 + ν) * (1 - 2ν))
	μ = E / (2(1 + ν))
	δ(i,j) = i == j ? 1.0 : 0.0
	f = (i,j,k,l) -> λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
	return SymmetricTensor{4,dim}(f)
end

# ╔═╡ e690864c-4e98-11eb-1a9a-e5b290db6135
md"""
### Assemblierungsfunktion
Wir wollen das Gleichungssystem folgender schwachen Form aufstellen:

$$\int_{\mathcal{B}} \delta\varepsilon \colon \mathbb{C} \colon \varepsilon \ \text{d} V = \int_{\partial\mathcal{B}} b \cdot \delta u \ \text{d} A + \int_{\partial\mathcal{B}_t} t \cdot \delta u\ \text{d} A \qquad \forall \delta u$$

mit Dirichlet Bedingung

$$u(\mathbf{x}) = g(\mathbf{x}) \qquad \forall \mathbf{x} \in \partial \mathcal{B}_D$$

**Verständnisfrage: Woher bekommen wir die Schwache Form mithilfe welcher Schritte?**

Anstatt über das gesamte Gebiet $\mathcal{B}$ zu integrieren, teilen wir das Gebiet in finite Elemente auf und bilden die Summe über alle Elemente (Elementsteifigkeit bzw. Elementlastvektor)

$$\sum_{E\in\mathcal{B}} \int_{E} \delta\varepsilon \colon \mathbb{C} \colon \varepsilon \ \text{d} V = \sum_{E\in\mathcal{B}} \int_{\partial E} b \cdot \delta u \ \text{d} A + \sum_{E\in\mathcal{B}} \int_{\partial E_t} t \cdot \delta u\ \text{d} A$$

Nun setzen wir unsere Ansatzfunktionen rein

$$\mathbf{u} = \sum_i^n \mathbf{N}_i \mathbf{d}_i$$

$$\varepsilon = \text{grad} \ \mathbf{u} = \sum_i^n \text{grad}(\mathbf{N}_i) \mathbf{d}_i$$

Da wir Galerkin FEM betreiben, setzen wir Ansatzfunktionen gleich Testfunktionen, daher

$$\delta\varepsilon = \text{grad} \ \delta u = \sum_i^n \text{grad}(\mathbf{N}_i)$$

Also können wir schreiben

$$\sum_{E\in\mathcal{B}} \int_{E} \text{grad}(\mathbf{N}_i) \colon \mathbb{C} \colon \text{grad}(\mathbf{N}_j) \mathbf{d_j} \ \text{d} V = \sum_{E\in\mathcal{B}} \int_{\partial E} b \cdot \mathbf{N}_i \ \text{d} A + \sum_{E\in\mathcal{B}} \int_{\partial E_t} t \cdot \mathbf{N}_i\ \text{d} A$$

da $\forall \delta u$ genau $N_i$ entspricht und somit 

$$\underbrace{\sum_{E\in\mathcal{B}} \int_{E} \text{grad}(\mathbf{N}_i) \colon \mathbb{C} \colon \text{grad}(\mathbf{N}_j) \ \text{d} V}_{\mathbf{K}_{ij}} \mathbf{d_j} = \underbrace{\sum_{E\in\mathcal{B}} \int_{\partial E} b \cdot \mathbf{N}_i \ \text{d} A + \sum_{E\in\mathcal{B}} \int_{\partial E_t} t \cdot \mathbf{N}_i\ \text{d} A}_{\mathbf{P_i}}$$

$$\mathbf{K}_{ij} \mathbf{d}_j = \mathbf{P}_i$$


Die Assemblierungsfunktion ist essentiell in folgende Schritte aufgeteilt:

1. Einmalige Initialisierung von globaler Steifigkeitsmatrix sowie globalen Lastvektor. Aber auch Elementsteifigkeitsmatrix und -lastvektor werden nur einmal initialisert, dafür aber in jedem Element überschrieben!
2. Definierung von Konstanten, in diesem Falle Volumenkraft `b` und Elastizitätstensor $ℂ$
3. **Schleife** über alle Elemente
    - Für jedes Element befüllen wir Elementsteifigkeitsmatrix sowie -lastvektor mit 0.0 `fill(Ke,0)`
    - Aktualisieren der aktuellen Elementwerte, die in `cell` gespeichert sind über `reinit!(cellvalues, cell)`
        - **Schleife** über alle Integrationspunkte `q_point in 1:getnquadpoints(cellvalues)`
            - **Schleife** über alle Ansatzfunktionen `i in 1:n_basefuncs`
            - Ansatzfunktionswert sowie den symmetrischen Gradienten an gegebener Stelle auswerten um die Lastvektor Anteile zu berechnen
                - **Schleife** über alle Ansatzfunktionen `j in 1:n_basefuncs`
                - Auswertung des symmetrischen Gradienten für Ansatzfunktion $j$ sowie Anteile auf Elementsteifigkeitsmatrix hinzu addieren
    - **Schleife** über alle `faces in 1:nfaces(cell)` 
        - **Schleife** über alle Quadraturpunkte `q_point` auf dem Rand
            - **Schleife** über alle Ansatzfunktionen `i` auf dem Rand
            - Anteile auf den Lastvektor hinzufügen
    - Elementkontributionen in das globale `K` und `f` assemblieren `assemble!(assembler, celldofs(cell), fe, Ke)`

"""

# ╔═╡ bd55c480-4908-11eb-36cc-0903839bc578
"""
	doassemble(elementvalues, dh) -> K, f
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

# ╔═╡ f95aa3f6-4912-11eb-38a6-b39e1aef5eb6
begin
	K, f = doassemble(elementvalues, randvalues, dh)
	apply!(K, f, ch)
	u = K \ f;
end

# ╔═╡ 739c6a38-519b-11eb-01f9-cfb2142ce91b
md"### Visualisierung der Verschiebung. Farbe entspricht der *vertikalen* Verschiebung."

# ╔═╡ 31a63ac0-503b-11eb-27c9-f1a0378daa8b
begin
	fig_u, ax_u, sp_u = JuAFEM.warp_by_vector(dh,u; colormap=:inferno, process=x->x[2])
	fig_u[1,2] = WGLMakie.Colorbar(fig_u; limits=(0,maximum(u)),colormap=:inferno, width=20)
	fig_u
end

# ╔═╡ a80c9800-4eb8-11eb-07b7-ebee935ed257
point = Vec(47.9,59.9)

# ╔═╡ ce2ed6f4-4eba-11eb-0a70-7f9670c20c03
measure_function(point,u,dh,ansatzfunktionen)

# ╔═╡ 8fb01998-50ee-11eb-344f-31d50aae98e3
md"""
## Aufgabe 

Gegeben ist das berechnete Verschiebungsfeld `u`, sowie unser gewählter Ansatz (lineare Dreieckslemente). Wie können Sie hieraus die Spannungen berechnen? Welche Schritte muss die Funktion `compute_stresses` durchführen?
"""

# ╔═╡ 2797b4fe-4f92-11eb-1cd3-5d87949ca4cd
σ = compute_stresses(elementvalues,dh,u,elasticity(210,0.3))

# ╔═╡ ae73a568-5010-11eb-0489-774cbf3956ce
projector = L2Projector(integration, ansatzfunktionen, grid);

# ╔═╡ c00cc0a2-5010-11eb-332f-3baf508ef8f8
σ_nodes = JuAFEM.project(σ, projector)

# ╔═╡ b977e69c-5199-11eb-30a3-e97c73476ff5
md" ### σ₁₁"

# ╔═╡ bb132868-503b-11eb-3070-bfe4ee2d9c40
begin
	plot_σ = [[σ_i[1],σ_i[4]] for σ_i in σ_nodes]
	plot_σ = vcat(plot_σ...)
	fig_σ11, ax_σ11, sp_σ11 = WGLMakie.mesh(dh,plot_σ, process=x->x[1],colormap=:inferno)
	fig_σ11[1,2] = WGLMakie.Colorbar(fig_σ11; limits=(minimum(plot_σ[1:2:end]),maximum(plot_σ[1:2:end])),colormap=:inferno, width=20)
	fig_σ11
end

# ╔═╡ 22e85d80-5197-11eb-123a-cb266d672b8b
md" ### σ₂₂"

# ╔═╡ 3e3ee84c-50ca-11eb-3185-1d1d0257c2b5
begin
	fig_σ22, ax_σ22, sp_σ22 = WGLMakie.mesh(dh, plot_σ,process=x->x[2],colormap=:inferno)
	fig_σ22[1,2] = WGLMakie.Colorbar(fig_σ22; limits=(minimum(plot_σ[2:2:end-1]),maximum(plot_σ[2:2:end-1])),colormap=:inferno, width=20)
	fig_σ22
end

# ╔═╡ 122af3e0-5197-11eb-20eb-1f8cd87988fd
md" ### σ₁₂"

# ╔═╡ 71d080a2-5053-11eb-01f5-01fe9b8a6772
begin
	plot_σ12 = [[σ_i[2],0] for σ_i in σ_nodes]
	plot_σ12 = vcat(plot_σ12...)
	fig_σ12, ax_σ12, sp_σ12 = WGLMakie.mesh(dh,plot_σ12, process=x->x[1],colormap=:inferno,scale_plot=false)
	fig_σ12[1,2] = WGLMakie.Colorbar(fig_σ12; limits=(minimum(plot_σ12[1:2:end]),maximum(plot_σ12[1:2:end])),colormap=:inferno, width=20)
	fig_σ12
end

# ╔═╡ 394fed0a-4f7a-11eb-3221-bb3bba5784d0
md"""
## Aufgabe

Führen sie eine Konvergenzstudie durch. Messen sie hierfür immer am Punkt A die y Komponente der Verschiebung. Überlegen sie sich hierfür eine geeignete Strategie die vorgestellten Schritte für die FE-Analyse mehrmals für die vorgegebenen `meshsizes` auszuführen. Die erste Spalte der Matrix korrepsondiert zu `nx`, die zweite zu `ny`.

Visualisieren sie Ihre Ergebnisse mit Plots!

"""

# ╔═╡ 9d2900a2-4f6e-11eb-0243-576fc03ecc4b
meshsizes = [10 5
			 20 10
			 40 20
			 80 40
		     160 80
			 320 160]

# ╔═╡ bbc5661e-519d-11eb-133c-857fa8735ce7
md"""
Zur Erinnerung Graphiken werden mit dem Paket Plots wie folgt erstellt 

```julia
plot(x,y)
```

Über die Zeilen einer Matrix kann man wie folgt iterieren
```julia
for size in eachrow(sizes)
...
end
```

"""

# ╔═╡ c6898046-51b7-11eb-1a8a-ebd4f9b71cc7
function convergence_study(sizes)
	u_study = [] #hier soll der Wert reingepusht werden über push!(u_study,u_measured)
	
	return u_study
end

# ╔═╡ e42fa276-51b7-11eb-1e3e-7bd6ad2854aa
study = convergence_study(meshsizes)

# ╔═╡ ca2ed45c-51f6-11eb-096f-cd20e0402bde
try
	if isapprox(study,[7.81095,8.38392,8.595,8.67185,8.69987,8.70929], atol=1e-4)
		correct()
	else
		keep_working()
	end
catch
	still_missing(md"Füllen sie die Funktion `convergence_study`")
end

# ╔═╡ 72b096f6-50ee-11eb-15a9-4519894c8957
plot(study, label = "u am Punkt A", legend=:bottomright)

# ╔═╡ Cell order:
# ╟─3940b9e4-4e75-11eb-29ec-4d1da14f089c
# ╠═fa20e096-48f1-11eb-1c40-f7c48ff8b2a4
# ╟─ee0b0c68-48ed-11eb-35e0-d513bd952297
# ╟─23b1189e-48ee-11eb-3652-37e022de561e
# ╟─1dc8fe14-48f1-11eb-37fb-ab620a99c16c
# ╟─6efe578a-4e98-11eb-33e3-d7edbb880dd5
# ╟─5ded3c3e-48f2-11eb-0883-5ddf347145ac
# ╠═2dfa3f64-48f2-11eb-1c16-0b6477748f43
# ╟─baf3f6b6-4e86-11eb-03a1-8d8d66bb944e
# ╠═262a7aee-48f7-11eb-0db3-3be29e5a9d6b
# ╟─a6be058e-4e88-11eb-0dfa-5110b45ac242
# ╟─c9abb7a8-4e8a-11eb-1ae5-11d947d575f4
# ╠═2e611812-48f7-11eb-29b9-c5a3e1651bfa
# ╟─62061634-4e8f-11eb-187c-177d497addd4
# ╠═39d110c4-48f7-11eb-3da9-37d44c9c73f3
# ╟─ca537870-4e90-11eb-0487-c7159de71eeb
# ╠═d9a7d7c0-490f-11eb-3c59-470d0a89297a
# ╠═04a20cfc-4910-11eb-15cd-31526230c2f6
# ╟─f3f77000-4e90-11eb-38d2-ff24eeabd75d
# ╠═9d7940d6-4905-11eb-3d1d-b36e18ec60d6
# ╟─9582f632-4e92-11eb-204e-71d11cc15f37
# ╠═91a2e654-4905-11eb-2fa1-430af651e6de
# ╟─cd0fcf3a-4e92-11eb-17ec-157c6c4da712
# ╠═b72722a8-4907-11eb-079f-119a41d0d91f
# ╟─c805bdfa-4907-11eb-32a8-f5c93f962213
# ╟─5d2ae784-50ef-11eb-316f-256a96f3bea8
# ╟─5f8f3e0e-4e93-11eb-216c-cfe9f4dd20cd
# ╠═1821efe8-4908-11eb-14f6-3f24ca70f1d1
# ╟─534ed35c-4e98-11eb-2631-d38c5677ae24
# ╠═4ea13f36-490c-11eb-3c24-c378588dde20
# ╟─e690864c-4e98-11eb-1a9a-e5b290db6135
# ╠═bd55c480-4908-11eb-36cc-0903839bc578
# ╠═f95aa3f6-4912-11eb-38a6-b39e1aef5eb6
# ╟─739c6a38-519b-11eb-01f9-cfb2142ce91b
# ╠═31a63ac0-503b-11eb-27c9-f1a0378daa8b
# ╠═a80c9800-4eb8-11eb-07b7-ebee935ed257
# ╠═ce2ed6f4-4eba-11eb-0a70-7f9670c20c03
# ╟─8fb01998-50ee-11eb-344f-31d50aae98e3
# ╠═2797b4fe-4f92-11eb-1cd3-5d87949ca4cd
# ╠═ae73a568-5010-11eb-0489-774cbf3956ce
# ╠═c00cc0a2-5010-11eb-332f-3baf508ef8f8
# ╟─b977e69c-5199-11eb-30a3-e97c73476ff5
# ╠═bb132868-503b-11eb-3070-bfe4ee2d9c40
# ╟─22e85d80-5197-11eb-123a-cb266d672b8b
# ╠═3e3ee84c-50ca-11eb-3185-1d1d0257c2b5
# ╟─122af3e0-5197-11eb-20eb-1f8cd87988fd
# ╠═71d080a2-5053-11eb-01f5-01fe9b8a6772
# ╟─394fed0a-4f7a-11eb-3221-bb3bba5784d0
# ╠═9d2900a2-4f6e-11eb-0243-576fc03ecc4b
# ╠═fce316a4-50ed-11eb-335d-7d24c4032a43
# ╟─bbc5661e-519d-11eb-133c-857fa8735ce7
# ╠═c6898046-51b7-11eb-1a8a-ebd4f9b71cc7
# ╠═e42fa276-51b7-11eb-1e3e-7bd6ad2854aa
# ╟─ca2ed45c-51f6-11eb-096f-cd20e0402bde
# ╠═72b096f6-50ee-11eb-15a9-4519894c8957
