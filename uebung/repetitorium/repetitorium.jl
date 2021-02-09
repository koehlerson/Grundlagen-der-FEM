### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# ╔═╡ ec6a786c-6b18-11eb-2dcb-07940b60575f
using JuAFEM

# ╔═╡ 061082f4-69ee-11eb-022d-938a3459b43c
md"""
# Repetitorium
"""

# ╔═╡ 8f064afa-69f6-11eb-31b6-47f70024daa7
md"""
## Aufgabe 1: Grundidee der FEM und Elementtypen
"""

# ╔═╡ 87cda91c-69f7-11eb-29f8-99958c09df92
md"""
Gegeben sei das elastische Potential

$\Pi = \frac{1}{2} \int_{\mathcal{B}}\varepsilon \cdot \mathbf{\sigma} \ \text{d}V - \int_{\mathcal{B}} \mathbf{u} \cdot \mathbf{b} \ \text{d}V - \int_{\partial\mathcal{B}_t} \mathbf{u} \cdot \mathbf{t} \ \text{d}A$

a) Wie lautet die zugehörige Variationsgleichung (schwache Form)?

b) Leiten sie aus der schachen Form die Euler Lagrange-Gleichungen (starke Form) für das Randwertproblem der Elastostatik her

c) Wie gelangt man von der schwachen Form zur Elementsteifigkeitsmatrix $\mathbf{k}_e$ und Elementlastmatrix $\mathbf{p}_e$?

d) Wie lassen sich die Integrale in der diskretisierten schwachen Form berechnen?

e) Wie gelangt man zur Lösung des globalen Problems?
"""

# ╔═╡ 3c920704-6afb-11eb-0dc7-49ee42e041ec
md"""
### Lösung 1a
$$\delta_u \Pi = \int_{\mathcal{B}} \delta \varepsilon \cdot \sigma \ \text{d}V - \int_{\mathcal{B}} \delta \mathbf{u} \cdot \mathbf{b} \ \text{d}V - \int_{\partial\mathcal{B}_t} \delta \mathbf{u} \cdot \mathbf{t} \ \text{d}A = 0$$

Herleitung in Übung 3 mit allgemeiner Formel die erste Variation ($\delta_u \Pi$) eines Funktionals ($\Pi$) zu bilden
"""

# ╔═╡ eb963a7a-6afd-11eb-0414-dd22488694d4
md"""
### Lösung 1b
Wir wissen dass $\sigma$ symmetrisch ist, daher gilt

$\delta\varepsilon \cdot \sigma = \nabla \delta \mathbf{u} \cdot \sigma$

Kann in Indexnotation bewiesen werden und wurde in der Vorlesung hergeleitet.
Einsetzen in die Schwache Form:

$\delta_u \Pi = \int_{\mathcal{B}} \nabla \delta \mathbf{u} \cdot \sigma \ \text{d}V - \int_{\mathcal{B}} \delta \mathbf{u} \cdot \mathbf{b} \ \text{d}V - \int_{\partial\mathcal{B}_t} \delta \mathbf{u} \cdot \mathbf{t} \ \text{d}A = 0.$

Notwendigen Schritte um nun auf die Euler-Lagrange Gleichung (starke Form) zu kommen

!!! note "Erforderliche Schritte"
    I) Partielle Integration

    $\int_{\mathcal{B}} \nabla \delta \mathbf{u} \cdot \sigma \ \text{d}V = \int_{\mathcal{B}} \text{div}(\delta \mathbf{u} \cdot \sigma) \ \text{d}V - \int_{\mathcal{B}} \delta \mathbf{u} \cdot \text{div}(\sigma) \ \text{d}V$

    II) Gauss'scher Integralsatz

    $\int_{\mathcal{B}} \text{div}(\delta \mathbf{u} \cdot \sigma) \ \text{d}V = \int_{\partial\mathcal{B}} \delta \mathbf{u} \cdot \sigma \mathbf{n} \ \text{d}A$
    
    III) Dirichlet/Neumann-Rand

    $\int_{\partial\mathcal{B}} \delta \mathbf{u} \cdot \sigma \mathbf{n} \ \text{d}A = \int_{\partial\mathcal{B_t}} \delta \mathbf{u} \cdot \sigma \mathbf{n} \ \text{d}A + \underbrace{\int_{\partial\mathcal{B_u}} \delta \mathbf{u} \cdot \sigma \mathbf{n} \ \text{d}A}_{ = 0, \text{ da } \delta \mathbf{u}|_{\partial\mathcal{B}_u}=0}$

$\implies \delta_u \Pi = \int_{\mathcal{B}} \nabla \delta \mathbf{u} \cdot [-\text{div}(\sigma - \mathbf{b})]\ \text{d}V + \int_{\partial\mathcal{B}_t} \delta \mathbf{u} \cdot [\sigma \mathbf{n} - \mathbf{t}] \ \text{d}A = 0$

!!! note "Randwertproblem der Elastostatik (starke Form)"
    Finde $\mathbf{u}(\mathbf{x}) \forall \mathbf{x} \in \mathcal{B}$, so dass
    
    $-\text{div} \sigma = \mathbf{b} \quad \text{in } \mathcal{B}$

    $\sigma \mathbf{n} = \mathbf{t} \quad \text{auf } \partial \mathcal{B}_t$

    $\mathbf{u} = \mathbf{u}_D \quad \text{auf } \partial \mathcal{B}_u$ 
"""

# ╔═╡ ee4975de-6afd-11eb-1d47-fdffdd4f7417
md"""
### Lösung 1c
Anstelle die Integrale direkt über das Gebiet $\mathcal{B}$ aufzulösen, teilen wir das Gebiet in eine endliche (finite) Anzahl von Elementen. Nun lösen wir die Integrale auf den Elementen und summieren diese auf.

$$\delta_u \Pi = \int_{\mathcal{B}} \delta \varepsilon \cdot \sigma \ \text{d}V - \int_{\mathcal{B}} \delta \mathbf{u} \cdot \mathbf{b} \ \text{d}V - \int_{\partial\mathcal{B}_t} \delta \mathbf{u} \cdot \mathbf{t} \ \text{d}A = 0$$

$$\implies$$

$$\sum_e \int_{\mathcal{B}_e} \delta \varepsilon \cdot \sigma \ \text{d}V - \int_{\mathcal{B}_e} \delta \mathbf{u} \cdot \mathbf{b} \ \text{d}V - \int_{\partial\mathcal{B}_{e_t}} \delta \mathbf{u} \cdot \mathbf{t} \ \text{d}A = 0$$

Nun setzen wir unsere Ansatzfunktionen ein für die unbekannte Lösung $\mathbf{u}$ und den Testfunktionen $\delta\mathbf{u}$

$$\sum_e \int_{\mathcal{B}_e} \delta \varepsilon^h \cdot \sigma^h \ \text{d}V - \int_{\mathcal{B}_e} \delta \mathbf{u}^h \cdot \mathbf{b} \ \text{d}V - \int_{\partial\mathcal{B}_{e_t}} \delta \mathbf{u}^h  \cdot \mathbf{t} \ \text{d}A = 0$$

Approximation der relevanten Lösungsgrößen

$\delta \varepsilon^h = \mathbf{B}_e \delta \mathbf{d}_e$
$\sigma^h = \mathbb{C} \mathbf{B}_e \mathbf{d}_e$
$\delta \mathbf{u}^h = \mathbf{N}_e \delta \mathbf{d}_e$

$$\mathbf{B}^e =
\begin{bmatrix}
N_{1,x} & 0 &  &  & N_{k,x} & 0 \\
0 & N_{1,y} & & ... & 0& N_{k,y} \\
N_{1,y} & N_{1,x} & &  & N_{k,y} & N_{k,x}
\end{bmatrix}$$

$$\mathbf{N}^e =
\begin{bmatrix}
N_{1} & 0 &  &...  & N_{k} & 0 \\
0 & N_{1} & & ... & 0 & N_{k} \\
\end{bmatrix}$$

In Detail in Übung 7, 9.

Elementsteifigkeitsmatrix:

$\mathbf{k}_e = \int_{\mathcal{B}_e}\mathbf{B}_e^T \mathbb{C} \mathbf{B}_e \ \text{d}V$

Elementlastvektor:

$\mathbf{p}_e = \int_{\mathcal{B}_e} \mathbf{N}_e^T \mathbf{b} \ \text{d}V + \int_{\mathcal{B}_{e_t}} \mathbf{N}_e^T  \cdot \mathbf{t} \ \text{d}A$

Eine andere Perspektive um auf das diskrete System von Steifigkeitsmatrix und Lastvektor zu kommen wurde in Übung 6 besprochen.
"""

# ╔═╡ 31a0e074-6b0d-11eb-0f1c-d7cf66382b7c
md"""
### Lösung 1d

**1) Verwenden von Referenzelement im Parameterraum**
"""

# ╔═╡ 29f25f92-6b0d-11eb-340e-17b92cebcecf
HTML(open(f->read(f, String), "assets/quadrilateral.svg"))

# ╔═╡ f45fac36-6afd-11eb-1cd3-d9a4f9d28a74
md"""
- Ansatzfunktionen für alle Elemente gleich
- numerische Integration möglich

$\mathbf{N}_e := \mathbf{N}_e(\xi,\eta)$
$\mathbf{B}_e := \mathbf{B}_e(\xi,\eta)$ (bei linearen Ansatzfunktionen ist $\mathbf{B}_e$ konstant)

**2) Isoparametrisches Konzept**

Gleiche Ansatzfunktionen für $\mathbf{x}^h$ und $\mathbf{u}^h$. Geometrische Ansatzfunktionen werden hier zur Unterscheidung $\mathbf{M}$ genannt, auch wenn sie die gleichen wie $\mathbf{N}$ sind.

$$\mathbf{x}^h = \mathbf{N}_e \mathbf{x}_e$$
$$\mathbf{u}^h = \mathbf{N}_e \mathbf{d}_e$$

- Jacobimatrix $\mathbf{J} = \dfrac{\partial \mathbf{x}^h}{\partial \mathbf{\xi}}$
- Jacobi Determinante $\det \mathbf{J}$

$$\implies \mathbf{k}_e = \int_{-1}^1 \int_{-1}^1 \mathbf{B}_e^T(\xi, \eta) \mathbb{C} \mathbf{B}_e(\xi, \eta) \det \mathbf{J}(\xi,\eta) \ \text{d}\xi\text{d}\eta$$
$$\implies \mathbf{p}_e = \int_{-1}^1 \int_{-1}^1 \mathbf{N}_e^T(\xi, \eta) \mathbf{b} \det \mathbf{J}(\xi,\eta) \ \text{d}\xi\text{d}\eta$$

**3) Numerische Integration**

- $\mathbf{k}_e = \sum_{i=1}^{n_{GP}}\mathbf{B}_e^T(\xi\_i, \eta\_i) \mathbb{C} \mathbf{B}_e(\xi_i, \eta_i) \det \mathbf{J}(\xi_i,\eta_i) \ w_{\xi_i} w_{\eta_i}$

- analog $\mathbf{p}_e$
"""

# ╔═╡ f7a3b1bc-6afd-11eb-3eef-77b4d527fe4c
md"""
### Lösung 1e
Assemblierung zu globalem Gleichungssystem

- globale Steifigkeitsmatrix $\mathbf{K}$
- Knotenverschiebungen $\mathbf{D}$
- globaler Lastvektor $\mathbf{P}$


$$\mathbf{K}\mathbf{D} = \mathbf{P}$$

Lösung

$$\mathbf{D} = \mathbf{K}^{-1}\mathbf{P}$$
"""

# ╔═╡ 77e8fa5a-69f8-11eb-3d0e-4b9ab8b89a9d
md"""
## Aufgabe 2: Assemblierung
Gegeben sei das in der unteren Abbildung dargestellte Mini-Randwertproblem bestehend aus linearen Dreieckselementen
"""

# ╔═╡ fc3a0c54-69fd-11eb-1581-57781dd5867e
HTML(open(f->read(f, String), "assets/kv_assem.svg"))

# ╔═╡ 398b8f42-69fe-11eb-093f-d96e8db372e7
md"""
a) Erstellen Sie eine Skizze in der alle globalen Knotenfreiheitsgrade (ohne Berücksichtigung der Randbedingungen eingetragen sind

b) Wie lauten die Einträge des globalen Lösungsvektors unter Berücksichtigung der Randbedingungen? Wie lautet der globale Lastvektor?

c) Wie lauten die Einträge der Elementsteifigkeitsmatrix $\mathbf{k}^1$ und $\mathbf{k}^2$, welche nach Berücksichtigung der Randbedingungen die globale Steifigkeitsmatrix eingehen? Bestimmen sie Einträge der globalen Steifigkeitsmatrix $\mathbf{K}$

d) Bilden sie die Jacobimatrix für Element 2. Nehmen sie dafür folgende Knotenkoordinaten an

| Knoten       |   x   | y |
| :------------- | :----------: | -----------: |
| Nummer 3    | 0 | 0 |
| Nummer 4    | 5 | 0 |
| Nummer 2    | 5 | 5 |
"""

# ╔═╡ 3cd7d928-6afb-11eb-2ee3-5b4da6501764
md"""
### Lösung 2a
"""

# ╔═╡ 8945eb08-6b17-11eb-1e65-3549cf210734
HTML(open(f->read(f, String), "assets/kv_assem_dof.svg"))

# ╔═╡ da65d508-6afd-11eb-31f7-2fc89137aad6
md"""
### Lösung 2b

$\mathbf{D} = \begin{bmatrix}
D_1 \\
D_3 \\
D_4 \\
D_7 
\end{bmatrix} \qquad
\mathbf{P} = \begin{bmatrix}
0 \\
-F \\
0 \\
-F
\end{bmatrix}$
"""

# ╔═╡ e06832b6-6afd-11eb-28a6-edc160ac2b4c
md"""
### Lösung 2c

$\mathbf{k}^1 = \begin{bmatrix}
k_{22}^1 & k_{23}^1 & k_{24}^1 \\
k_{32}^1 & k_{33}^1 & k_{34}^1 \\
k_{42}^1 & k_{43}^1 & k_{44}^1
\end{bmatrix} \qquad
\mathbf{k}^2 = \begin{bmatrix}
k_{33}^2 & k_{34}^2 & k_{37}^2 \\
k_{43}^2 & k_{44}^2 & k_{47}^2 \\
k_{73}^2 & k_{74}^2 & k_{77}^2
\end{bmatrix}$

$\mathbf{K} = \begin{bmatrix}
k_{22}^1 & k_{23}^1 & k_{24}^1 & 0 \\
k_{32}^2 & k_{33}^1+k_{33}^2 & k_{34}^1+k_{34}^2 & k_{37}^2 \\
k_{42}^1 & k_{43}^1+k_{43}^2 & k_{44}^1+k_{44}^2 & k_{47}^2 \\
0 & k_{73}^2 & k_{74}^2 & k_{77}^2
\end{bmatrix}$

"""

# ╔═╡ e3f8a578-6afd-11eb-0b97-5961a1564ef3
md"""
### Lösung 2d
Jacobi Berechnungen besprochen in Übung 9
"""

# ╔═╡ e20fee72-6b18-11eb-04c8-1f286580c5aa
begin
	n1 = Node((0.0,0.0))
	n2 = Node((5.0, 0.0))
	n3 = Node((5.0,5.0))
	element = Triangle((1,2,3))
	one_element_grid = Grid([element],[n1,n2,n3])
	cv = CellVectorValues(QuadratureRule{2, RefTetrahedron}(1), Lagrange{2, RefTetrahedron, 1}())
	n_geom_basefuncs = JuAFEM.getngeobasefunctions(cv) # 3 Stück, jeweils eine pro Knoten
	
	J = zero(Tensor{2,2})
	x = getcoordinates(one_element_grid,1)
    for j in 1:n_geom_basefuncs #Loop über alle geometrischen Ansatzfunktionen
        global J += x[j] ⊗ cv.dMdξ[j]
    end
	detJ = det(J)
    Jinv = inv(J)
	J
end

# ╔═╡ fe175474-6b19-11eb-19f1-815b8090ca58
cv.dMdξ

# ╔═╡ c3748d46-6b19-11eb-05e5-5ba0e4f3e793
J

# ╔═╡ e9cc575c-6b18-11eb-2bda-1f0c754e2ac4
det(J)/2

# ╔═╡ 25d035ee-6b1a-11eb-3e75-f95fca75d23c
x[1] ⊗ cv.dMdξ[1]

# ╔═╡ 5cf5c9ee-6b1a-11eb-0316-0557ba410ab7
x[2] ⊗ cv.dMdξ[2]

# ╔═╡ 6ce1ebbe-6b1a-11eb-161f-693555832e4d
x[3] ⊗ cv.dMdξ[3]

# ╔═╡ Cell order:
# ╟─061082f4-69ee-11eb-022d-938a3459b43c
# ╟─8f064afa-69f6-11eb-31b6-47f70024daa7
# ╟─87cda91c-69f7-11eb-29f8-99958c09df92
# ╟─3c920704-6afb-11eb-0dc7-49ee42e041ec
# ╟─eb963a7a-6afd-11eb-0414-dd22488694d4
# ╟─ee4975de-6afd-11eb-1d47-fdffdd4f7417
# ╟─31a0e074-6b0d-11eb-0f1c-d7cf66382b7c
# ╟─29f25f92-6b0d-11eb-340e-17b92cebcecf
# ╟─f45fac36-6afd-11eb-1cd3-d9a4f9d28a74
# ╟─f7a3b1bc-6afd-11eb-3eef-77b4d527fe4c
# ╟─77e8fa5a-69f8-11eb-3d0e-4b9ab8b89a9d
# ╟─fc3a0c54-69fd-11eb-1581-57781dd5867e
# ╟─398b8f42-69fe-11eb-093f-d96e8db372e7
# ╟─3cd7d928-6afb-11eb-2ee3-5b4da6501764
# ╟─8945eb08-6b17-11eb-1e65-3549cf210734
# ╟─da65d508-6afd-11eb-31f7-2fc89137aad6
# ╟─e06832b6-6afd-11eb-28a6-edc160ac2b4c
# ╠═ec6a786c-6b18-11eb-2dcb-07940b60575f
# ╟─e3f8a578-6afd-11eb-0b97-5961a1564ef3
# ╠═e20fee72-6b18-11eb-04c8-1f286580c5aa
# ╠═fe175474-6b19-11eb-19f1-815b8090ca58
# ╠═c3748d46-6b19-11eb-05e5-5ba0e4f3e793
# ╠═e9cc575c-6b18-11eb-2bda-1f0c754e2ac4
# ╠═25d035ee-6b1a-11eb-3e75-f95fca75d23c
# ╠═5cf5c9ee-6b1a-11eb-0316-0557ba410ab7
# ╠═6ce1ebbe-6b1a-11eb-161f-693555832e4d
