### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# ╔═╡ 47277816-518d-11eb-244f-fd99d21c2a6e
include("../../definitions/def.jl")

# ╔═╡ 96d1e7b6-596c-11eb-1305-716af2ea8390
import Plots

# ╔═╡ 1d876196-518c-11eb-0144-d3ec17e099b4
md""" 
# Ansatzfunktionen
Willkommen in der achten Übung im Modul Grundlagen der FEM. In dieser Übung werden wir uns mit Ansatzfunktionen beschäftigen.
"""

# ╔═╡ f78f5792-518c-11eb-1f4d-5bdf07c2ae20
md"""
# Aufgabe 8.1 - Verständnisfragen
"""

# ╔═╡ 8fe076c0-518d-11eb-10ae-2f725b2f0a05
md"""
## Aufgabe 8.1a) 
Wie lassen sich grundsätzlich FE-Ansatzfunktionen bestimmen?
"""

# ╔═╡ 50078ece-59a0-11eb-2945-4376d914232c
higherorder1d=HTML(open(f->read(f, String), "assets/5thorder.svg"));

# ╔═╡ 86dc8bde-59a2-11eb-1e90-bd4a56e53b43
datatable =
md"""
| Eigenschaft   |   Wert     |
|---------------|------------|
|Polynomgrad p  | 5          |
|Knoten k       | 1          |
|sonstige Knoten i| $i\in\{0,2,3,4,5\}$ |
""";

# ╔═╡ dccd7022-518e-11eb-2d7d-990c9ad72a32
md"""
## Aufgabe 8.1b) 
Welche Eigenschaft der Ansatzfunktionen wird dabei genutzt?
"""

# ╔═╡ f952cc6a-518e-11eb-3cfa-2be6123c5271
solution(
md"""
$N_I(\boldsymbol{x}_J) = \delta_{IJ}$
""", blur=false)

# ╔═╡ 6c07fbe0-518f-11eb-351f-1b138a4d859c
md"""
## Aufgabe 8.1c) 
Was ist die Idee des isoparametrischen Konzeptes?
"""

# ╔═╡ 44703a2e-59b5-11eb-0176-61aee94d5348
firstorderiso=HTML(open(f->read(f, String), "assets/firstorder1d-iso.svg"));

# ╔═╡ 833ffd44-518f-11eb-261f-359001e33f7d
solution(
md"""
- Wir interpolieren die Geometrie mit den gleichen Ansatzfunktionen wie unsere Lösung (Verschiebung)
- Mithilfe dieser Idee kann jeder Punkt innerhalb eines verzerrten Elements durch die Knotenkoordinaten des Elements im physikalischen Raum und die Knotenpunktansätze im Einheits- oder isoparametrischen Raum beschrieben werden.
	
$firstorderiso
	
""", blur=false)

# ╔═╡ 95d3074e-518f-11eb-09a7-1b8cf299caf0
md"""
## Aufgabe 8.1d) 
Welche Vorteile hat es die FE-Approximation über das Referenzelement im Parameterraum durchuführen?
"""

# ╔═╡ af4b8480-518f-11eb-2b47-334c93cf76ad
solution(
md"""
a) Ansatzfunktion wird einmal für das Referenzelement aufgestellt und kann für jedes Element verwendet werden 
	
b) Weiterer Vorteil: numerische Integration möglich
	
Beispiel: lineares Element 
	
$\boxed{1} \quad k = 0 \quad \xi_k = \xi_1 \quad i \in \{1\} \quad \xi_i \in \{x_2\} \Rightarrow N_1 = l_0^1(\xi) = \frac{\xi - 1}{-1-1} = \frac{1}{2}(1-\xi)$
	
$\boxed{2} \quad k = 1 \quad \xi_k = \xi_2 \quad i =\{0\} \quad \xi_i \in \{\xi_1\} \Rightarrow N_2 = l_1^1(\xi) = \frac{\xi + 1}{1+1} = \frac{1}{2}(1+\xi)$
""", blur=false)

# ╔═╡ d1784dc6-519f-11eb-0d92-43feaf202c19
md"""
# Aufgabe 8.2 - Isoparametrisches Viereckselement
"""

# ╔═╡ a5eae37e-51a1-11eb-248e-0ff3ac1a4218
HTML(open(f->read(f, String), "assets/8.2.svg"))

# ╔═╡ c62c5492-51a1-11eb-259c-f5cb7cbf35cf
md"""
## Aufgabe 8.2a) 
 Bestimmen Sie die Ansatzfunktionen $N_1$ und $N_2$ für das lineare Viereckselement aus der Abbildung über diesen Text
"""

# ╔═╡ 65c1c186-51a2-11eb-224a-b3c2e521b61b
hint(
md"""
Benutzen Sie Lagrange-Polynome für beide Koordinaten
	
$$l_k^p(\xi, \eta) = l_k^p(\xi)\ l_k^p(\eta) = 
	\prod_{i=0, i\neq k}^p \dfrac{(\xi-\xi_i)}{(\xi_k-\xi_i)}\prod_{i=0, i\neq k}^p \dfrac{(\eta-\eta_i)}{(\eta_k-\eta_i)}$$
	
Betrachten sie dazu die Linien "separat" als würden sie jeweils eine 1D Funktion aufstellen und anschließend miteinander multiplizieren
"""; blur=false)

# ╔═╡ 2976907e-51a2-11eb-03c9-5597c9dff33d
function N₁l(ξ,η)
    return 1/4*(1-η)*(1-ξ)
end

# ╔═╡ facff0ca-51a4-11eb-0918-99ec861bfde4
begin
	ξ1 = collect(-1:0.1:1)
	η1 = collect(-1.:0.1:1.)
	Plots.surface(ξ1,η1,N₁l,camera = (30, 30))
end

# ╔═╡ 57fa8e8a-59a8-11eb-31b2-db3166b9ebe2
firstorder =HTML(open(f->read(f, String), "assets/firstorder1d.svg"));

# ╔═╡ b6d4fe8a-518e-11eb-1783-179c0af183d0
solution(
md"""
$\boldsymbol{u}(\boldsymbol{x}) = N_I(\boldsymbol{x})\boldsymbol{d}_I$
	
Minimialbeispiel: lineares 1D - Stabelement
	
|     Fig.8.1.1 |            |
|---------------|------------|
|$(firstorder)  |linear Ansatz $N_I = a_I + b_Ix$|
	
Bestimmen der Koeffizienten:
	
$I = 1\Rightarrow \begin{gather}
	a_1 + b_1x_1 = 1 \\
	a_1 + b_1x_2 = 0
	\end{gather} \Leftrightarrow  \underbrace{\begin{bmatrix}
	1 & 0 \\
	1 & l
	\end{bmatrix}}_{\boldsymbol{C}} \begin{bmatrix}
	a_1\\b_1 \end{bmatrix} = \begin{bmatrix}
	1\\0 \end{bmatrix}$
	
$\boldsymbol{C}^{-1} = \frac{1}{l} \begin{bmatrix}
	l & 0 \\
	-1 & 1
	\end{bmatrix}$
	
$\begin{bmatrix}
	a_1\\b_1 \end{bmatrix} = \begin{bmatrix}
	1 & 0 \\
	-\frac{1}{l} & \frac{1}{l}
	\end{bmatrix} \begin{bmatrix}
	1\\0 \end{bmatrix} = \begin{bmatrix}
	1\\-\frac{1}{l} \end{bmatrix} \Rightarrow N_1 = 1 - \frac{1}{l}x$
	
$I = 2\Rightarrow  \begin{bmatrix}
	a_2\\b_2 \end{bmatrix} = \begin{bmatrix}
	1 & 0 \\
	-\frac{1}{l} & \frac{1}{l}
	\end{bmatrix} \begin{bmatrix}
	0\\1 \end{bmatrix} = \begin{bmatrix}
	0\\\frac{1}{l} \end{bmatrix} \Rightarrow N_2 = \frac{1}{l}x$
	
	
Methode kann für jede Art von Element verwendet werden.
	
**Probleme der Methode**
	
  	1. bei größeren Elementen hoher Rechenaufwand zur Berechnung von $\boldsymbol{C}^{-1}$
	
 	2. muss für jedes Element neu aufgestellt werden
	
Problem 1 kann mit sogenannten Lagrange-Polynomen umgangen werden
	
$l_k^p(x) =
	\prod_{i=0, i\neq k}^p \dfrac{(x-x_i)}{(x_k-x_i)}$
	

|     Fig.8.1.2   |            |
|---------------|------------|
|$(higherorder1d) |$(datatable)|


	
$k = 1 \quad i \in \{0,2,3,4,5 \} \quad x_k = x_2 \quad x_i \in \{x_1,x_3,x_4,x_5,x_6 \}$
	
Beispiel: p = 1 
	
$k = 0 \quad x_k = x_1 \quad i \in \{1\} \quad x_i \in \{x_2\} \Rightarrow N_1 = l_0^1(x) = \frac{x-l}{l-0} =1 - \frac{x}{l}$
	
$k = 1 \quad x_k = x_2 \quad i = \{0\} \quad x_i \in \{x_1\} \Rightarrow N_2 = l_1^1(x) = \frac{x-0}{l-0} = \frac{x}{l}$	
""", blur=false)

# ╔═╡ 1e393458-51a6-11eb-3e4f-1dc7b554177b
function N₂l(ξ,η)
    return 1/4*(1-η)*(1+ξ)
end

# ╔═╡ 418267fe-51a6-11eb-2afd-3748ac675115
begin
	ξ2 = collect(-1:0.1:1)
	η2 = collect(-1.:0.1:1.)
	Plots.surface(ξ2,η2,N₂l,camera = (20, 40))
end

# ╔═╡ 9d9c4c0c-51a2-11eb-0350-83efe127a32f
solution(md"""
$\boxed{1} \quad \xi: \quad k = 0, \quad i \in \{1\}, \quad \xi_k = -1, \quad \xi_i \in \{1\}$
	
$\qquad \eta: \quad k = 0, \quad i \in \{1\}, \quad \eta_k = -1, \quad \eta_i \in \{1\}$
	
$\Rightarrow N_1 = l_0^1(\xi)l_0^1(\eta) = \frac{(\xi - 1)(\eta - 1)}{(-1-1)(-1-1)} = \frac{1}{4}(1-\xi)(1-\eta)$
	
$\boxed{2} \quad \xi: \quad k = 1, \quad i \in \{0\}, \quad \xi_k = 1, \quad \xi \in \{-1\}$
	
$\qquad \eta: \quad k = 0, \quad i \in \{1\}, \quad \eta_k = -1, \quad \eta_i \in \{1\}$
	
$\Rightarrow N_2 = l_1^1(\xi)l_0^1(\eta) = \frac{(\xi + 1)(\eta - 1)}{(1+1)(-1-1)} = \frac{1}{4}(1+\xi)(1-\eta)$	
""",blur=false)

# ╔═╡ 99a2d9c8-51a6-11eb-2e69-2b5ca5565d69
md"""
## Aufgabe 8.2b) 
 Bestimmen Sie die Ansatzfunktionen $N_2$ und $N_6$ für das quadratische Viereckselement aus der Abbildung
"""

# ╔═╡ f6ae9760-51c4-11eb-3861-87e7030628a8
function N₂q(ξ,η)
    return 1/4*ξ*(1+ξ)*η*(-1+η)
end

# ╔═╡ 57a7f00c-51c5-11eb-216e-1f2c676d81ea
begin
	ξ3 = collect(-1:0.1:1)
	η3 = collect(-1.:0.1:1.)
	Plots.surface(ξ3,η3,N₂q,camera = (30, 40))
end

# ╔═╡ a50129e0-51c5-11eb-3126-6345cf39d287
function N₆q(ξ,η)
    return -1/2*ξ*(1+ξ)*(-1+η^2)
end

# ╔═╡ ca71f52e-51c5-11eb-0597-097eaed7b7b2
begin
	ξ4 = collect(-1:0.1:1)
	η4 = collect(-1.:0.1:1.)
	Plots.surface(ξ4,η4,N₆q,camera = (40, 40))
end

# ╔═╡ 53c5d182-51a9-11eb-3ce0-59bc172b982d
solution(md"""
$\boxed{2} \quad \xi: \quad k = 2, \quad i \in \{0,1\}, \quad \xi_k = 1, \quad \xi_i \in \{-1,0\}$
	
$\qquad \eta: \quad k = 0, \quad i \in \{1,2\}, \quad \eta_k = -1, \quad \eta_i \in \{0,1\}$
	
$\Rightarrow N_2 = l_2^2(\xi)l_0^2(\eta) = \frac{(\xi + 1)(\xi-0)(\eta-0)(\eta - 1)}{(1+1)(1-0)(-1-0)(-1-1)} = \frac{1}{4}\xi(1+\xi)\eta(-1+\eta)$

$\boxed{6} \quad \xi: \quad k = 2, \quad i \in \{0,1\}, \quad \xi_k = 1, \quad \xi_i \in \{-1,0\}$
	
$\qquad \eta: \quad k = 1, \quad i \in \{0,2\}, \quad \eta_k = 0, \quad \eta_i \in \{-1,1\}$
	
$\Rightarrow N_6 = l_2^2(\xi)l_1^2(\eta) = \frac{(\xi + 1)(\xi-0)(\eta+1)(\eta - 1)}{(1+1)(1-0)(0+1)(0-1)} = -\frac{1}{2}\xi(1+\xi)(-1+\eta^2)$
""",blur=false)

# ╔═╡ 5b81e7e8-51aa-11eb-1af1-8b08c6d9b441
md"""
# Aufgabe 8.3 - Referenzdreieckselement
"""

# ╔═╡ 0c007e40-51ab-11eb-03c2-afb37f8f292e
HTML(open(f->read(f, String), "assets/8.3.svg"))

# ╔═╡ cd0c8f48-51ab-11eb-12f1-09e5c6402e21
md"""
## Aufgabe 8.3a) 
 Bestimmen Sie die Ansatzfunktionen für das lineare Dreieckselement aus der Abbildung mithilfe des allgemeinen Ansatzes
"""

# ╔═╡ 138b87f4-51ba-11eb-1933-7901c56f9b24
hint(md"""
allgemeine Ansatz 
$$N_I = a_I + b_I\xi + c_I\eta$$
	
Spalten stehen für $1, \xi_J, \eta_J$
	
$\begin{bmatrix}
	1 & 1 & 0 \\
	1 & 0 & 1 \\
	1 & 0 & 0 \end{bmatrix} \begin{bmatrix}
	a_I \\ b_I \\ c_I \end{bmatrix} = \delta_{IJ}$
""", blur=false)

# ╔═╡ 61a0467e-5669-11eb-1387-19be37b3cc37
md"""
### N₁ lineares Dreieck:
"""

# ╔═╡ a4b48c88-5669-11eb-0dc4-1506e071610a
[1 1 0
 1 0 1
 1 0 0]\[1,0,0]	

# ╔═╡ 2d3c53fa-51c8-11eb-064d-cde446cc7d72
function N₁tl(ξ,η)
    return ξ
end

# ╔═╡ 9b77a9c2-51c9-11eb-1d15-1f104be60183
begin
	ξ5 = collect(0:0.02:1)
	η5 = collect(0.:0.02:1)
	Plots.surface(ξ5,η5, (x1,x2) -> (1 - x1 - x2) > 0.0 ? N₁tl(x1,x2) : 0.0,camera = (10, 40)) 
end

# ╔═╡ 3db4ba64-566a-11eb-2230-bbe82a13b8ca
md"""
### N₂ lineares Dreieck:
"""

# ╔═╡ 4654300a-566a-11eb-313f-03cdb24885fa
[1 1 0
 1 0 1
 1 0 0]\[0,1,0]	

# ╔═╡ b5b78574-51c8-11eb-0958-b53010124aa1
function N₂tl(ξ,η)
    return η
end

# ╔═╡ 58147a8e-566a-11eb-3311-dd04ee91c1f7
begin
	ξ6 = collect(0:0.02:1)
	η6 = collect(0.:0.02:1)
	Plots.surface(ξ6,η6, (x1,x2) -> (1 - x1 - x2) > 0.0 ? N₂tl(x1,x2) : 0.0,camera = (60, 50)) 
end

# ╔═╡ a2d0df5e-566a-11eb-1c58-d5631e9c7635
md"""
### N₃ lineares Dreieck:
"""

# ╔═╡ b07170ba-566a-11eb-3244-0d0fd670f8df
[1 1 0
 1 0 1
 1 0 0]\[0,0,1]	

# ╔═╡ c3c89a42-51c8-11eb-2b9d-b1bcdc3f6690
function N₃tl(ξ,η)
    return 1-ξ-η
end

# ╔═╡ bb652a84-566a-11eb-25f8-83c449659164
begin
	ξ7 = collect(0:0.05:1)
	η7 = collect(0.:0.05:1)
	Plots.surface(ξ7,η7, (x1,x2) -> (1 - x1 - x2) > 0.0 ? N₃tl(x1,x2) : 0.0,camera = (40, 50)) 
end

# ╔═╡ 18d09e7e-51ac-11eb-0dfc-5ffefbcafdc0
solution(md"""

	
Direktes einsetzen: 
	
$\boxed{1} \qquad \begin{bmatrix}
	1 & 1 & 0 \\
	1 & 0 & 1 \\
	1 & 0 & 0 \end{bmatrix} \begin{bmatrix}
	a_1 \\ b_1 \\ c_1 \end{bmatrix} = \begin{bmatrix}
	1 \\ 0 \\ 0 \end{bmatrix} \Rightarrow a_1 = 0, \quad c_1 = 0, \quad b_1 = 1 \Rightarrow N_1 = \xi$
	
$\boxed{2} \qquad \begin{bmatrix}
	1 & 1 & 0 \\
	1 & 0 & 1 \\
	1 & 0 & 0 \end{bmatrix} \begin{bmatrix}
	a_2 \\ b_2 \\ c_2 \end{bmatrix} = \begin{bmatrix}
	0 \\ 1 \\ 0 \end{bmatrix} \Rightarrow a_2 = 0, \quad b_2 = 0, \quad c_2 = 1 \Rightarrow N_2 = \eta$
	
$\boxed{3} \qquad \begin{bmatrix}
	1 & 1 & 0 \\
	1 & 0 & 1 \\
	1 & 0 & 0 \end{bmatrix} \begin{bmatrix}
	a_3 \\ b_3 \\ c_3 \end{bmatrix} = \begin{bmatrix}
	0 \\ 0 \\ 1 \end{bmatrix} \Rightarrow a_3 = 1, \quad b_3 = -1, \quad c_3 = -1 \Rightarrow N_3 = 1-\xi-\eta$
	
""", blur=false)

# ╔═╡ f5de8324-51ae-11eb-19c3-334aa018ec50
md"""
## Aufgabe 8.3b) 
 Wie lauten die Ansatzfunktionen $N_1$ und $N_6$ für das quadratische Dreieckselement aus der Abbildung?
"""

# ╔═╡ 9e27e4f8-51af-11eb-044b-51a6179e9f02
hint(md"""
Benutzen Sie den modifizierter Lagrange-Ansatz
	
$N = l_{k^1}^{k^1}(\lambda^1)l_{k^2}^{k^2}(\lambda^2)l_{k^3}^{k^3}(\lambda^3)$
	
$\lambda^1 = \xi$
	
$\lambda^2 = \eta$
	
$\lambda^3 = 1-\xi-\eta$
""", blur=false)

# ╔═╡ 8df86c16-566b-11eb-3cb8-5542f33ef968
md"""
### N₁ quadratisches Dreieck:
"""

# ╔═╡ db4d5714-51c8-11eb-11f6-8df484766689
function N₁tq(ξ,η)
    return ξ*(2*ξ-1)
end

# ╔═╡ 6e6da7ae-51d0-11eb-291c-619325c4dd33
begin
	ξ8 = collect(0:0.01:1)
	η8 = collect(0.:0.01:1)
	Plots.surface(ξ8,η8, (x1,x2) -> (1 - x1 - x2) > 0.0 ? N₁tq(x1,x2) : 0.0,camera = (10, 40)) 
end

# ╔═╡ a5335596-566b-11eb-0e12-8f81e4a31037
md"""
### N₆ quadratisches Dreieck:
"""

# ╔═╡ dc6ad310-51c8-11eb-1b7f-3509217ae646
function N₆tq(ξ,η)
    return 4*ξ*(1-ξ-η)
end

# ╔═╡ 134636ee-566b-11eb-27f1-a190f410d987
begin
	ξ9 = collect(0:0.05:1)
	η9 = collect(0.:0.05:1)
	Plots.surface(ξ9,η9, (x1,x2) -> (1 - x1 - x2) > 0.0 ? N₆tq(x1,x2) : 0.0,camera = (40, 50)) 
end

# ╔═╡ 28a643a0-51b4-11eb-1385-db4581abfd7f
solution(md"""
$\boxed{1} \quad \lambda^1: \quad k_1 = 2, \quad i \in \{0,1\}, \quad \lambda_k^1 = 1, \quad \lambda_i^1 \in \{0,\frac{1}{2}\}$

$\qquad \lambda^2: \quad k_2 = 0 \Rightarrow l_0^0(\lambda^2) = 1$
	
$\qquad \lambda^3: \text{wie } \lambda^2$
	
$N_1 = l_2^2(\lambda^1)l_0^0(\lambda^2)l_0^0(\lambda^3) = \frac{(\lambda^1 - 0)(\lambda^1 - \frac{1}{2})}{(1-0)(1-\frac{1}{2})}.1.1 = \lambda^1(2\lambda^1-1)$
	
$\boxed{6} \quad \lambda^1: \quad k_1 = 1, \quad i \in \{0\}, \quad \lambda_k^1 = \frac{1}{2}, \quad \lambda_i^1 \in \{0\}$

$\qquad \lambda^2: \quad k_2 = 0 \Rightarrow l_0^0(\lambda^2) = 1$
	
$\qquad \lambda^3: \text{wie } \lambda^1$
	
$N_1 = l_1^1(\lambda^1)l_0^0(\lambda^2)l_1^1(\lambda^3) = \frac{(\lambda^1 - 0)(\lambda^3 - 0)}{(\frac{1}{2}-0)(\frac{1}{2}-0)} = 4\lambda^1\lambda^3$

""",blur=false)

# ╔═╡ f7b23924-51b9-11eb-069a-052683356115
md"""
## Aufgabe 8.3c) 
 Wie lauten die Ansatzfunktionen für ein dreidimensionales lineares Tetraederelement?
"""

# ╔═╡ cf604c22-59ac-11eb-159b-2f9e83ff0b8b
HTML(open(f->read(f, String), "assets/Tet.svg"))

# ╔═╡ 51e54382-51ba-11eb-2641-63d832fc16df
hint(md"""
allgemeine Ansatz 
$$N_I = a_I + b_I\xi + c_I\eta + d_I\delta$$
	
$\begin{bmatrix}
	1 & 1 & 0 & 0 \\
	1 & 0 & 1 & 0 \\
	1 & 0 & 0 & 1 \\
	1 & 0 & 0& 0 \end{bmatrix} \begin{bmatrix}
	a_I \\ b_I \\ c_I \\ d_I \end{bmatrix} = \delta_{IJ}$
""", blur=false)

# ╔═╡ 4435ca94-566b-11eb-25f6-e7bca5ba7e11
md"""
### N₁ Tetrahedron:
"""

# ╔═╡ b8bf9e4e-566b-11eb-3fb2-83e31f53284f
[1 1 0 0
 1 0 1 0
 1 0 0 1
 1 0 0 0]\[1,0,0,0]

# ╔═╡ 3842cf58-51c9-11eb-0714-17cd79007491
function N₁t³(ξ,η,δ)
    return ξ
end

# ╔═╡ e256f95a-566b-11eb-2547-cb02887363b6
md"""
### N₂ Tetrahedron:
"""

# ╔═╡ e7ed2556-566b-11eb-1f41-cda5daf72b49
[1 1 0 0
 1 0 1 0
 1 0 0 1
 1 0 0 0]\[0,1,0,0]

# ╔═╡ 580cd77a-51c9-11eb-23fe-03dd3735b4d6
function N₂t³(ξ,η,δ)
    return η
end

# ╔═╡ f2785b14-566b-11eb-131e-a71652a0ac35
md"""
### N₃ Tetrahedron:
"""

# ╔═╡ fc549526-566b-11eb-217c-1b0e941c88b1
[1 1 0 0
 1 0 1 0
 1 0 0 1
 1 0 0 0]\[0,0,1,0]

# ╔═╡ 6058709a-51c9-11eb-03fb-b13e0d47ab46
function N₃t³(ξ,η,δ)
    return δ
end

# ╔═╡ 0937d9a4-566c-11eb-24c2-7742e5a2f3df
md"""
### N₄ Tetrahedron:
"""

# ╔═╡ 12e78a58-566c-11eb-33c5-1dba9e71b67e
[1 1 0 0
 1 0 1 0
 1 0 0 1
 1 0 0 0]\[0,0,0,1]

# ╔═╡ 6b82b4f0-51c9-11eb-1d0e-c5266583c2e7
function N₄t³(ξ,η,δ)
    return 1-ξ-η-δ
end

# ╔═╡ 8a4a5af0-51ba-11eb-1660-7b44765dfce7
solution(md"""
$\boxed{1} \qquad \begin{bmatrix}
	1 & 1 & 0 & 0 \\
	1 & 0 & 1 & 0 \\
	1 & 0 & 0 & 1 \\
	1 & 0 & 0& 0 \end{bmatrix} \begin{bmatrix}
	a_1 \\ b_1 \\ c_1 \\ d_1 \end{bmatrix} = \begin{bmatrix}
	1 \\ 0 \\ 0 \\ 0 \end{bmatrix} \Rightarrow a_1 = c_1 = d_1 = 0, \quad b_1 = 1 \Rightarrow N_1 = \xi := \lambda^1$
	
$\boxed{2} \qquad \begin{bmatrix}
	1 & 1 & 0 & 0 \\
	1 & 0 & 1 & 0 \\
	1 & 0 & 0 & 1 \\
	1 & 0 & 0& 0 \end{bmatrix} \begin{bmatrix}
	a_2 \\ b_2 \\ c_2 \\ d_2 \end{bmatrix} = \begin{bmatrix}
	0 \\ 1 \\ 0 \\ 0 \end{bmatrix} \Rightarrow a_2 = b_2 = d_2 = 0, \quad c_2 = 1 \Rightarrow N_2 = \eta := \lambda^2$
	
$\boxed{3} \qquad \begin{bmatrix}
	1 & 1 & 0 & 0 \\
	1 & 0 & 1 & 0 \\
	1 & 0 & 0 & 1 \\
	1 & 0 & 0& 0 \end{bmatrix} \begin{bmatrix}
	a_3 \\ b_3 \\ c_3 \\ d_3 \end{bmatrix} = \begin{bmatrix}
	0 \\ 0 \\ 1 \\ 0 \end{bmatrix} \Rightarrow a_3 = b_3 = c_3 = 0, \quad d_3 = 1 \Rightarrow N_3 = \delta := \lambda^3$
	
$\boxed{4} \quad \begin{bmatrix}
	1 & 1 & 0 & 0 \\
	1 & 0 & 1 & 0 \\
	1 & 0 & 0 & 1 \\
	1 & 0 & 0& 0 \end{bmatrix} \begin{bmatrix}
	a_4 \\ b_4 \\ c_4 \\ d_4 \end{bmatrix} = \begin{bmatrix}
	0 \\ 0 \\ 0 \\ 1 \end{bmatrix}$

$\Rightarrow a_4 = 1, \: b_4 = c_4 = d_4 = -1, \Rightarrow N_4 =1-\xi-\eta- \delta := \lambda^4$	
""", blur=false)

# ╔═╡ Cell order:
# ╟─47277816-518d-11eb-244f-fd99d21c2a6e
# ╟─96d1e7b6-596c-11eb-1305-716af2ea8390
# ╟─1d876196-518c-11eb-0144-d3ec17e099b4
# ╟─f78f5792-518c-11eb-1f4d-5bdf07c2ae20
# ╟─8fe076c0-518d-11eb-10ae-2f725b2f0a05
# ╟─50078ece-59a0-11eb-2945-4376d914232c
# ╟─b6d4fe8a-518e-11eb-1783-179c0af183d0
# ╟─86dc8bde-59a2-11eb-1e90-bd4a56e53b43
# ╟─dccd7022-518e-11eb-2d7d-990c9ad72a32
# ╟─f952cc6a-518e-11eb-3cfa-2be6123c5271
# ╟─6c07fbe0-518f-11eb-351f-1b138a4d859c
# ╟─833ffd44-518f-11eb-261f-359001e33f7d
# ╟─44703a2e-59b5-11eb-0176-61aee94d5348
# ╟─95d3074e-518f-11eb-09a7-1b8cf299caf0
# ╟─af4b8480-518f-11eb-2b47-334c93cf76ad
# ╟─d1784dc6-519f-11eb-0d92-43feaf202c19
# ╟─a5eae37e-51a1-11eb-248e-0ff3ac1a4218
# ╟─c62c5492-51a1-11eb-259c-f5cb7cbf35cf
# ╟─65c1c186-51a2-11eb-224a-b3c2e521b61b
# ╟─facff0ca-51a4-11eb-0918-99ec861bfde4
# ╠═2976907e-51a2-11eb-03c9-5597c9dff33d
# ╟─57fa8e8a-59a8-11eb-31b2-db3166b9ebe2
# ╟─418267fe-51a6-11eb-2afd-3748ac675115
# ╠═1e393458-51a6-11eb-3e4f-1dc7b554177b
# ╟─9d9c4c0c-51a2-11eb-0350-83efe127a32f
# ╟─99a2d9c8-51a6-11eb-2e69-2b5ca5565d69
# ╟─57a7f00c-51c5-11eb-216e-1f2c676d81ea
# ╠═f6ae9760-51c4-11eb-3861-87e7030628a8
# ╟─ca71f52e-51c5-11eb-0597-097eaed7b7b2
# ╠═a50129e0-51c5-11eb-3126-6345cf39d287
# ╟─53c5d182-51a9-11eb-3ce0-59bc172b982d
# ╟─5b81e7e8-51aa-11eb-1af1-8b08c6d9b441
# ╟─0c007e40-51ab-11eb-03c2-afb37f8f292e
# ╟─cd0c8f48-51ab-11eb-12f1-09e5c6402e21
# ╟─138b87f4-51ba-11eb-1933-7901c56f9b24
# ╟─61a0467e-5669-11eb-1387-19be37b3cc37
# ╠═a4b48c88-5669-11eb-0dc4-1506e071610a
# ╟─9b77a9c2-51c9-11eb-1d15-1f104be60183
# ╠═2d3c53fa-51c8-11eb-064d-cde446cc7d72
# ╟─3db4ba64-566a-11eb-2230-bbe82a13b8ca
# ╠═4654300a-566a-11eb-313f-03cdb24885fa
# ╟─58147a8e-566a-11eb-3311-dd04ee91c1f7
# ╠═b5b78574-51c8-11eb-0958-b53010124aa1
# ╟─a2d0df5e-566a-11eb-1c58-d5631e9c7635
# ╠═b07170ba-566a-11eb-3244-0d0fd670f8df
# ╟─bb652a84-566a-11eb-25f8-83c449659164
# ╠═c3c89a42-51c8-11eb-2b9d-b1bcdc3f6690
# ╟─18d09e7e-51ac-11eb-0dfc-5ffefbcafdc0
# ╟─f5de8324-51ae-11eb-19c3-334aa018ec50
# ╟─9e27e4f8-51af-11eb-044b-51a6179e9f02
# ╟─8df86c16-566b-11eb-3cb8-5542f33ef968
# ╟─6e6da7ae-51d0-11eb-291c-619325c4dd33
# ╠═db4d5714-51c8-11eb-11f6-8df484766689
# ╟─a5335596-566b-11eb-0e12-8f81e4a31037
# ╟─134636ee-566b-11eb-27f1-a190f410d987
# ╠═dc6ad310-51c8-11eb-1b7f-3509217ae646
# ╟─28a643a0-51b4-11eb-1385-db4581abfd7f
# ╟─f7b23924-51b9-11eb-069a-052683356115
# ╟─cf604c22-59ac-11eb-159b-2f9e83ff0b8b
# ╟─51e54382-51ba-11eb-2641-63d832fc16df
# ╟─4435ca94-566b-11eb-25f6-e7bca5ba7e11
# ╠═b8bf9e4e-566b-11eb-3fb2-83e31f53284f
# ╠═3842cf58-51c9-11eb-0714-17cd79007491
# ╟─e256f95a-566b-11eb-2547-cb02887363b6
# ╠═e7ed2556-566b-11eb-1f41-cda5daf72b49
# ╠═580cd77a-51c9-11eb-23fe-03dd3735b4d6
# ╟─f2785b14-566b-11eb-131e-a71652a0ac35
# ╠═fc549526-566b-11eb-217c-1b0e941c88b1
# ╠═6058709a-51c9-11eb-03fb-b13e0d47ab46
# ╟─0937d9a4-566c-11eb-24c2-7742e5a2f3df
# ╠═12e78a58-566c-11eb-33c5-1dba9e71b67e
# ╠═6b82b4f0-51c9-11eb-1d0e-c5266583c2e7
# ╟─8a4a5af0-51ba-11eb-1660-7b44765dfce7
