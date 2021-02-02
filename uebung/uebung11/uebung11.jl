### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# ╔═╡ e6b03cba-6565-11eb-1696-753338d1dc37
using JuAFEM

# ╔═╡ 91f2c318-5b1b-11eb-3ec0-536aff8a0ea7
include("../../definitions/def.jl")

# ╔═╡ d58d0ee6-659e-11eb-1dd1-55cf264c91da
html"""<iframe width="100%" height="400" src="https://www.youtube.com/embed/_L-xAbwYuwo" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>"""

# ╔═╡ fd4a9bb6-5b09-11eb-00e1-65a2a2bfa41b
md"""
# Gauss-Quadratur
Willkommen in der letzten regulären Übung im Modul Grundlagen der FEM. In dieser Übung werden wir uns die Gauss-Quadratur genauer anschauen.

Die Formel für die numerische Integration mittels Gauss-Quadratur lautet

$\int_{-1}^{1}f(\xi)d\xi \approx \sum_{i=1}^n w_if(\xi_i)$

wobei wir für FE Elemente häufig folgende Relation nutzen

$\int_{\Omega}f(\mathbf{x})\ \text{d}V \approx \sum_{i=1}^n w_if(\xi_i) \det \mathbf{J}$

Dabei sind $w_i$ die Wichungsfaktoren und $\xi_i$ die Gausspunkte.
"""

# ╔═╡ 776d1262-5b0b-11eb-2512-f348a02d41ac
md"""
## Aufgabe 11.1 - Quadratisches Viereckselement
Für das quadratische isoparametrische Viereckselement wird die Gaussintegrationsordnung $n=3$ gewählt.
"""

# ╔═╡ fc8ef012-5b1a-11eb-21a2-8b9036e30fc5
HTML(open(f->read(f, String), "assets/gauss-points.svg"))

# ╔═╡ 12c052a8-5b0d-11eb-061c-cb1a70bbf0c9
md"""
### Aufgabe 1.a) 
Wie viele Gausspunkte hat das Element dann?
"""

# ╔═╡ 831eabf8-5b1c-11eb-2c9b-6bf730631a18
GP = missing

# ╔═╡ a6c179dc-5b1c-11eb-1ea7-2b795958277f
try
if GP == 9
	 correct1 = true
	 correct()
else
	 correct1 = false
	 warning(md"""Falsche oder nicht gegebene Antwort 🦥""")
end
catch
	still_missing()
end

# ╔═╡ 6b76f6a0-5b1b-11eb-3191-09b38930e4b7
solution(md"""
n Gauss-points in $\xi$-richtung 
	
n Gauss-points in $\eta$-richtung 
	
Insgesamt	
$\Rightarrow$ nxn = 9 Gauss-points
""")

# ╔═╡ 36a3cc9a-5b0d-11eb-0869-978659cccfb3
md"""
### Aufgabe 1.b)

Ist die Gaussintegration der Ordnung $n=3$ für eine quadratische Funktion $f(\xi)$ (Polynomgrad $p=2$) exakt? 
Welche Bedingung muss erfüllt sein?
          
Hinweis: Da bei Viereckselementen die Übertragung von einer auf zwei Koordinatenrichtungen direkt möglich ist reicht es in dieser Aufgabe eine Funktion $f(\xi)$ mit nur einer Koordinate zu betrachten.
"""

# ╔═╡ d79955c2-5b1f-11eb-2c7a-15c846bc640c
function min_gauss_order(polynomgrad::Int64)
	return ceil((polynomgrad+1)/2)
end

# ╔═╡ 2f1abbe4-5b20-11eb-2555-e98ebfe040b5
min_gauss_order()

# ╔═╡ e52fdfc4-5b1c-11eb-2a0b-7b5b7dfdb234
solution(md"""
Betrachtung einzelner Koordinate $f(\xi)$
	
f: Polynom mit Polynomgrad $p=2$
	
	
$p \le 2n-1$ 
	
$2 \le 2\cdot 3 - 1 = 5$
	
welche Integrationsordung für welchen Polynomgrad? 
	
$n=1 \Rightarrow 2\cdot 1-1 = 1 \Rightarrow p \le 1$
	
$n=2 \Rightarrow 2\cdot 2-1 = 3 \Rightarrow p \le 3$
	
$n=3 \Rightarrow 2\cdot 3-1 = 5 \Rightarrow p \le 5$
	
| p (Polynomgrad)       |      | n (Gauss Integrationsordnung)    |
| :-------------: | :----------: | :-----------: |
|  1 |    | 1    |
|  2 |    | 2    |
|  3 |    | 2    |
|  4 |    | 3    |
|  5 |    | 3   |

""")

# ╔═╡ 8ac8756e-5b0d-11eb-347b-f7344b7592d7
md"""
### Aufgabe 1.c)
Gibt es Gründe eine höhere Integrationsordnung als nötig zu verwenden? Erklären Sie ihre Antwort anhand dieses Beispiels
"""

# ╔═╡ d9fe9d9c-5b27-11eb-1c7f-9f109507856f
HTML(open(f->read(f, String), "assets/fall.svg"))

# ╔═╡ 7bfd01fe-5b20-11eb-262a-11e0d0928d51
solution(md""" 
* Kein Gewinn an Genauigkeit, Rechenaufwand steigt 
* Kann dennoch sinnvoll sein um einen Rangabfall der globalen Steifigkeitsmatrix zu verhindern
	
Rangabfall:

$s < n_{dof}$
-$s$: Linear unabhängige Dehnungsgleichungen
	
-$n_{dof}$: Anzahl globalen Freiheitsgrade

In dem gegebenen Beispiel 
	
$$n_{dof} = 2\cdot 9-3 = 15$$ 
	
$$s = 4 \cdot 3 = 12 < 15$$
	
wenn man höherer Gauss Integrationordnung benutzt ergibt das für $s$

$s = 9\cdot 3 = 27 > 15$
""")

# ╔═╡ bd169d52-5b0d-11eb-1c14-35b8a694e69c
md"""
Findet die Gaussintegration über das Intervall $[-1,1]$ statt so lassen sich die Gausspunke $\xi_i\ (i=1,\ldots,n)$ als Nullstellen des n-ten Legendrepolynoms bestimmen:

$P_n(\xi)=\frac{1}{2^n\,n!}\ \frac{{d}^n}{{d}\xi^n}(\xi^2-1)^n$
"""

# ╔═╡ f855016a-5b0d-11eb-1ff8-e5208dd3f631
md"""
### Aufgabe 1.d)
Wie lautet das Legendrepolynom 3. Ordnung?
"""

# ╔═╡ db8b12ba-5b29-11eb-1f24-63d8f5aafbd9
p3 = md"""
$$1$$
"""

# ╔═╡ 5803d2d2-5b2a-11eb-3ec4-332aad488d5e
begin
true_p3 = md"""
$\frac{1}{2}(5\xi^2-3)\xi$
"""
	nothing
end

# ╔═╡ fb9952ec-5b29-11eb-2226-9d8bf2487c5e
try
if p3 == true_p3
	 correct2 = true
	 correct()
else
	 correct2 = false
	 warning(md"""Falsche oder nicht gegebene Antwort 🦥""")
end
catch
	still_missing()
end

# ╔═╡ b1822af4-5b28-11eb-2e36-e95d365494e5
solution(md"""
$p_3 = \frac{1}{2^3 3!}\frac{d^3}{d\xi^3}(\xi^2-1)^3$

$\frac{d^3}{d\xi^3}(\xi^2-1)^3 = 120\xi^3 - 72\xi$
	
$\Rightarrow p_3 = \frac{1}{2}(5\xi^2-3)\xi$
""")

# ╔═╡ 0ed3b2ce-5b0e-11eb-0e5b-1bd44c706fd8
md"""
### Aufgabe 1.e)
Bestimmen Sie die entsprechenden Gausspunktkoordinaten $\xi_i$
"""

# ╔═╡ 90fedb7a-5b2c-11eb-05e7-45752e237bd6
ξ = Set([])

# ╔═╡ 9d70be1e-5b2c-11eb-08b6-ab31d6976cd8
try
if ξ == Set([-sqrt(3/5),0.0,sqrt(3/5)])
	 correct3 = true
	 correct()
else
	 correct3 = false
	 warning(md"""Falsche oder nicht gegebene Antwort 🦥""")
end
catch
	still_missing()
end

# ╔═╡ 0411230a-5b2c-11eb-2e37-9f2bad9fe71f
solution(md"""
$\frac{1}{2}(5\xi^2-3)\xi = 0$
	
$\xi = 0, \quad 5\xi^2 = 3 \Rightarrow \xi = \pm \sqrt{\frac{3}{5}}$
	
$\Rightarrow \xi \in \{-\sqrt{\frac{3}{5}},0,\sqrt{\frac{3}{5}}\}$
""")

# ╔═╡ 1dfa74e0-5b0e-11eb-104b-6f14cd7c3a7d
md"""
Die Zugehörigen Wichtungsfaktoren $w_i$ lassen sich über folgende Beziehung bestimmen:

$w_i=\int_{-1}^1 l_i(\xi)\ d\xi$

Dabei ist $l_i$ das i-te Lagrangepolynom (Polynomgrad: $n-1$) mit den Gausspunkten $\xi_i$ als Stützstellen.
"""

# ╔═╡ 4259dac4-5b0e-11eb-20c2-6b63027093c2
md"""
### Aufgabe 1.f)
Bestimmen Sie die zu den Gausspunktkoordinaten gehörigen Wichtungsfaktoren $w_i$
"""

# ╔═╡ 7b8693b4-5b2f-11eb-2639-bffc95f7117a
w₁ = missing

# ╔═╡ 858d0bce-5b2f-11eb-17d0-e37d7613e09f
w₂ = missing

# ╔═╡ 89f7f8cc-5b2f-11eb-3856-33ae5716a418
w₃ = missing

# ╔═╡ 8fefa20c-5b2f-11eb-1819-61de241a247d
try
if w₁ == 5/9 && w₂ == 8/9 && w₃ == 5/9
	 correct4 = true
	 correct()
else
	 correct4 = false
	 warning(md"""Falsche oder nicht gegebene Antwort 🦥""")
end
catch
	still_missing()
end

# ╔═╡ 5f9997c0-5b2d-11eb-2d0a-934dff5a1f50
solution(md"""
$l_1 = \frac{(\xi-0)(\xi-\sqrt{\frac{3}{5}})}{(-\sqrt{\frac{3}{5}} - 0)(-\sqrt{\frac{3}{5}} - \sqrt{\frac{3}{5}})} = \frac{5}{6}\xi^2 - \frac{5}{6}\sqrt{\frac{3}{5}}\xi$
	
$w_1 = \int_{-1}^{1} l_1dx = \left[ \frac{5}{18}\xi^3 - \frac{5}{12}\sqrt{\frac{3}{5}}\xi^2 \right]_{-1}^1 = \frac{5}{9}$
	
$l_2 = \frac{(\xi+\sqrt{\frac{3}{5}})(\xi-\sqrt{\frac{3}{5}})}{(0-\sqrt{\frac{3}{5}})(0+\sqrt{\frac{3}{5}})} = \frac{\xi^2 - \frac{3}{5}}{-\frac{3}{5}} = -\frac{5}{3}\xi^2 + 1$
	
$w_2 = \int_{-1}^{1} l_2dx = \left[ -\frac{5}{9}\xi^3 + \xi \right]_{-1}^1 = \frac{8}{9}$
	
$l_3 = (\xi+\sqrt{\frac{3}{5}})\frac{(\xi-0)}{(\sqrt{\frac{3}{5}} + \sqrt{\frac{3}{5}})(\sqrt{\frac{3}{5}} - 0)} = \frac{5}{6}\xi^2 + \frac{5}{6}\sqrt{\frac{3}{5}}\xi$
	
$w_3 = \int_{-1}^{1} l_3dx = \left[ \frac{5}{18}\xi^3 + \frac{5}{12}\sqrt{\frac{3}{5}}\xi^2 \right]_{-1}^1 = \frac{5}{9}$
""")

# ╔═╡ db053212-6565-11eb-2841-57c36434f811
QuadratureRule{1,RefCube}(3)

# ╔═╡ 20e8a5f2-6566-11eb-1245-f7bc2e3b4b7a
(5/9)*(8/9)

# ╔═╡ 14f9fd36-6566-11eb-3e91-a1aa7beddd66
QuadratureRule{2,RefCube}(3)

# ╔═╡ Cell order:
# ╟─91f2c318-5b1b-11eb-3ec0-536aff8a0ea7
# ╟─d58d0ee6-659e-11eb-1dd1-55cf264c91da
# ╟─fd4a9bb6-5b09-11eb-00e1-65a2a2bfa41b
# ╟─776d1262-5b0b-11eb-2512-f348a02d41ac
# ╟─fc8ef012-5b1a-11eb-21a2-8b9036e30fc5
# ╟─12c052a8-5b0d-11eb-061c-cb1a70bbf0c9
# ╠═831eabf8-5b1c-11eb-2c9b-6bf730631a18
# ╟─a6c179dc-5b1c-11eb-1ea7-2b795958277f
# ╟─6b76f6a0-5b1b-11eb-3191-09b38930e4b7
# ╟─36a3cc9a-5b0d-11eb-0869-978659cccfb3
# ╠═d79955c2-5b1f-11eb-2c7a-15c846bc640c
# ╠═2f1abbe4-5b20-11eb-2555-e98ebfe040b5
# ╟─e52fdfc4-5b1c-11eb-2a0b-7b5b7dfdb234
# ╟─8ac8756e-5b0d-11eb-347b-f7344b7592d7
# ╟─d9fe9d9c-5b27-11eb-1c7f-9f109507856f
# ╟─7bfd01fe-5b20-11eb-262a-11e0d0928d51
# ╟─bd169d52-5b0d-11eb-1c14-35b8a694e69c
# ╟─f855016a-5b0d-11eb-1ff8-e5208dd3f631
# ╠═db8b12ba-5b29-11eb-1f24-63d8f5aafbd9
# ╟─5803d2d2-5b2a-11eb-3ec4-332aad488d5e
# ╟─fb9952ec-5b29-11eb-2226-9d8bf2487c5e
# ╟─b1822af4-5b28-11eb-2e36-e95d365494e5
# ╟─0ed3b2ce-5b0e-11eb-0e5b-1bd44c706fd8
# ╠═90fedb7a-5b2c-11eb-05e7-45752e237bd6
# ╟─9d70be1e-5b2c-11eb-08b6-ab31d6976cd8
# ╟─0411230a-5b2c-11eb-2e37-9f2bad9fe71f
# ╟─1dfa74e0-5b0e-11eb-104b-6f14cd7c3a7d
# ╟─4259dac4-5b0e-11eb-20c2-6b63027093c2
# ╠═7b8693b4-5b2f-11eb-2639-bffc95f7117a
# ╠═858d0bce-5b2f-11eb-17d0-e37d7613e09f
# ╠═89f7f8cc-5b2f-11eb-3856-33ae5716a418
# ╟─8fefa20c-5b2f-11eb-1819-61de241a247d
# ╟─5f9997c0-5b2d-11eb-2d0a-934dff5a1f50
# ╠═e6b03cba-6565-11eb-1696-753338d1dc37
# ╠═db053212-6565-11eb-2841-57c36434f811
# ╠═20e8a5f2-6566-11eb-1245-f7bc2e3b4b7a
# ╠═14f9fd36-6566-11eb-3e91-a1aa7beddd66
