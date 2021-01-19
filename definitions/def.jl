using PlutoUI
using StringDistances
using Luxor

hint(text;blur=true) = blur ? Markdown.MD(Markdown.Admonition("hint", "Hilfe", [text])) : Markdown.MD(Markdown.Admonition("note", "Hilfe", [text])) 
warning(text) = Markdown.MD(Markdown.Admonition("warning", "Warnung", [text]))
yays = [md"Sehr gut! 🐣", md"Yay ❤", md"Genau so! 🎉", md"Gut gemacht! 🐦", md"Weiter so! 🐤", md"Klasse! 🐧", md"Korrekt! 🐖", md"Sehr schön! 🐿"]
correct(text=rand(yays)) = Markdown.MD(Markdown.Admonition("correct", "Richtig!", [text]))
still_missing(text=md"Ersetzen sie `missing` mit ihrer Antwort") = Markdown.MD(Markdown.Admonition("warning", "Here we go! 🦦", [text]))
keep_working(text=md"Noch nicht die richtige Antwort, noch ein Versuch! 🦥") = Markdown.MD(Markdown.Admonition("danger", "Falsch", [text]))
solution(text; blur=true) = blur ? Markdown.MD(Markdown.Admonition("hint", "Lösung", [text])) : Markdown.MD(Markdown.Admonition("note", "Lösung", [text]))

function stringmatrix_to_mdmatrix(stringmatrix; latextype="bmatrix")
	prefix = "\\begin{$latextype}"
	suffix = "\\end{$latextype}"
	result = ""
	result *= prefix
	rows = lastindex(stringmatrix,1)
	cols = lastindex(stringmatrix,2)
	for i in 1:rows
		row = "" 
		for j in 1:cols
			row *= isa(stringmatrix[i,j], String) ? strip(stringmatrix[i,j]) : "$(stringmatrix[i,j])"
			if j != cols 
				row *= "&" 
			end
		end
		row *= "\\\\" 
		result *= row
	end
	return result *= suffix
end

function stringmatrix_vector(stringmatrix, stringvector;latextype="bmatrix")
   matrix_string = stringmatrix_to_mdmatrix(stringmatrix; latextype)
   vector_string = stringmatrix_to_mdmatrix(stringvector; latextype)
   return matrix_string * "\\qquad" * vector_string
end

function stringmatrix_vector(stringmatrix, stringvector1, stringvector2, operator;latextype="bmatrix")
   matrix_string = stringmatrix_to_mdmatrix(stringmatrix; latextype)
   vector_string1 = stringmatrix_to_mdmatrix(stringvector1; latextype)
   vector_string2 = stringmatrix_to_mdmatrix(stringvector2; latextype)
   return matrix_string * "\\cdot " * vector_string1 * operator * vector_string2
end

show_matrix(stringmatrix;prefix="", suffix="")=Markdown.MD(Markdown.LaTeX(prefix*stringmatrix_to_mdmatrix(stringmatrix)*suffix))
show_matrix(stringmatrix,stringvector) = Markdown.MD(Markdown.LaTeX(stringmatrix_vector(stringmatrix, stringvector)))
show_matrix(stringmatrix,stringvector1, stringvector2, operator) = Markdown.MD(Markdown.LaTeX(stringmatrix_vector(stringmatrix, stringvector1, stringvector2, operator)))

"""
    draw_structure(elements::Array) -> svg
Zeichnet die Struktur anhand eines `Array{Element,1}`

# Argument
- `elements::Array`: die Elemente die gezeichnet werden sollen

# Keyword Argumente
- `elebox_size=15`
- `drawing_size=(650,200)`
- `origin_point=(50,150)`
- `scale=100`
- `node_label_offset=15`
- `ele_label_offset=-3.5`
- `support_offset_x_dir = Point(-10,0)`
- `support_offset_y_dir = Point(0,10)`
- `support_size=10`
"""
function draw_structure(elements::Array; elebox_size=15, drawing_size=(650,200),origin_point=(50,150), 
    scale=100, node_label_offset=15, ele_label_offset=-3.5, support_offset_x_dir = Point(-10,0), support_offset_y_dir = Point(0,10), support_size=10)
    labeled = zeros(Bool, length(elements),2)
    Drawing(drawing_size[1],drawing_size[2], :svg)
    origin(origin_point[1],origin_point[2])
    for (id,element) in enumerate(elements)
        @layer begin
            p1 = Point(element.n1.coord[1]*scale, -element.n1.coord[2]*scale)
            p2 = Point(element.n2.coord[1]*scale, -element.n2.coord[2]*scale)
            line(p1, p2, :stroke)
            circle(p1,5,:fill)
            circle(p2,5,:fill)
        end
        
        if !labeled[id,1]
            label("$(element.n1.id)", :N, p1, offset=node_label_offset)
            labeled[id,1] = true
        end
        if !labeled[id,2]
            label("$(element.n2.id)", :N, p2, offset=node_label_offset)
            labeled[id,2] = true
        end
        sethue("white")
        middle = (p1 + p2) / 2
        box(middle, elebox_size, elebox_size, :fill)
        sethue("black")
        box(middle, elebox_size, elebox_size, :stroke)
        label("$(element.enum)", :N, middle, offset=ele_label_offset)
        
        support_offset_x = support_offset_x_dir
        support_offset_y = support_offset_y_dir
        if false in element.n1.bcs
            if !element.n1.bcs[1] && element.n1.bcs[2]
                ngon(p1+support_offset_x, support_size, 3, 0, :stroke)
            elseif !element.n1.bcs[2] && element.n1.bcs[1]
                ngon(p1+support_offset_y, support_size, 3, -pi/2, :stroke)
            else
                ngon(p1+support_offset_y, support_size, 3, -pi/2, :stroke)
                ngon(p1+support_offset_x, support_size, 3, 0, :stroke)
            end
        end
        if false in element.n2.bcs
            if !element.n2.bcs[1] && element.n2.bcs[2]
                ngon(p2+support_offset_x, support_size, 3, 0, :stroke)
            elseif !element.n2.bcs[2] && element.n2.bcs[1]
                ngon(p2+support_offset_y, support_size, 3, -pi/2, :stroke)
            else
                ngon(p2+support_offset_y, support_size, 3, -pi/2, :stroke)
                ngon(p2+support_offset_x, support_size, 3, 0, :stroke)
            end
        end
    end
    finish()
    preview()
end

nothing
