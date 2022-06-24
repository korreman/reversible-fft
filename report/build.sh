pdflatex -halt-on-error -shell-escape master.tex;
bibtex master.aux;
pdflatex -halt-on-error -shell-escape master.tex;
pdflatex -halt-on-error -shell-escape master.tex;
