
all: sample

sample: sample.tex
	xelatex sample.tex
	xelatex sample.tex
	rm -f *.aux *.dvi *.log
	open sample.pdf
