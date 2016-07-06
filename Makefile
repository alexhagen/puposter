
all: figures sample

figures: figure1 figure3 figure4 figure5 figure7

figure1: py/mctmfd_eff.py
	cd py; \
	python mctmfd_eff.py

figure3: py/panel_analysis.py
	cd py; \
	python panel_analysis.py

figure4: py/ef15_eff.py
	cd py; \
	python ef15_eff.py

figure5: py/datmfd_analysis.py
	cd py; \
	python datmfd_analysis.py

figure7: py/photoneutron_eff.py
	cd py; \
	python photoneutron_eff.py
	
figure8: py/clinac_plot.py
	cd py; \
	python clinac_plot.py

figure9: py/ddaa.py
	cd py; \
	python ddaa.py

sample: sample.tex
	xelatex sample.tex
	xelatex sample.tex
	rm -f *.aux *.dvi *.log *.mw *.cpsheadings
	/Applications/Adobe\ Acrobat\ Reader\ DC.app/Contents/MacOS/AdobeReader sample.pdf
