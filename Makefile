# Common tasks. Override the interpreter with: make PY=python3 <target>
PY ?= .venv/bin/python

# LaTeX: latexmk if installed, otherwise pdflatex + bibtex + pdflatex x2
ifneq ($(shell command -v latexmk 2>/dev/null),)
TEX = latexmk -pdf -interaction=nonstopmode
else
TEX = sh -c 'pdflatex -interaction=nonstopmode "$$1" && bibtex "$${1%.tex}" && pdflatex -interaction=nonstopmode "$$1" && pdflatex -interaction=nonstopmode "$$1"' _
endif

.PHONY: help venv test results figures paper report docs slides all clean check-tex

help:
	@echo "venv     create .venv and install simulation/requirements.txt"
	@echo "test     fast regression tests (~10 s)"
	@echo "results  rerun all studies -> simulation/results (~45 min)"
	@echo "figures  regenerate paper figures from results (seconds)"
	@echo "paper    build paper/manuscript/main.pdf"
	@echo "report   build paper/background_report/report.pdf"
	@echo "docs     build simulation/docs/code_guide.pdf (code documentation)"
	@echo "slides   build presentation/radial_eit_cooling_talk.pptx"
	@echo "clean    remove LaTeX build products"

venv:
	python3 -m venv .venv && .venv/bin/pip install -r simulation/requirements.txt

test:
	$(PY) -m pytest simulation/tests -q

results:
	cd simulation/studies && ../../$(PY) run_all.py && ../../$(PY) report_scaling.py

figures:
	cd simulation/plotting && ../../$(PY) make_figures.py

check-tex:
	@command -v pdflatex >/dev/null 2>&1 || { echo "pdflatex not found. Install TeX Live:"; \
	 echo "  sudo apt install texlive-latex-extra texlive-pictures texlive-publishers texlive-fonts-recommended latexmk"; \
	 exit 1; }

paper: check-tex
	cd paper/manuscript && $(TEX) main.tex

report: check-tex
	cd paper/background_report && $(TEX) report.tex
	cd paper/background_report && $(TEX) core_facts.tex

docs: check-tex
	cd simulation/docs && $(TEX) code_guide.tex

slides:
	cd presentation && ../$(PY) make_assets.py && node build_deck.js

all: figures paper report docs

clean:
	cd paper/manuscript && rm -f *.aux *.bbl *.blg *.log *.out *.fls *.fdb_latexmk main.pdf
	cd paper/background_report && rm -f *.aux *.bbl *.blg *.log *.out *.toc *.fls *.fdb_latexmk report.pdf core_facts.pdf
	cd simulation/docs && rm -f *.aux *.bbl *.blg *.log *.out *.toc *.fls *.fdb_latexmk code_guide.pdf
