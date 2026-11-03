# Common tasks. Override the interpreter with: make PY=python3 <target>
PY ?= .venv/bin/python
LATEX ?= latexmk -pdf -interaction=nonstopmode

.PHONY: help venv test results figures paper report all clean

help:
	@echo "venv     create .venv and install simulation/requirements.txt"
	@echo "test     fast regression tests (~10 s)"
	@echo "results  rerun all studies -> simulation/results (~15 min)"
	@echo "figures  regenerate paper figures from results (seconds)"
	@echo "paper    build paper/manuscript/main.pdf"
	@echo "report   build paper/background_report/report.pdf"
	@echo "clean    remove LaTeX build products"

venv:
	python3 -m venv .venv && .venv/bin/pip install -r simulation/requirements.txt

test:
	$(PY) -m pytest simulation/tests -q

results:
	cd simulation/studies && ../../$(PY) run_all.py && ../../$(PY) report_scaling.py

figures:
	cd simulation/plotting && ../../$(PY) make_figures.py

paper:
	cd paper/manuscript && $(LATEX) main.tex

report:
	cd paper/background_report && $(LATEX) report.tex

all: figures paper report

clean:
	cd paper/manuscript && latexmk -C main.tex || true
	cd paper/background_report && latexmk -C report.tex || true
