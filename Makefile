.PHONY: all test clean figure install supp docs
PY := python3
install: ; $(PY) -m pip install -r requirements.txt
all: results/01_meta_analysis.txt results/02_sanity_checks.txt \
     results/03_ledger_bias.txt results/04_diagnostic_utility.txt \
     results/06_extended_analyses.txt supp figure
supp: src/07_emit_supp_tables.py src/eaat/*.py data/*.csv
	@mkdir -p results && $(PY) src/07_emit_supp_tables.py > results/supp_S7_S12.md
	@echo 'wrote results/supp_S7_S12.md (paste verbatim into the supplementary)'
results/%.txt: src/%.py src/eaat/*.py data/*.csv
	@mkdir -p results && $(PY) $< | tee $@
figure: src/05_figure1.py src/08_figure2.py src/09_figure3_prisma.py src/eaat/*.py
	@mkdir -p results && $(PY) src/05_figure1.py && $(PY) src/08_figure2.py && $(PY) src/09_figure3_prisma.py
test: ; $(PY) -m pytest tests/ -q
clean: ; rm -rf results/*.txt results/Figure1.* __pycache__ src/__pycache__ \
	src/eaat/__pycache__ tests/__pycache__ .pytest_cache

docs: figure
	cd manuscript && $(PY) build_manuscript.py --build manuscript_src.md refs.json Manuscript_submission.md
	cd manuscript && node build_docx.js
