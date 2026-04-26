.PHONY: all install test coverage docs reference-vectors \
        c-build c-test cpp-build cpp-test \
        fortran-build fortran-test julia-test rust-test \
        js-test octave-test \
        validate figures clean sync-git help

PYTHON = python3
PIP    = pip3

help:
	@echo "UniversalFFT build targets:"
	@echo "  install           Install Python package (editable)"
	@echo "  test              Run Python unit tests"
	@echo "  coverage          Run tests with HTML coverage report"
	@echo "  reference-vectors Generate cross-language reference .npz files"
	@echo "  c-build           Build C library and demo"
	@echo "  c-test            Build and run C demo"
	@echo "  cpp-build         Build C++ demo"
	@echo "  cpp-test          Build and run C++ demo"
	@echo "  fortran-build     Build Fortran demo"
	@echo "  fortran-test      Build and run Fortran demo"
	@echo "  julia-test        Run Julia tests and demo"
	@echo "  rust-test         Run Rust tests and demo"
	@echo "  js-test           Run JavaScript tests and demo (Node ≥ 18)"
	@echo "  octave-test       Run Octave demo"
	@echo "  validate          Run cross-language validation (after demos)"
	@echo "  figures           Generate documentation figures (PNG)"
	@echo "  docs              Build MkDocs HTML documentation"
	@echo "  clean             Remove build artefacts"
	@echo "  sync-git          Stage, commit, and push to origin/master"
	@echo "                    Optional: MSG=\"commit message\""

install:
	$(PIP) install -e "python/.[dev]"

test: install
	pytest tests/python/ -v --tb=short

coverage: install
	pytest tests/python/ \
	  --cov=python/universalfft \
	  --cov-report=term-missing \
	  --cov-report=html:htmlcov \
	  --cov-fail-under=90
	@echo "HTML report: htmlcov/index.html"

reference-vectors: install
	$(PYTHON) tests/cross_language/generate_reference.py

# ── C ─────────────────────────────────────────────────────────────────────
c-build:
	$(MAKE) -C c all

c-test: c-build
	mkdir -p tests/data
	./c/ufft_demo tests/data

# ── C++ ───────────────────────────────────────────────────────────────────
cpp-build:
	$(MAKE) -C cpp all

cpp-test: cpp-build
	$(MAKE) -C cpp test

# ── Fortran ───────────────────────────────────────────────────────────────
fortran-build:
	$(MAKE) -C fortran all

fortran-test: fortran-build
	$(MAKE) -C fortran test

# ── Julia ─────────────────────────────────────────────────────────────────
julia-test:
	cd julia && julia --project=. -e 'import Pkg; Pkg.instantiate()'
	cd julia && julia --project=. test/runtests.jl
	mkdir -p tests/data
	cd julia && julia --project=. demo.jl ../tests/data

# ── Rust ──────────────────────────────────────────────────────────────────
rust-test:
	cd rust && cargo test
	mkdir -p tests/data
	cd rust && cargo run -- ../tests/data

# ── JavaScript ────────────────────────────────────────────────────────────
js-test:
	cd js && node --test test.js
	mkdir -p tests/data
	cd js && node demo.js ../tests/data

# ── Octave ────────────────────────────────────────────────────────────────
octave-test:
	mkdir -p tests/data
	cd octave && octave --no-gui ufft_demo_octave.m ../tests/data

# ── Validation ────────────────────────────────────────────────────────────
validate: install
	$(PYTHON) tests/cross_language/validate_all.py

figures: install
	$(PYTHON) docs/gen_figures.py

docs: install figures
	pip install -q -r docs/requirements.txt
	mkdocs build --strict

clean:
	$(MAKE) -C c clean
	$(MAKE) -C cpp clean
	$(MAKE) -C fortran clean
	rm -rf htmlcov .coverage coverage.xml
	rm -rf docs/_build site/
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true

# ── Git sync ──────────────────────────────────────────────────────────────
## Usage:
##   make sync-git                     # commit everything with auto message
##   make sync-git MSG="your message"  # commit with a custom message
sync-git:
	git add -A
	git diff --cached --quiet || git commit -m "$(if $(MSG),$(MSG),chore: sync local changes [$(shell date '+%Y-%m-%d %H:%M')])"
	git push origin master
