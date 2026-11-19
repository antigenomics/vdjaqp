PYTHON ?= .venv/bin/python
CPUS ?= all

.PHONY: analyze test
analyze:
	$(PYTHON) scripts/recount_publication.py --cpus $(CPUS)

test:
	$(PYTHON) -m pytest tests -q
