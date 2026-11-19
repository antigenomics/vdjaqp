# VDJ-AQP development

VDJ-AQP is scripts, a Python package and a Nextflow workflow. Reuse the existing
Evidence APIs and seqtree batch mapping; do not create a separate mapper.

- Stable branch: `master`; integration branch: `dev`. Branch features from `dev`.
- Use `uv sync --frozen --extra test --extra docs --extra data`.
- Check `.venv/bin/python -m pytest tests -q`, the toy example and
  `make -C docs html` (warnings are errors).
- Software tests do not authorize scientific recomputation. Preserve scientific
  defaults, populations, provenance, measured sequence extent and ambiguity.
- Preserve CDR3 versus junction semantics, I/L equivalence, calibrated error,
  shared-observation competition, residue caps and deterministic ties.
- Keep peptide placements separate from credited evidence and accepted identity.
- Use one native query batch; bound CPU/memory and avoid nested oversubscription.
- Never commit private inputs, stakeholder identifiers, host paths, credentials or
  generated scientific outputs. Preserve data access and redistribution limits.
- Documentation follows Diataxis and retains stable page URLs.
