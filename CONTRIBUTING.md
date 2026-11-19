# Contributing

Branch from `dev` with a descriptive feature or fix branch. Submit changes to
`dev`; promote a verified integration to stable `master` through a pull request.

Recreate the pinned environment with
`uv sync --frozen --extra test --extra docs --extra data`. Run
`.venv/bin/python -m pytest tests -q`,
`.venv/bin/python examples/toy/run_toy.py` and `make -C docs html`.
CI runs these data-free checks on pushes and pull requests. Optional native-tool
tests require the documented isolated runtimes. Do not run scientific datasets
or change scientific defaults as incidental software verification.

Use the existing scripts and workflow interfaces. Retain schema, confidence,
missingness and receptor-boundary semantics, and document behavior changes.
Keep private data and host configuration outside the source tree.
