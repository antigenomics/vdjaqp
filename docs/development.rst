Development and review
======================

Recreate and verify the locked environment from the repository root:

.. code-block:: sh

   uv sync --frozen --extra test --extra docs --extra data
   .venv/bin/python -m pytest tests -q
   .venv/bin/python examples/toy/run_toy.py
   make -C docs html

The tests use invented inputs. The unit-test interpreter starts with a serial
Polars pool; subprocess checks exercise other CPU modes and process cleanup.
Optional native-runtime tests skip when their documented prerequisites are absent.
A passing unit suite does not qualify a scientific dataset or every native host.

CPU allocation
--------------

Select the process budget before importing numerical libraries:

.. code-block:: python

   from vdjaqp import resolve_cpu_budget
   resolve_cpu_budget("all")  # or a smaller positive allocation
   from vdjaqp.pipeline import Evidence

Positive oversized requests are capped to the usable allocation. An initialized
Polars pool cannot shrink; requesting a smaller ceiling then raises without
changing the declared budget. Start a fresh process instead. Raising a ceiling
also cannot enlarge an existing smaller Polars pool. Set it early for utilization.

Input and model boundaries
--------------------------

``Evidence.from_frames`` scores supplied tables using an explicit ``error_model``;
it does not load datasets or fit calibration. ``Evidence.build`` loads the named
dataset and calibrates before calling the same frame route. Use the former for
independently calibrated workflow inputs; see :doc:`pipeline` and :doc:`input-output`.

Confidence calibration pools observations with identical confidence before
weighted isotonic regression. This corrects an order-dependent fit at ties and
can change future newly fitted models. Frozen model files, loading/interpolation,
scientific thresholds and accepted historical results are unchanged. No scientific
recalibration is part of software verification.

Simplification review
---------------------

The repository review covered package modules, scripts, workflow dispatch,
packaging, environment recipes and documentation. Changes favor existing APIs:

* Remove an unreferenced historical runtime-cleanup script from public software.
* Reuse streaming SHA-256 and same-run baseline scores instead of duplicating work.
* Use the existing file-lock dependency for destination publication.
* Let uv read ``.python-version`` instead of passing an unsupported action input.

Residue caps, shared-observation competition and ambiguity handling are deliberate
model behavior. Their complexity is not a reason to remove them. Scientific
recomputation and optional native-platform qualification require their own scope.

Agent guidance is maintained in ``skills/vdjaqp/SKILL.md`` at the repository root.
Branch from ``dev``, submit a reviewable change there, then promote verified
integration to ``master``. Keep required-review protections enabled.
