"""Tiny unit tables use one immutable native pool; subprocess tests cover all CPU modes."""
import vdjaqp

# Select before test modules import/query Polars. Individual CLI tests must not
# pretend that the test interpreter can resize an already initialized native pool.
vdjaqp.resolve_cpu_budget('serial')
