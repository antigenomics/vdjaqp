"""The public synthetic tutorial runs the unchanged scoring path and shows shared versus distinguishing evidence."""
import importlib.util
from pathlib import Path


def test_toy_example_separates_target_from_relatives():
    path = Path(__file__).resolve().parents[1] / 'examples/toy/run_toy.py'
    spec = importlib.util.spec_from_file_location('run_toy', path)
    toy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(toy)
    clones, peptides, evidence, table = toy.run()
    assert table.height == clones.height == 303
    rows = {r['clone_id']: r for r in table.iter_rows(named=True)}
    assert table['clone_id'][0] == 'target' and table.filter('accepted')['clone_id'].to_list() == ['target']
    assert rows['target']['junction_nats'] > 20
    # Every candidate shares framework peptides, but shared evidence is worth little.
    assert table.filter(table['peptide_matches'] > 0).height == 303
    assert rows['junction_relative']['confirmed_nats'] < 5
    assert max(r['confirmed_nats'] for k, r in rows.items() if k.startswith('background_')) < 1


def test_database_error_rate_is_the_e_db_keyword():
    """Database hits carry no residue confidence, so only ``e_db`` sets their error rate.

    The fourth positional ``ErrorModel`` argument is a substitution count; passing an error
    rate there leaves ``e_db`` at its default. Guard the keyword route used by the tutorial.
    """
    import numpy as np
    from vdjaqp.errormodel import ErrorModel, mass_prior
    from vdjaqp.pipeline import Evidence
    path = Path(__file__).resolve().parents[1] / 'examples/toy/run_toy.py'
    spec = importlib.util.spec_from_file_location('run_toy', path)
    toy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(toy)
    clones, peptides = toy.repertoire(), toy.peptides()
    target = []
    for e_db in (0.02, 0.10, 0.30):
        model = ErrorModel(np.array([0., 1.]), np.array([.1, .01]), mass_prior(), 0, e_db=e_db)
        assert model.e_db == e_db
        scores = Evidence.from_frames(clones, peptides, error_model=model).score(use_prior=False)
        target.append(scores.filter(scores['clone_id'] == 'target')['confirmed_nats'][0])
    assert target[0] > target[1] > target[2] > 20
