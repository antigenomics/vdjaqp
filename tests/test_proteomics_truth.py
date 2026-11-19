"""Native MS2 truth must retain coisolated contributors before spectra cleanup."""
import importlib.util
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location("smoke_proteomics", Path(__file__).parents[1] / "scripts/smoke_proteomics.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(tmp_path, *, rt="4", mz="500", two_truth=False):
    spectra = tmp_path / "spectra.mzML"
    spectra.write_text('''<mzML xmlns="urn:mzml"><run><spectrumList><spectrum id="scan=2"><cvParam accession="MS:1000511" value="2"/><scanList><scan><cvParam accession="MS:1000016" value="4" unitAccession="UO:0000010"/></scan></scanList><precursorList><precursor><selectedIonList><selectedIon><cvParam accession="MS:1000744" value="500"/></selectedIon></selectedIonList></precursor></precursorList></spectrum></spectrumList></run></mzML>''')
    truth = tmp_path / "truth.idXML"
    event = f'''<PeptideIdentification RT="{rt}" MZ="{mz}"><PeptideHit sequence="PEPTIDEK" charge="2" score="0.6" protein_refs="P1 P2"><UserParam name="MZ" value="500"/></PeptideHit><PeptideHit sequence="ANOTHERK" charge="3" score="0.4" protein_refs="P3"><UserParam name="MZ" value="600"/></PeptideHit></PeptideIdentification>'''
    truth.write_text('<IdXML><IdentificationRun><ProteinIdentification><ProteinHit id="P1" accession="a"/><ProteinHit id="P2" accession="b"/><ProteinHit id="P3" accession="c"/></ProteinIdentification>'+event*(2 if two_truth else 1)+'</IdentificationRun></IdXML>')
    return spectra, truth


def test_retains_all_native_contributors(tmp_path):
    assert hasattr(module, "link_native_truth"), "native truth linkage is not implemented"
    spectra, truth = fixture(tmp_path)
    rows, ids = module.link_native_truth(spectra, truth)
    assert ids == {"scan=2"}
    assert rows[0]["native_spectrum_reference"] == "scan=2"
    assert [c["peptide_aa"] for c in rows[0]["contributors"]] == ["PEPTIDEK", "ANOTHERK"]
    assert rows[0]["contributors"][0]["protein_accessions"] == ["a", "b"]
    assert [c["intensity_fraction"] for c in rows[0]["contributors"]] == [0.6, 0.4]


@pytest.mark.parametrize("change", [{"rt":"5"}, {"mz":"501"}, {"two_truth":True}])
def test_rejects_native_truth_metadata_or_count_mismatch(tmp_path, change):
    spectra, truth = fixture(tmp_path, **change)
    with pytest.raises(ValueError):
        module.link_native_truth(spectra, truth)


def test_native_nan_fractions_preserve_identity_without_inventing_weights(tmp_path):
    spectra,truth=fixture(tmp_path)
    truth.write_text(truth.read_text().replace('score="0.6"','score="-nan"').replace('score="0.4"','score="-nan"'))
    rows,_=module.link_native_truth(spectra,truth)
    assert rows[0]['intensity_fraction_status']=='undefined_native_nan'
    assert all(c['intensity_fraction'] is None for c in rows[0]['contributors'])
    assert rows[0]['contributors'][0]['peptide_aa']=='PEPTIDEK'
    truth.write_text(truth.read_text().replace('score="-nan"','score="inf"',1))
    with pytest.raises(ValueError):module.link_native_truth(spectra,truth)
