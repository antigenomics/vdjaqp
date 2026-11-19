"""Pinned native patch contract; ordinary tests never rebuild or simulate MS."""
import importlib.util
import json
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('digest_build', REPO/'scripts/build_openms_digest_fix.py')
build = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build)


def test_patch_batches_only_accession_metadata():
    source = """    std::map<AASequence, Feature> generated_features;
unchanged intensity accumulation
        // add current protein accession
old per-contribution metadata rebuilding
        generated_features[dp].setPeptideIdentifications(pep_idents);
unchanged feature iteration
      // round up intensity
unchanged rounding and output
"""
    patched = build.patch_digest(source)
    assert 'generated_accessions[dp].insert(protein_hit.getAccession());' in patched
    assert 'old per-contribution' not in patched
    assert 'std::map<AASequence, std::set<String>> generated_accessions;' in patched
    assert 'generated_accessions.at(it_gf->first)' in patched
    for unchanged in ('unchanged intensity accumulation', 'unchanged feature iteration',
                      'unchanged rounding and output', '// round up intensity'):
        assert patched.count(unchanged) == source.count(unchanged) == 1
    for invalid in ('missing', source+source):
        with pytest.raises(ValueError):
            build.patch_digest(invalid)


def test_built_fixture_manifest_matches_current_native_artifacts():
    binary = REPO/'.tool-envs/openms/bin/MSSimulator-digestfix'
    path = binary.with_name(binary.name+'.manifest.json')
    if not path.exists():
        pytest.skip('optional scoped native build absent')
    report = json.loads(path.read_text())
    assert report['executable_sha256'] == build.sha256(binary)
    assert report['fixture_sha256'] == build.sha256(REPO/'tests/openms_digest_fixture.cpp')
    original = next(line for line in report['fixture_results']['original'].splitlines() if line.startswith('PEPTIDER intensity='))
    patched = next(line for line in report['fixture_results']['patched'].splitlines() if line.startswith('PEPTIDER intensity='))
    for result in report['fixture_results'].values():
        assert 'overlap_order_checks=passed' in result
    assert original == 'PEPTIDER intensity=3000 unique_accessions=3 evidence_rows=6'
    assert patched == 'PEPTIDER intensity=3000 unique_accessions=3 evidence_rows=3'
    assert report['upstream'] == json.loads(build.MANIFEST.read_text())
    source = REPO/'.tool-envs/openms/native_digest_fix/source/DigestSimulation.cpp'
    import hashlib
    assert hashlib.sha256(build.patch_digest(source.read_text()).encode()).hexdigest() == report['patched_digest_sha256']


def test_linux_flags_and_prewrite_abi_validation(tmp_path, monkeypatch):
    import struct
    monkeypatch.setattr(build.platform, 'system', lambda:'Linux')
    monkeypatch.setattr(build.platform, 'machine', lambda:'x86_64')
    monkeypatch.setattr(build.shutil, 'which', lambda name, **kw:'/usr/bin/'+name)
    with pytest.raises(ValueError, match='missing'):
        build.build(tmp_path)
    assert not (tmp_path/'native_digest_fix').exists()
    elf = bytearray(20)
    elf[:6] = b'\x7fELF\x02\x01'
    struct.pack_into('<H', elf, 18, 62)
    for relative in ('bin/MSSimulator','lib/libOpenMS.so',
                     'include/OpenMS/SIMULATION/DigestSimulation.h','include/qt/QtCore/qstring.h'):
        path = tmp_path/relative
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(elf)
    monkeypatch.setattr(build.subprocess, 'check_output', lambda *a, **kw:
        '_ZN6OpenMS6StringC1ERKNSt7__cxx1112basic_string')
    compiler, flags, link, pattern = build.build_flags(tmp_path)
    assert compiler == '/usr/bin/g++' and pattern == '*.so*'
    assert '-fPIC' in flags and '-std=c++17' in flags
    assert '-D_GLIBCXX_USE_CXX11_ABI=1' in flags
    assert '-arch' not in flags and not any('mmacosx' in flag for flag in flags)
    assert '-lOpenMS' in link and '-lxerces-c' in link
    assert link.index('-lxerces-c') > link.index('-lOpenMS')
    monkeypatch.setattr(build.subprocess, 'check_output', lambda *a, **kw:'unsupported old ABI')
    with pytest.raises(ValueError, match='string ABI'):
        build.build(tmp_path)
    assert not (tmp_path/'native_digest_fix').exists()

    monkeypatch.setattr(build.platform, 'system', lambda:'Darwin')
    (tmp_path/'lib/libOpenMS.dylib').write_bytes(bytes.fromhex('cffaedfe07000001'))
    _, _, mac_link, _ = build.build_flags(tmp_path)
    assert '-lxerces-c' not in mac_link
