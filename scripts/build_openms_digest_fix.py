"""Build a scoped native OpenMS metadata fix; preserve installed MSSimulator.

Only raw pinned upstream sources are retained; every invocation recompiles all
three translation units and rebuilds the technical shared-peptide fixture.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import struct
import subprocess
import tempfile
from urllib.request import urlopen

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO/'envs/openms/native_digest_fix.json'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def patch_digest(source):
    """Batch native accession metadata; leave biological/numerical operations intact."""
    declaration = '    std::map<AASequence, Feature> generated_features;'
    begin = '        // add current protein accession\n'
    finish = '        generated_features[dp].setPeptideIdentifications(pep_idents);\n'
    final = '      // round up intensity\n'
    if any(source.count(site) != 1 for site in (declaration, begin, finish, final)):
        raise ValueError('pinned native accession patch sites differ')
    start = source.index(begin)
    end = source.index(finish, start) + len(finish)
    source = source[:start] + '        generated_accessions[dp].insert(protein_hit.getAccession());\n' + source[end:]
    source = source.replace(declaration, declaration + '\n    std::map<AASequence, std::set<String>> generated_accessions;')
    materialize = """      std::vector<PeptideIdentification> pep_idents = it_gf->second.getPeptideIdentifications();
      std::vector<PeptideHit> pep_hits = pep_idents[0].getHits();
      for (const String& accession : generated_accessions.at(it_gf->first))
      {
        PeptideEvidence pe;
        pe.setProteinAccession(accession);
        pep_hits[0].addPeptideEvidence(pe);
      }
      pep_idents[0].setHits(pep_hits);
      it_gf->second.setPeptideIdentifications(pep_idents);

"""
    return source.replace(final, materialize + final)


def build_flags(prefix):
    """Validate installed ABI and tools before creating any build directory."""
    system = platform.system()
    if system not in ('Darwin', 'Linux'):
        raise ValueError('scoped native build supports Darwin and Linux only')
    library = prefix/'lib'/('libOpenMS.dylib' if system == 'Darwin' else 'libOpenMS.so')
    for path in (prefix/'bin/MSSimulator', library,
                 prefix/'include/OpenMS/SIMULATION/DigestSimulation.h',
                 prefix/'include/qt/QtCore/qstring.h'):
        if not path.is_file():
            raise ValueError(f'required installed OpenMS input missing: {path}')
    compiler = shutil.which('clang++', path='/usr/bin') if system == 'Darwin' else shutil.which('g++')
    if compiler is None:
        raise ValueError('required native C++ compiler unavailable')
    if system == 'Darwin':
        with library.open('rb') as handle:
            header = handle.read(8)
        if header != bytes.fromhex('cffaedfe07000001'):
            raise ValueError('Darwin build requires installed thin x86_64 OpenMS ABI')
        flags = [compiler, '-arch', 'x86_64', '-std=c++17', '-O2', '-mmacosx-version-min=11.0']
    else:
        if platform.machine() != 'x86_64':
            raise ValueError('Linux build requires x86_64 host')
        for path in (prefix/'bin/MSSimulator', library):
            with path.open('rb') as handle:
                header = handle.read(20)
            if (len(header) != 20 or header[:6] != b'\x7fELF\x02\x01'
                    or struct.unpack_from('<H', header, 18)[0] != 62):
                raise ValueError('Linux build requires installed ELF64 little-endian x86_64 ABI')
        nm = shutil.which('nm')
        if nm is None:
            raise ValueError('nm is required to verify installed Linux C++ ABI')
        symbols = subprocess.check_output([nm, '-D', '--defined-only', str(library)], text=True, timeout=30)
        if '_ZN6OpenMS6StringC1ERKNSt7__cxx1112basic_string' not in symbols:
            raise ValueError('Linux build requires verified libstdc++ C++11 string ABI')
        flags = [compiler, '-std=c++17', '-O2', '-fPIC', '-D_GLIBCXX_USE_CXX11_ABI=1']
    flags += ['-I'+str(prefix/'include'/part) for part in
              ('', 'qt', 'qt/QtCore', 'qt/QtNetwork', 'eigen3', 'OpenMS/thirdparty/evergreen')]
    link = ['-L'+str(prefix/'lib'), '-Wl,-rpath,'+str(prefix/'lib'), '-lOpenMS', '-lQt5Core', '-lQt5Network']
    if system == 'Linux':
        link += ['-lxerces-c']  # MSSimulator TU directly references SAXParseException RTTI.
    return compiler, flags, link, '*.dylib' if system == 'Darwin' else '*.so*'


def build(prefix):
    prefix = Path(prefix).resolve()
    compiler, flags, link, library_pattern = build_flags(prefix)
    manifest = json.loads(MANIFEST.read_text())
    original_hash = sha256(prefix/'bin/MSSimulator')
    source_dir = prefix/'native_digest_fix/source'
    source_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, record in manifest['sources'].items():
        path = source_dir/name
        # Only original raw upstream source downloads may be reused.
        if not path.exists():
            with urlopen(record['url'], timeout=30) as response:
                path.write_bytes(response.read())
        if sha256(path) != record['sha256']:
            raise ValueError(f'pinned source SHA256 mismatch: {name}')
        paths.append(path)
    commands = []
    def run(command):
        commands.append(list(map(str, command)))
        return subprocess.run(command, check=True, text=True, capture_output=True, timeout=120)
    with tempfile.TemporaryDirectory(prefix='openms-digest-build-', dir=source_dir.parent) as temp:
        work = Path(temp)
        patched = work/'DigestSimulation.cpp'
        patched.write_text(patch_digest(paths[0].read_text()))
        binary = prefix/'bin/MSSimulator-digestfix'
        staged = work/'MSSimulator-digestfix'
        command = flags + [str(patched), str(paths[1]), str(paths[2])] + link + ['-o', str(staged)]
        try:
            run(command)
            fixture = REPO/'tests/openms_digest_fixture.cpp'
            results = {}
            for mode, sources in [('original', [str(fixture)]), ('patched', [str(fixture), str(patched)])]:
                exe = work/mode
                run(flags + sources + link + ['-o', str(exe)])
                result = run([str(exe), '6' if mode == 'original' else '3'])
                results[mode] = result.stdout.strip()
            run([str(staged), '--help'])
        except subprocess.CalledProcessError as error:
            raise RuntimeError(error.stderr or error.stdout) from error
        if sha256(prefix/'bin/MSSimulator') != original_hash:
            raise ValueError('original MSSimulator changed during scoped build')
        staged.replace(binary)
        report = dict(original_executable_sha256=original_hash, upstream=manifest, patched_digest_sha256=sha256(patched),
            compiler=subprocess.check_output([compiler, '--version'], text=True).strip(),
            compiler_sha256=sha256(compiler), executable=str(binary), executable_sha256=sha256(binary),
            linked_library_sha256={p.name:sha256(p) for p in (prefix/'lib').glob(library_pattern) if not p.is_symlink()},
            fixture_sha256=sha256(fixture), fixture_results=results, commands=commands)
        target = binary.with_name(binary.name+'.manifest.json')
        target.write_text(json.dumps(report, indent=2)+'\n')
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, default=REPO/'.tool-envs/openms')
    print(build(parser.parse_args().prefix))
