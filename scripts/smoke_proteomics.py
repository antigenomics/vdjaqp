"""Purpose: verify published OpenMS proteomics plumbing on a small synthetic fixture.

Context: short acquisition settings are software checks, not manuscript parameters.
Input FASTA is the simulated physical protein mixture; search FASTA must come from
observed repertoire sequences for an end-to-end check. No private cohort is read.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import os
import resource
import subprocess
import time
import xml.etree.ElementTree as ET


def link_native_truth(spectra, truth):
    """Link published out_id MS2-order truth to actual native IDs, never guesses.

    OpenMS 3.0 getIdentifications emits one row per MS2 in experiment order;
    each row contains all coisolated precursor hits with intensity fractions.
    Validate count, RT and first-precursor m/z before retaining this association.
    """
    ids, events = set(), []
    for _, elem in ET.iterparse(spectra, events=("end",)):
        if elem.tag.rsplit("}", 1)[-1] != "spectrum":
            continue
        native_id = elem.attrib["id"]
        if native_id in ids:
            raise ValueError("duplicate native spectrum ID")
        ids.add(native_id)
        cv = [e for e in elem.iter() if e.tag.rsplit("}", 1)[-1] == "cvParam"]
        def value(accession):
            return next(e for e in cv if e.get("accession") == accession)
        if value("MS:1000511").get("value") == "2":
            rt = value("MS:1000016")
            scale = {"UO:0000010": 1., "UO:0000031": 60.}.get(rt.get("unitAccession"))
            if scale is None:
                raise ValueError("unknown native RT unit")
            events.append({"native_spectrum_reference": native_id,
                           "rt_seconds": float(rt.get("value")) * scale,
                           "first_precursor_mz": float(value("MS:1000744").get("value"))})
        elem.clear()
    tree = ET.parse(truth)
    proteins = {e.attrib["id"]: e.attrib["accession"] for e in tree.iter("ProteinHit")}
    identifications = list(tree.iter("PeptideIdentification"))
    if len(identifications) != len(events):
        raise ValueError("native MS2/truth count mismatch")
    for ordinal, (event, identification) in enumerate(zip(events, identifications)):
        for key, attribute in (("rt_seconds", "RT"), ("first_precursor_mz", "MZ")):
            if not math.isclose(event[key], float(identification.attrib[attribute]), rel_tol=1e-10, abs_tol=1e-6):
                raise ValueError("native MS2/truth metadata mismatch")
        event["truth_ordinal"] = ordinal
        event["contributors"] = []
        for hit in identification.findall("PeptideHit"):
            fraction = float(hit.attrib["score"])
            if not math.isnan(fraction) and (not math.isfinite(fraction) or not 0 <= fraction <= 1):
                raise ValueError("invalid native precursor intensity fraction")
            mz = next(e.get("value") for e in hit.findall("UserParam") if e.get("name") == "MZ")
            event["contributors"].append({
                "peptide_aa": hit.attrib["sequence"], "charge": int(hit.attrib["charge"]),
                "precursor_mz": float(mz), "intensity_fraction": None if math.isnan(fraction) else fraction,
                "protein_accessions": sorted({proteins[x] for x in hit.get("protein_refs", "").split()}),
                "native_hit_attributes": dict(hit.attrib),
                "native_hit_user_params": [dict(e.attrib) for e in hit.findall("UserParam")],
            })
        fractions = [c['intensity_fraction'] for c in event['contributors']]
        # OpenMS 3 getMS2Identifications divides by the total without a zero
        # guard. Preserve all-NaN native identities, never invent their weights.
        undefined = bool(fractions) and all(f is None for f in fractions)
        if not undefined and (not fractions or any(f is None for f in fractions) or not math.isclose(sum(fractions), 1., abs_tol=1e-6)):
            raise ValueError("native contributor fractions do not sum to one")
        event['intensity_fraction_status'] = 'undefined_native_nan' if undefined else 'defined'
    return events, ids


def executable_identity(path):
    """Measure executable bytes without inferring a package or platform identity."""
    path = Path(path).resolve()
    with path.open('rb') as stream:
        return dict(path=str(path), sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("proteins", type=Path)
    parser.add_argument("database", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--keep-spectra", action="store_true", help="retain until caller validates downstream linkage")
    parser.add_argument("--enzyme", choices=("Trypsin", "Chymotrypsin"), default="Trypsin")
    parser.add_argument("--gradient-seconds",type=float,default=120.)
    parser.add_argument("--ms2-per-bin",type=int,default=5)
    parser.add_argument("--min-detectability",type=float,default=.5,
                        help="published peptide detectability-score cutoff, in [0,1]")
    parser.add_argument("--dynamic-exclusion-seconds",type=float,default=0.,
                        help="native precursor exclusion duration in seconds; zero disables exclusion")
    parser.add_argument("--mssimulator",type=Path,help="explicit MSSimulator executable; other tools use pinned environment")
    parser.add_argument("--stage-timeout",type=float,default=300.)
    args = parser.parse_args()
    if not math.isfinite(args.min_detectability) or not 0<=args.min_detectability<=1:
        parser.error("--min-detectability must be finite and in [0,1]")
    if not math.isfinite(args.dynamic_exclusion_seconds) or args.dynamic_exclusion_seconds<0:
        parser.error("--dynamic-exclusion-seconds must be finite and nonnegative")
    if not math.isfinite(args.stage_timeout) or args.stage_timeout<=0:
        parser.error("--stage-timeout must be finite and positive")
    if args.mssimulator is not None:
        args.mssimulator=args.mssimulator.resolve()
        if not args.mssimulator.is_file() or not os.access(args.mssimulator,os.X_OK):
            parser.error("--mssimulator must be an executable file")
    if args.threads < 1:
        parser.error("threads must be positive")
    if not math.isfinite(args.gradient_seconds) or args.gradient_seconds<=0 or args.ms2_per_bin<1:
        parser.error("finite positive gradient and positive MS2 capacity are required")
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0,str(root))
    from vdjaqp import resolve_cpu_budget
    args.threads=resolve_cpu_budget(args.threads)
    binaries = root / ".tool-envs/openms/bin"
    selected_simulator=args.mssimulator or binaries/"MSSimulator-digestfix"
    if args.mssimulator is None and not selected_simulator.is_file():
        raise FileNotFoundError("build the pinned simulator first: .venv/bin/python scripts/build_openms_digest_fix.py")
    simulator_hash=executable_identity(selected_simulator)['sha256']
    artifact_path=Path(str(selected_simulator)+".manifest.json")
    simulator_artifact=None
    if args.mssimulator is None and not artifact_path.is_file():
        raise FileNotFoundError("pinned simulator build manifest is required")
    if artifact_path.is_file():
        artifact_bytes=artifact_path.read_bytes()
        simulator_artifact=dict(path=str(artifact_path),sha256=hashlib.sha256(artifact_bytes).hexdigest(),
                                contents=json.loads(artifact_bytes))
        if simulator_artifact['contents'].get('executable_sha256')!=simulator_hash:
            raise ValueError('selected simulator does not match its build manifest')
    args.output.mkdir(parents=True, exist_ok=False)
    out = args.output.resolve()
    stages = []
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads))

    def run(name, command):
        begin = time.perf_counter()
        with (out / (name + ".log")).open("w") as log:
            executable=selected_simulator if command[0]=="MSSimulator" else binaries/command[0]
            subprocess.run([str(executable), *command[1:]], env=env,
                           stdout=log, stderr=subprocess.STDOUT, check=True, timeout=args.stage_timeout)
        stages.append({"stage": name, "seconds": time.perf_counter() - begin})

    ini = out / "mssimulator.ini"
    run("ini", ["MSSimulator", "-write_ini", str(ini)])
    tree = ET.parse(ini)
    settings = {
        "algorithm:MSSim:Digestion:enzyme": args.enzyme,
        "algorithm:MSSim:RT:total_gradient_time": str(args.gradient_seconds),
        "algorithm:MSSim:RT:scan_window:min": "0.0",
        "algorithm:MSSim:RT:scan_window:max": str(args.gradient_seconds),
        "algorithm:MSSim:RT:sampling_rate": "2.0",
        "algorithm:MSSim:RawSignal:enabled": "true",
        "algorithm:MSSim:RawTandemSignal:status": "precursor",
        "algorithm:MSSim:RawTandemSignal:tandem_mode": "0",
        "algorithm:MSSim:RawTandemSignal:Precursor:ms2_spectra_per_rt_bin": str(args.ms2_per_bin),
        "algorithm:MSSim:Detectability:dt_simulation_on": "true",
        "algorithm:MSSim:Detectability:min_detect": str(args.min_detectability),
        "algorithm:MSSim:RawTandemSignal:Precursor:Exclusion:use_dynamic_exclusion":
            "true" if args.dynamic_exclusion_seconds>0 else "false",
        "algorithm:RandomNumberGenerators:biological": "reproducible",
        "algorithm:RandomNumberGenerators:technical": "reproducible",
    }
    if args.dynamic_exclusion_seconds>0:
        settings["algorithm:MSSim:RawTandemSignal:Precursor:Exclusion:exclusion_time"]=str(args.dynamic_exclusion_seconds)
    found = set()
    def visit(node, parts=()):
        for child in node:
            path = (*parts, child.get("name", ""))
            key = ":".join(path[2:])
            if child.tag == "ITEM" and key in settings:
                child.set("value", settings[key]); found.add(key)
            visit(child, path)
    visit(tree.getroot())
    assert found == set(settings), set(settings) - found
    tree.write(ini, encoding="utf-8", xml_declaration=True)
    spectra, ids = out / "spectra.mzML", out / "psms.idXML"
    run("simulate", ["MSSimulator", "-ini", str(ini), "-in", str(args.proteins.resolve()),
        "-out", str(spectra), "-out_fm", str(out / "truth.featureXML"), "-out_id", str(out / "truth.idXML"), "-threads", str(args.threads)])
    search_ini = out / "search.ini"
    run("search_ini", ["SimpleSearchEngine", "-write_ini", str(search_ini)])
    search_tree = ET.parse(search_ini)
    for node in search_tree.iter("ITEMLIST"):
        if node.get("name") in {"fixed", "variable"}:
            for item in list(node): node.remove(item)
    search_tree.write(search_ini, encoding="utf-8", xml_declaration=True)
    run("search", ["SimpleSearchEngine", "-ini", str(search_ini), "-in", str(spectra), "-database", str(args.database.resolve()),
        "-out", str(ids),
        "-Search:enzyme", args.enzyme,
        "-threads", str(args.threads)])
    run("export", ["TextExporter", "-in", str(ids), "-out", str(out / "peptides.tsv"), "-id:peptides_only"])
    truth_rows, observed = link_native_truth(spectra, out / "truth.idXML")
    (out / "ms2_truth.json").write_text(json.dumps(truth_rows, indent=2) + "\n")
    references = [e.get("spectrum_reference") for e in ET.parse(ids).iter("PeptideIdentification")]
    assert set(references) <= observed, "PSM spectrum references do not link to native spectra"
    record = {"stages": stages, "psm_groups": len(references), "native_spectra": len(observed),
              "settings": settings, "threads": args.threads, "enzyme": args.enzyme,
              "min_detectability":args.min_detectability,
              "dynamic_exclusion_seconds":args.dynamic_exclusion_seconds,
              "native_ms2_truth_groups": len(truth_rows),
              "source_inspection_pin": "OpenMS Release3.0.0 / 60ce8966d577b43a6419403dec171539e08fe513",
              "search_export_executables": {name: executable_identity(binaries/name)
                                            for name in ("SimpleSearchEngine", "TextExporter")},
              "selected_mssimulator":str(selected_simulator),"selected_mssimulator_sha256":simulator_hash,
              "selected_mssimulator_artifact":simulator_artifact,
              "mssimulator_override":args.mssimulator is not None,"stage_timeout_seconds":args.stage_timeout,
              "peak_child_rss_bytes": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * (1 if sys.platform == "darwin" else 1024)}
    (out / "summary.json").write_text(json.dumps(record, indent=2) + "\n")
    if not args.keep_spectra:
        spectra.unlink()
    print(json.dumps(record))


if __name__ == "__main__":
    main()
