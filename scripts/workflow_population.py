"""Stage checksum-pinned raw constants, then call the existing population driver."""
import hashlib
import json
from pathlib import Path
import sys
import subprocess
import shutil
from tempfile import TemporaryDirectory
from urllib.request import urlopen
from smoke_metrics import main

if __name__ == '__main__':
    wrapper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    manifest=Path(__file__).resolve().parents[1]/'envs/generation/constant_candidates.json'
    with TemporaryDirectory(prefix='vdjaqp-raw-constants-') as temporary:
        supplied='--additional-reference-dir' in sys.argv
        root=Path(sys.argv[sys.argv.index('--additional-reference-dir')+1]) if supplied else Path(temporary)
        records={} if supplied else json.loads(manifest.read_text())['downloaded_records']
        for name,record in records.items():
            with urlopen(record['url'],timeout=30) as source:
                content=source.read(2_000_001)
            metadata = name in ('ENST00000390559.lookup.json','ENST00000390559.cds-features.json')
            if len(content)>2_000_000 or (not metadata and hashlib.sha256(content).hexdigest()!=record['sha256']):
                raise ValueError(f'pinned public raw reference changed: {name}')
            (root/name).write_bytes(content)
        repo=manifest.parents[2]
        subprocess.run([str(repo/'.tool-envs/generation/bin/python'),
                        str(repo/'scripts/validate_constant_candidates.py'),
                        '--reference-dir',str(root)],check=True,timeout=60)
        if not supplied:sys.argv.extend(['--additional-reference-dir',str(root)])
        import importlib.util
        recipe=json.loads((manifest.parents[2]/'envs/arda/source-manifest.json').read_text())
        actual_root=Path(importlib.util.find_spec('arda').origin).parent
        modules={name:hashlib.sha256((actual_root/name).read_bytes()).hexdigest() for name in recipe['modules']}
        if any(modules[name]!=pin['patched_sha256'] for name,pin in recipe['modules'].items()):
            raise ValueError('the pinned ARDA source overlay must be active; run envs/arda/setup_overlay.sh')
        main()
        # Evaluation-only labels are read after the complete observed scoring run.
        import polars as pl
        from vdjaqp.bench import classification_metrics
        output=Path(sys.argv[sys.argv.index('--output')+1])
        (output/'workflow_producer.json').write_text(json.dumps(dict(workflow_wrapper_sha256=wrapper_sha256,arda_distribution='2.36.0',arda_backport=recipe['backport_commit'],native_read_membership_patch_sha256=recipe['read_membership_patch']['sha256'],module_sha256=modules,cpu_budget=int(sys.argv[sys.argv.index('--cpus')+1]),read_error_seed=int(sys.argv[sys.argv.index('--read-seed')+1]) if '--read-seed' in sys.argv else 2468,native_ms_randomness='reproducible'),indent=2)+'\n')
        statistics=[]
        for mode in ('amplicon','rnaseq'):
            candidates=pl.read_parquet(output/(mode+'_candidate_identity.parquet'))
            scores=pl.read_parquet(output/(mode+'_scores.parquet'))
            positives=set(candidates.filter(pl.col('full_class_count')==1)['clone_id'])
            if positives:
                values=classification_metrics(scores,positives)
            else:
                n_calls=scores.filter(pl.col('confirmed_nats')>=20).height
                values=dict(n_candidates=scores.height,n_positives=0,
                            n_benchmark_negatives=scores.height,true_positive=0,
                            false_positive=n_calls,false_negative=0,
                            true_negative=scores.height-n_calls,n_clonotypes=scores.height,threshold_confirmed_nats=20,
                            median_rank=None,worst_rank=None,top_50=0,top_500=0,f1_at20=0. if n_calls else None,
                            precision_at20=0. if n_calls else None,recall_at20=None,
                            average_precision=None,roc_auc=None)
            for quantity,value in values.items():
                count=quantity.startswith('n_') or quantity in ('true_positive','false_positive','false_negative','true_negative','worst_rank','top_50','top_500')
                statistics.append(dict(geometry=mode,quantity=quantity,
                    value=None if value is None else str(int(value)) if count else str(value),
                    units='rank' if 'rank' in quantity else 'count' if count else 'nats' if quantity=='threshold_confirmed_nats' else 'fraction',
                    n_candidates=scores.height,n_positives=len(positives),
                    positive_definition='candidate exactly compatible with one full true heavy-variable mass class',
                    negative_definition='all remaining observed candidates for the full-variable identity endpoint; partial compatibility retained separately',
                    source=mode+'_candidate_identity.parquet + '+mode+'_scores.parquet'))
        pl.DataFrame(statistics).write_csv(output/'classification_statistics.tsv',separator='\t')

        # Retain requested final tables/proof, never native bulk intermediates.
        if output.is_symlink():raise ValueError('refuse retention through symlink output')
        receipt=json.loads((output/'metric_integration.json').read_text())
        removed=[]
        for path in list(output.iterdir()):
            keep_dir=path.name in ('amplicon_observed','rnaseq_observed') or (path.name=='references' and '--reference-dir' not in sys.argv)
            keep_file=path.suffix in ('.parquet','.tsv','.json')
            if (path.is_dir() and keep_dir) or (path.is_file() and keep_file):continue
            if path.is_dir():
                summary=path/'summary.json'
                if summary.is_file():shutil.copy2(summary,output/(path.name+'_native_summary.json'))
                removed.extend(str(p.relative_to(output)) for p in path.rglob('*') if p.is_file())
                shutil.rmtree(path)
            else:
                removed.append(path.name);path.unlink()
        (output/'retention.json').write_text(json.dumps(dict(
            removed_artifact_sha256={name:receipt['artifact_sha256'].get(name) for name in removed},
            meaning='native inputs/intermediates removed after validated final exports; raw public reference mirror remains external',
            retained_sha256={str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in output.rglob('*') if p.is_file()}
        ),indent=2)+'\n')
