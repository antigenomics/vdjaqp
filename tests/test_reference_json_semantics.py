"""JSON byte drift cannot relax transcript/CDS or sequence integrity pins."""
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory


def test_json_invariants(tmp_path):
    if __name__ != '__main__':
        import pytest
        pytest.importorskip('Bio', reason='pinned generation runtime required')
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
    from validate_constant_candidates import validated_source_manifest
    row=dict(gene='IGHM',nucleotide_accession='ENST00000390559.6',assembly='GRCh38',
             cds_features=[dict(Parent='ENST00000390559',start=2,end=5,strand=-1)])
    meta=dict(id='ENST00000390559',version=6,assembly_name='GRCh38',strand=-1,
              seq_region_name='14',source='havana')
    names=['ENST00000390559.lookup.json','ENST00000390559.cds-features.json','sequence.fasta']
    contents=[json.dumps(meta).encode(),json.dumps(row['cds_features']).encode(),b'>raw\nACGT\n']
    records={}
    for name,content in zip(names,contents):
        (tmp_path/name).write_bytes(content)
        records[name]=dict(sha256=hashlib.sha256(content).hexdigest())
    manifest=tmp_path/'manifest.json'
    manifest.write_text(json.dumps(dict(candidates=[row],downloaded_records=records)))
    (tmp_path/names[0]).write_text(json.dumps(meta,sort_keys=True,indent=2))
    validated=validated_source_manifest(tmp_path,manifest)
    assert validated['downloaded_records'][names[0]]['original_sha256']==records[names[0]]['sha256']
    def reject():
        try:validated_source_manifest(tmp_path,manifest)
        except ValueError:return
        raise AssertionError('changed biological source accepted')
    meta['version']=7
    (tmp_path/names[0]).write_text(json.dumps(meta));reject()
    meta['version']=6
    (tmp_path/names[0]).write_text(json.dumps(meta))
    changed=[dict(row['cds_features'][0],end=6)]
    (tmp_path/names[1]).write_text(json.dumps(changed));reject()
    (tmp_path/names[1]).write_bytes(contents[1])
    (tmp_path/names[2]).write_bytes(b'>raw\nACTT\n');reject()


if __name__=='__main__':
    with TemporaryDirectory() as temporary:test_json_invariants(Path(temporary))
    print('JSON order accepted; changed version/CDS/sequence rejected')
