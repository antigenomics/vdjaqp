"""Adapt native OpenMS idXML PSM metadata without inventing confidence values."""
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import polars as pl

_SCHEMA={'peptide_id':pl.String,'spectrum_id':pl.String,'native_spectrum_reference':pl.String,
    'seq':pl.String,'local_confidence':pl.List(pl.Float32),'search_type':pl.String,
    'alc':pl.Float64,'alc_rank':pl.Float64,'n_reads':pl.Int64,'charge':pl.Int64,
    'rt':pl.Float64,'mz':pl.Float64,'hyperscore':pl.Float64}


def read_idxml(path: str | Path, *, run_id: str) -> pl.DataFrame:
    """Read all native PSM hits; spectrum identity and scores remain distinct.

    Supports the unmodified peptide sequences used by the bounded fixture;
    modified sequences fail explicitly until a published modification parser
    is configured. No PSM threshold or posterior-probability fit is added.
    """
    if not run_id.strip():
        raise ValueError('run_id must identify the actual simulated acquisition')
    rows=[]
    for _, node in ET.iterparse(path,events=('end',)):
        if node.tag!='PeptideIdentification':
            continue
        metadata={p.attrib.get('name'):p.attrib.get('value') for p in node.findall('UserParam')}
        native=metadata.get('spectrum_reference') or node.attrib.get('spectrum_reference')
        if not native:
            raise ValueError('PSM lacks native spectrum_reference')
        if node.attrib.get('score_type','').lower()!='hyperscore':
            raise ValueError('expected native SimpleSearchEngine hyperscore')
        for hit in node.findall('PeptideHit'):
            seq=hit.attrib['sequence']
            if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+',seq):
                raise ValueError('modified/noncanonical peptide requires a published sequence parser')
            rows.append({'peptide_id':f'{run_id}:psm:{len(rows)}','spectrum_id':f'{run_id}:{native}',
                'native_spectrum_reference':native,'seq':seq,'local_confidence':[float('nan')]*len(seq),
                'search_type':'db','alc':float('nan'),'alc_rank':1.,'n_reads':1,
                'charge':int(hit.attrib['charge']),'rt':float(node.attrib['RT']),
                'mz':float(node.attrib['MZ']),'hyperscore':float(hit.attrib['score'])})
        node.clear()
    return pl.DataFrame(rows,schema=_SCHEMA)
