"""Build a tiny synthetic shared-family fixture using published GenAIRR and real C references.

Software integration fixture only: supplied abundances/SHM are not experiment parameters.
"""
import argparse
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import translate
from GenAIRR import Experiment


def reference_hashes(reference_dir, all_loci=False):
    names = ['P01857-1.fasta','P01834.fasta','J00228.embl','J00241.embl']
    if all_loci:
        names += ['P0DOY2.fasta', 'J00253.embl']
    paths = [Path(reference_dir)/name for name in names]
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def load_constants(reference_dir, all_loci=False):
    constants = {}
    specifications = [('IGH','J00228','P01857-1',330),('IGK','J00241','P01834',107)]
    if all_loci:
        specifications.append(('IGL','J00253','P0DOY2',106))
    for locus, ena, protein, length in specifications:
        record = SeqIO.read(Path(reference_dir)/f'{ena}.embl', 'embl')
        features = [f for f in record.features if f.type == 'CDS']
        assert len(features) == 1
        feature = features[0]
        nt = str(feature.extract(record.seq))
        aa = str(next(SeqIO.parse(Path(reference_dir)/f'{protein}.fasta','fasta')).seq)
        assert feature.qualifiers['codon_start'] == ['3']
        assert feature.qualifiers['translation'][0] == aa[1:]
        assert len(aa) == length
        constants[locus] = dict(nt=nt, aa=aa, cds_accession=record.id,
                               protein_accession=protein, cds_location=str(feature.location),
                               cds_sha256=hashlib.sha256(nt.encode()).hexdigest())
    return constants


def generate(output, all_loci=False):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    existing = [p.name for p in output.iterdir() if p.name != "references"]
    if existing:
        raise FileExistsError(f"generated fixture output already exists: {sorted(existing)}")
    references = reference_hashes(output/'references', all_loci)
    constants = load_constants(output/'references', all_loci)
    rows = []
    specifications = [(0,'IGH','IGHG1','human_igh','hh_s5f',123),
                      (0,'IGK','IGKC','human_igk','hkl_s5f',456)]
    if all_loci:
        specifications += [(1,'IGH','IGHG1','human_igh','hh_s5f',10123),
                           (1,'IGL','IGLC2','human_igl','hkl_s5f',10456)]
    for family, locus, gene, cfg, model, root_seed in specifications:
        experiment = Experiment.on(cfg)
        if locus == 'IGL':
            # A declared sequence-compatible J2/3–C2 fixture, not unrestricted
            # genomic lambda J/C sampling or an inferred unambiguous J gene.
            experiment = experiment.restrict_alleles(j=['IGLJ2/3-HDWH*01','IGLJ2/3-JWS7*02'])
        base = experiment.recombine().productive_only().compile()
        root = base.run(n=1, seed=root_seed)[0]
        mutation = Experiment.on(cfg).mutate(model='s5f', rate=0.01, s5f_model=model).productive_only().compile()
        seeds = [789+10000*family,790+10000*family]
        outcomes = [root]
        for seed in seeds:
            parent = outcomes[-1] if all_loci else root
            outcomes.append(mutation.simulator.run_from(parent.final_simulation(), seed, strict=True))
        for i, outcome in enumerate(outcomes):
            record = outcome.airr_record(base.refdata)
            nt, aa = record['sequence'].upper(), record['sequence_aa']
            assert record['productive'] and '*' not in aa
            # Engine annotation is authoritative; check frame-zero translation where complete.
            coding_nt = nt[:len(nt)//3*3]
            untranslated_tail = nt[len(coding_nt):]
            translated = str(translate(coding_nt, to_stop=False))
            assert translated == aa, (locus, i, 'translation mismatch')
            full_coding_nt = nt + constants[locus]['nt']
            assert len(full_coding_nt) % 3 == 0
            protein_aa = aa + constants[locus]['aa']
            assert str(translate(full_coding_nt)) == protein_aa+'*'
            parent_node = i-1 if all_loci else 0
            rows.append(dict(chain_id=f'family{family}_node{i}_{locus}',antibody_id=f'family{family}_node{i}',family_id=f'family{family}',parent_id=None if i==0 else f'family{family}_node{parent_node}_{locus}',locus=locus,variable_nt=nt,coding_variable_nt=coding_nt,untranslated_terminal_nt=untranslated_tail,variable_aa=aa,protein_aa=protein_aa,constant_accession=constants[locus]['protein_accession'],constant_cds=constants[locus]['cds_accession'],constant_nt_gene=gene,amplicon_nt=nt+constants[locus]['nt'][:29],rnaseq_nt=full_coding_nt,rna_pairs=[1000,600,400][i],secretion_factor=[1.,2.,0.5][i],airr_truth=record))
    for name, field in [('proteins.fasta','protein_aa'),('amplicon_templates.fasta','amplicon_nt'),('rnaseq_templates.fasta','rnaseq_nt')]:
        (output/name).write_text(''.join(f'>{r["chain_id"]}'+(f' [# intensity={r["rna_pairs"]*r["secretion_factor"]} #]' if name=='proteins.fasta' else '')+f'\n{r[field]}\n' for r in rows))
    (output/'read_counts.tsv').write_text(''.join(f'{r["chain_id"]}\t{2*r["rna_pairs"]}\n' for r in rows))
    (output/'truth.json').write_text(json.dumps({'purpose':'software-only synthetic fixture','seed':123,'chains':rows,'reference_sha256':references,'constant_references':constants, 'family_topology':'root-child-grandchild' if all_loci else 'root-two-children', 'lambda_scope':'native J2/3 sequence-compatible C2; gene alias remains ambiguous' if all_loci else 'not included', 'constant_nt_scope':'full spliced constant CDS extracted by Biopython from ENA annotated features; coding RNA includes terminal stop; no UTR or signal peptide model'},indent=2)+'\n')
    return rows


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--all-loci',action='store_true')
    a=p.parse_args();print('Generated',len(generate(a.output, a.all_loci)),'chains')
