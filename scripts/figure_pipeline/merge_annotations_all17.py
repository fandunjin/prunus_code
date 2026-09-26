#!/usr/bin/env python3
"""merge_annotations_all17.py - OG-level annotation from all 17 tracks.

Sources (0903/annotation/functional/):
  * pzi, psp, ppa        - annotated on original gene IDs; the six frozen
                           haplotype tracks are joined through id_maps
  * pcam pav pspe pmal par pmu pdu psal ppe gtr rch
                         - annotated directly on the frozen stable IDs

Outputs (0903/figures/inflorescence_evolution_review/data/):
  og_annotation_all17.tsv    og -> symbol/description/GO/source tracks/TE flag
  clade_family_scan_v3.tsv   clade_family_scan.tsv + the all-17 annotation
"""
import csv
import re
from collections import Counter
from pathlib import Path

_CANDIDATES = [Path(r'${PROJECT_ROOT}'), Path('${PROJECT_ROOT}')]
B = next(p for p in _CANDIDATES if p.exists())
FIG = B / '0903/figures/inflorescence_evolution_review'
ANN = B / '0903/annotation/functional'
FROZEN = B / '0903/analysisphylogenetics/01_frozen_inputs/main'
ORTHO = B / '0903/analysisphylogenetics/02_orthofinder/main/Results_Sep14/Orthogroups'
OUT = FIG / 'data'

TRACKS = ['pzi_A', 'pzi_B', 'psp_A', 'psp_B', 'ppa_A', 'ppa_B',
          'pcam', 'pav', 'pspe', 'pmal', 'par', 'pmu', 'pdu', 'psal', 'ppe', 'gtr', 'rch']
RAC = {'pzi_A', 'pzi_B', 'psp_A', 'psp_B', 'ppa_A', 'ppa_B'}
TE_RE = re.compile(r'transposon|retrotransposon|retrovirus|reverse transcriptase|'
                   r'integrase|transposase|gag-polyprotein|copia|gypsy|'
                   r'pol polyprotein|rnase h', re.I)


def load_id_map(track):
    m = {}
    f = FROZEN / 'id_maps' / f'{track}.id_map.tsv'
    if f.exists():
        with f.open(errors='ignore') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                m[r['stable_id']] = r['original_id']
    return m


def load_annotation(species):
    """species -> {protein_id: (symbol, description, set(go_terms))}"""
    d = {}
    f = ANN / species / 'swissprot_functional_annotation.tsv'
    if not f.exists():
        print(f'  WARNING: missing annotation table: {f}')
        return d
    with f.open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            gn = (r.get('gene_names') or '').strip()
            sym = ''
            for tok in re.split(r'[;\s]+', gn):
                if tok and not tok.startswith('{'):
                    sym = tok
                    break
            desc = (r.get('protein_description') or '').split(';')[0][:150].strip()
            gos = set(re.findall(r'GO:\d{7}', r.get('go_terms') or ''))
            d[r['protein_id']] = (sym, desc, gos)
    return d


def main():
    idmaps = {t: load_id_map(t) for t in RAC}
    lookup = {t: load_annotation(t.split('_')[0] if t in RAC else t) for t in TRACKS}
    for t in TRACKS:
        print(f'  {t}: {len(lookup[t]):,} annotated proteins'
              + (f", id_map {len(idmaps[t]):,}" if t in RAC else ''))

    merged = {}
    with (ORTHO / 'Orthogroups.tsv').open(errors='ignore') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            og = r['Orthogroup']
            symbols, descs, gos, src = Counter(), [], set(), set()
            rac_symbols, rac_descs = Counter(), []
            for tr in TRACKS:
                genes = [x.strip() for x in (r.get(tr) or '').split(',') if x.strip()]
                for g in genes:
                    key = g
                    if tr in RAC:
                        key = idmaps[tr].get(g)
                        if key is None:
                            continue
                    a = lookup[tr].get(key)
                    if not a:
                        continue
                    sym, desc, go = a
                    if sym:
                        symbols[sym] += 1
                        if tr in RAC:
                            rac_symbols[sym] += 1
                    if desc:
                        if desc not in descs:
                            descs.append(desc)
                        if tr in RAC and desc not in rac_descs:
                            rac_descs.append(desc)
                    gos.update(go)
                    src.add(tr)
            # naming policy: racemose annotation first, all-17 tracks as fallback
            symbol = (rac_symbols.most_common(1)[0][0] if rac_symbols
                      else (symbols.most_common(1)[0][0] if symbols else ''))
            desc = (' | '.join(rac_descs[:3]) if rac_descs else ' | '.join(descs[:3]))[:150]
            symbol_all17 = symbols.most_common(1)[0][0] if symbols else ''
            desc_all17 = ' | '.join(descs[:3])[:150]
            te = bool(TE_RE.search(symbol + ' ' + desc)) if (symbol or desc) else False
            merged[og] = dict(og=og, symbol=symbol, description=desc,
                              symbol_all17=symbol_all17, description_all17=desc_all17,
                              go_terms=';'.join(sorted(gos)),
                              source_tracks=';'.join(sorted(src)),
                              n_tracks_annotated=len(src), te=te)

    with (OUT / 'og_annotation_all17.tsv').open('w', newline='', encoding='utf-8') as fh:
        cols = ['og', 'symbol', 'description', 'symbol_all17', 'description_all17',
                'go_terms', 'source_tracks', 'n_tracks_annotated', 'te']
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=cols)
        w.writeheader()
        for og in sorted(merged):
            w.writerow(merged[og])
    print(f'og_annotation_all17.tsv: {len(merged):,} OGs, '
          f'{sum(1 for d in merged.values() if d["symbol"]):,} with symbol, '
          f'{sum(1 for d in merged.values() if d["go_terms"]):,} with GO')

    # join onto the scan table
    scan = list(csv.DictReader((OUT / 'clade_family_scan.tsv').open(), delimiter='\t'))
    stat = Counter()
    with (OUT / 'clade_family_scan_v3.tsv').open('w', newline='', encoding='utf-8') as fh:
        cols = list(scan[0].keys()) + ['symbol17', 'description17', 'go17',
                                       'src_tracks17', 'n_tracks17', 'te17']
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in scan:
            m = merged.get(r['og'], {})
            r2 = dict(r)
            r2.update(symbol17=m.get('symbol', ''), description17=m.get('description', ''),
                      go17=m.get('go_terms', ''), src_tracks17=m.get('source_tracks', ''),
                      n_tracks17=m.get('n_tracks_annotated', 0), te17=m.get('te', False))
            w.writerow(r2)
            old_ann = bool((r['symbol'] or '').strip() or (r['description'] or '').strip())
            new_ann = bool(r2['symbol17'] or r2['description17'])
            if old_ann:
                stat['old_annotated'] += 1
            if new_ann:
                stat['new_annotated'] += 1
            if new_ann and not old_ann:
                stat['gained'] += 1
            if int(r['rac']) > 0 and not new_ann:
                stat['rac_pos_still_unannotated'] += 1
            if int(r['rac']) == 0 and new_ann:
                stat['rac0_now_annotated'] += 1
    print(f'clade_family_scan_v3.tsv: {len(scan):,} families; '
          f'annotated old={stat["old_annotated"]:,} new={stat["new_annotated"]:,} '
          f'(+{stat["gained"]:,} gained); rac=0 now annotated {stat["rac0_now_annotated"]:,}; '
          f'rac>0 still unannotated {stat["rac_pos_still_unannotated"]:,}')


if __name__ == '__main__':
    main()
