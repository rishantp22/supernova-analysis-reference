#!/usr/bin/env python3
"""Copy an explicit lightweight research reference into a fresh local repository."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

SCRATCH = Path('/scratch/rpraka15')
DEST = SCRATCH/'supernova-analysis-reference'
HERE = Path(__file__).resolve().parent
SMALL_TABLES = {'hydro_diagnostics.csv','particle_id_join_demo.csv','prism_diagnostics.csv',
                'requested_field_map.csv','snapshot_catalog.csv','unit_crosscheck.csv'}


def choose(rel):
    p = rel.parts
    if p[0] == 'supernova-grain-analysis':
        if p[1] == 'code':
            return rel.suffix in {'.py','.md','.txt'}
        if p[1] == 'reports':
            return rel.suffix == '.md'
        if p[1] == 'debug':
            return rel.name != 'report_embedded_data.json' and rel.suffix in {'.json','.txt','.js'}
        if p[1] == 'inventory':
            return rel.name in {'inventory.json','sample_ids.json','archive_log_excerpts.txt'}
        if p[1] == 'tables':
            return rel.name in SMALL_TABLES
    if p[0] == 'final_chunks_export_20261018':
        if 'parts' in p or '__pycache__' in p:
            return False
        if rel.suffix == '.py':
            return True
        if rel.suffix in {'.txt','.json'}:
            return rel.name != 'chunk_validation.json'
        return rel.name in {'hydro_ids_without_chunks.csv','mass_fraction_summary.csv','column_map.csv'}
    return False


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison-results',type=Path)
    args=parser.parse_args()
    DEST.mkdir(parents=True,exist_ok=True)
    (DEST/'provenance').mkdir(exist_ok=True)
    copied=[]; excluded=[]
    def copy(source, relative, notebook=False):
        before=source.stat()
        data=source.read_bytes()
        if (source.stat().st_size,source.stat().st_mtime_ns)!=(before.st_size,before.st_mtime_ns):
            raise RuntimeError(f'Source changed: {source}')
        prepared=data
        if notebook:
            nb=json.loads(data)
            nb.get('metadata',{}).pop('widgets',None)
            for cell in nb['cells']:
                if cell.get('cell_type')=='code':
                    cell['outputs']=[];cell['execution_count']=None
            prepared=(json.dumps(nb,indent=1,ensure_ascii=False)+'\n').encode()
        if len(prepared)>1024*1024:
            raise RuntimeError(f'Unexpected large upload candidate: {source}')
        target=DEST/relative;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(prepared)
        copied.append(dict(source_path=str(source),source_bytes=len(data),source_sha256=digest(data),
                           upload_path=str(relative),upload_bytes=len(prepared),upload_sha256=digest(prepared),
                           transformation='notebook outputs/execution counts cleared in copy' if notebook else 'exact byte copy'))
    for root_name in ['supernova-grain-analysis','final_chunks_export_20261018']:
        root=SCRATCH/root_name
        result=subprocess.run(['git','-C',str(root),'rev-parse','--show-toplevel'],capture_output=True,text=True)
        if result.returncode==0:
            raise RuntimeError(f'Source unexpectedly belongs to Git; review tracking first: {root}')
        for source in sorted(root.rglob('*')):
            if not source.is_file():continue
            rel=source.relative_to(SCRATCH)
            if choose(rel) and source.stat().st_size<=1024*1024:
                copy(source,rel)
            else:
                excluded.append(dict(source_path=str(source),bytes=source.stat().st_size,
                                     reason='dataset, cache, bulky generated output or redundant intermediate',already_git_tracked=False))
    for source,target,notebook in [
        (SCRATCH/'analyze.ipynb',Path('notebooks/analyze.ipynb'),True),
        (SCRATCH/'output_data/analyze.ipynb',Path('notebooks/output_data_analyze.ipynb'),True),
        (SCRATCH/'move_final_chunks.sh',Path('slurm/historical_move_final_chunks.sh'),False),
        (Path('/home/rpraka15/supernova-grain-analysis/tests/test_readers.py'),Path('tests/test_readers.py'),False),
        (SCRATCH/'List of isotopes_elements_variables.txt',Path('schema/isotopes_elements_variables.txt'),False)]:
        copy(source,target,notebook)
    copy(Path(__file__).resolve(),Path('tools/prepare_reference.py'))
    comparison=SCRATCH/'file_comparison_20261018'
    for name in ['compare_outputs.py','test_compare_outputs.py','compare_outputs.sbatch','README.txt','environment.json','test-results.txt']:
        copy(comparison/name,Path('comparison')/name)
    if args.comparison_results:
        for source in sorted(args.comparison_results.iterdir()):
            if source.suffix in {'.json','.csv'} and source.stat().st_size<=1024*1024:
                copy(source,Path('comparison/validation')/source.name)
    for name in ['final_chunks.out','50Am_plotting_all_isotopes.out','50Am_plotting_complete.out',
                 'run3g_50Am6_sph.2000.out','run3g_50Am_jet3b_sph_2000_selected_fields.csv','50Am_UCF_OneDrive.tar.gz']:
        source=SCRATCH/name
        if source.exists():
            excluded.append(dict(source_path=str(source),bytes=source.stat().st_size,
                                 reason='large dataset or packaged research data',already_git_tracked=False))
    for name,rows in [('source_manifest.csv',copied),('excluded_files.csv',excluded)]:
        with (DEST/'provenance'/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=dict(copied_files=len(copied),copied_bytes=sum(r['upload_bytes'] for r in copied),
                 excluded_files=len(excluded),excluded_bytes=sum(r['bytes'] for r in excluded),
                 source_files_unchanged=True,destination=str(DEST),
                 proposed_github_repository='rishantp22/supernova-analysis-reference',
                 github_repository_created=False,uploaded=False)
    (DEST/'provenance/preparation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
