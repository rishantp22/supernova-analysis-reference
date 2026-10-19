#!/usr/bin/env python3
"""Bounded directory discovery, explicit omissions, and small source snapshots."""
import os,re,json,subprocess,hashlib
from pathlib import Path
from common import *
BASES=['/scratch/rpraka15','/home/rpraka15/supernova-grain-analysis','/scratch/svierlin/prism-v02','/scratch/pyoung3/runsnsph','/data/pyoung3/supernova_data/50Am','/data/pyoung3/gsvance/data_mining','/data/pyoung3/groupdata/svierlin/50Am_data','/data/pyoung3/burn','/data/pyoung3/burnpatrick','/data/pyoung3/burnpkg','/data/pyoung3/nucleosetup.dir','/data/pyoung3/setraj','/data/pyoung3/se_analysis','/data/pyoung3/snsphtools']
def main():
    files=[];omissions=[];sources=[];repos=[];seen=set()
    extensions={'.py','.sh','.bash','.f','.f90','.c','.h','.md','.json','.ctl','.sbatch','.slurm','.ipynb','.f95'}
    for base in BASES:
        log('Discover '+base)
        def error(e):omissions.append(dict(path=e.filename,reason=str(e)))
        for root,dirs,names in os.walk(base,followlinks=False,onerror=error):
            p=Path(root);depth=len(p.relative_to(base).parts)
            if str(p)==str(ROOT):dirs[:]=[];continue
            if '.git' in dirs:
                def git(*a):
                    r=subprocess.run(['git','-C',str(p),*a],capture_output=True,text=True)
                    return r.stdout.strip() if r.returncode==0 else r.stderr.strip()
                repos.append(dict(path=str(p),branch=git('branch','--show-current'),head=git('rev-parse','HEAD'),recent=git('log','-5','--oneline'),remotes=re.sub(r'https://[^/@\s]+:[^/@\s]+@','https://REDACTED@',git('remote','-v'))))
            keep=[]
            for d in dirs:
                reason=None
                if d.startswith('.') or d in ['__pycache__','software','node_modules']:reason='hidden/cache/software directory'
                elif d.isdecimal():reason='per-particle directory; selected by ID in separate audit'
                elif depth>=5:reason='discovery depth bound; targeted follow-up required'
                elif d in ['jet3b_outputs','cco2_outputs','vconv_outputs','jet3b_traj','cco2_traj','vconv_traj','50Am_abund','final_chunks','final_chunks_B','parts','slurm_log','logs','sorted_queries','queries','slurm_logs']:reason='large collection; separately audited where relevant'
                if reason:omissions.append(dict(path=str(p/d),reason=reason))
                else:keep.append(d)
            dirs[:]=keep
            for name in names:
                q=p/name
                if name.startswith('.') or str(q) in seen:continue
                seen.add(str(q))
                try:
                    s=q.stat(); row=dict(path=str(q),bytes=s.st_size,readable=os.access(q,os.R_OK),symlink=q.is_symlink())
                    files.append(row)
                    if s.st_size<=2*1024*1024 and (q.suffix in extensions or name in ['Makefile','README','INSTRUCTIONS','params','input']) and not re.search(r'credential|token|secret|password',name,re.I):
                        try:data=q.read_bytes()
                        except OSError as e:omissions.append(dict(path=str(q),reason=str(e)));continue
                        if q.suffix=='.ipynb':
                            nb=json.loads(data);data=('\n\n'.join(''.join(c['source']) for c in nb.get('cells',[]) if c.get('cell_type')=='code')).encode()
                        dest=ROOT/'provenance/source_snapshots'/q.relative_to('/')
                        if q.suffix=='.ipynb':dest=dest.with_suffix('.code.py')
                        dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
                        sources.append(dict(source=str(q),snapshot=str(dest.relative_to(ROOT)),snapshot_sha256=hashlib.sha256(data).hexdigest(),source_bytes=s.st_size))
                except OSError as e:omissions.append(dict(path=str(q),reason=str(e)))
    csvsave('provenance/discovered_files.csv',files);csvsave('provenance/search_omissions.csv',omissions);csvsave('provenance/source_snapshots.csv',sources);save('provenance/git_histories.json',repos)
    matches=[r for r in files if re.sub('[^a-z0-9]','',Path(r['path']).name.lower())=='50amplottingout']
    save('provenance/legacy_filename_search.json',dict(matches=matches,bases=BASES,scope='bounded five-level search, pruned collections documented separately; not a claim about inaccessible directories'))
    log(f'Discovery complete: {len(files)} files, {len(sources)} source snapshots, {len(omissions)} documented omissions')
if __name__=='__main__':main()
