#!/usr/bin/env python3
"""Read-only streaming metadata inventory of complete requested group trees.

No directory-depth pruning. Symlinks are recorded and not followed to avoid
cycles or escaping the requested roots. File contents and credentials are not read.
"""
import csv,gzip,os,stat,time,json
from pathlib import Path
from collections import Counter
from common import ROOT,save,log
BASES=['/scratch/svierlin','/scratch/pyoung3','/data/pyoung3','/home/pyoung3']
def main():
    counts=Counter();errors=[];start=time.time();candidates=0
    with gzip.open(ROOT/'provenance/all_group_files.csv.gz','wt',newline='') as out,open(ROOT/'provenance/additional_50AM_candidates.csv','w',newline='') as cf:
        fields=['root','path','type','bytes','uid','gid','mode','mtime_ns','inode'];w=csv.DictWriter(out,fields);w.writeheader();cw=csv.DictWriter(cf,fields);cw.writeheader()
        for base in BASES:
            pending=[base];log('Complete metadata inventory '+base)
            while pending:
                directory=pending.pop()
                try:
                    with os.scandir(directory) as scan:
                        for e in scan:
                            try:
                                s=e.stat(follow_symlinks=False);kind='symlink' if stat.S_ISLNK(s.st_mode) else ('directory' if stat.S_ISDIR(s.st_mode) else 'file')
                                row=dict(root=base,path=e.path,type=kind,bytes=s.st_size,uid=s.st_uid,gid=s.st_gid,mode=oct(stat.S_IMODE(s.st_mode)),mtime_ns=s.st_mtime_ns,inode=s.st_ino);w.writerow(row);counts[kind]+=1
                                if kind=='directory':pending.append(e.path)
                                lower=e.path.lower()
                                if '50am' in lower or ('svierlin' in lower and (e.name.endswith(('.sh','.py','.json','.ipynb')) or e.name in ['AGENTS.md','README.md'])):
                                    cw.writerow(row);candidates+=1
                                if sum(counts.values())%250000==0:log(f'{dict(counts)}; {directory}')
                            except OSError as ex:errors.append(dict(path=e.path,error=str(ex)))
                except OSError as ex:errors.append(dict(path=directory,error=str(ex)))
    save('provenance/full_group_inventory_summary.json',dict(roots=BASES,counts=dict(counts),candidate_rows=candidates,elapsed_seconds=time.time()-start,errors=errors,scope='Complete recursive directory traversal of accessible entries; symlink targets not followed; no source writes or source content reads'))
    log('Complete group metadata inventory finished')
if __name__=='__main__':main()
