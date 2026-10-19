#!/usr/bin/env python3
"""Read bounded original-output archive prefixes and objective tracer inputs."""
import tarfile,zipfile,re,hashlib,json,sys
from pathlib import Path
from common import *

PRISM=Path('/scratch/svierlin/prism-v02')
def main():
    records=[];controls=[];requested=set(pd.read_csv(ROOT/'tables/objective_representative_particles.csv').id.astype(int));requested.update([7753,157594,191136,0,193])
    for p in sorted(Path('/data/pyoung3/groupdata/svierlin/50Am_data').glob('eoc_batch*.tar.gz')):
        log('Bounded prefix '+p.name);rec=dict(source=stamp(p),scope='first particle directory, at most 256 MiB uncompressed archive prefix',files=[]);first=None
        try:
            with tarfile.open(p,'r|gz') as tf:
                for m in tf:
                    parts=Path(m.name).parts;numeric=[s for s in parts if s.isdecimal()]
                    if not numeric:continue
                    ident=int(numeric[-1])
                    if first is None:first=ident;requested.add(ident);rec['particle_id']=ident
                    if ident!=first:break
                    if m.offset_data+m.size>256*1024**2:rec['limit_reached']=True;break
                    name=Path(m.name).name
                    if not m.isfile() or name not in ['control.json','x.txt','y.txt','xz.txt','xtime.txt','prism.log','log.txt']:continue
                    if m.size>64*1024**2:rec.setdefault('oversized_members_skipped',[]).append(name);continue
                    data=tf.extractfile(m).read();dest=ROOT/'provenance/raw_archive_prefixes'/p.stem/str(ident)/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
                    rec['files'].append(dict(member=m.name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),saved=str(dest.relative_to(ROOT))))
                    if name=='control.json':
                        ctrl=json.loads(data);controls.append(dict(archive=str(p),particle=ident,control=ctrl))
            rec['status']='readable_prefix_only'
        except Exception as e:rec.update(status='unavailable_or_corrupt_prefix',error=repr(e))
        records.append(rec)
    save('provenance/archived_control_records.json',controls);save('provenance/archive_prefix_checks.json',records)
    zip_records=[]
    for name in ['50Am_output.zip','50Am_output_tempcorr.zip','50Am_abund.zip']:
        p=PRISM/'input/conditions'/name;log('Inspect zip '+name);rec=dict(source=stamp(p));copied=[]
        try:
            with zipfile.ZipFile(p) as z:
                info=z.infolist();rec['members']=len(info);rec['uncompressed_bytes']=sum(i.file_size for i in info);id_to_name={}
                for i in info:
                    if '__MACOSX' in i.filename:continue
                    m=re.fullmatch(r'(\d+)_(?:trajectory|abund)\.txt',Path(i.filename).name)
                    if m:id_to_name[int(m[1])]=i.filename
                rec['particle_files']=len(id_to_name);rec['requested_missing']=sorted(requested-set(id_to_name))
                for ident in sorted(requested&set(id_to_name)):
                    member=id_to_name[ident];data=z.read(member);dest=ROOT/'provenance/tracer_inputs'/p.stem/Path(member).name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
                    copied.append(dict(id=ident,member=member,sha256=hashlib.sha256(data).hexdigest(),saved=str(dest.relative_to(ROOT))))
            rec['status']='selected_members_CRC_verified';rec['copied']=copied
        except Exception as e:rec.update(status='error',error=repr(e))
        zip_records.append(rec)
    save('provenance/input_zip_inventory.json',zip_records)
    # Template extent is an available configuration, not proof of all historical runs.
    extent=[]
    for line in (PRISM/'input/extent/default.dat').read_text().splitlines():
        z,lo,hi=map(int,line.split())
        for n in range(lo,hi+1):extent.append(dict(Z=z,N=n,A=z+n,species=('n' if z==0 and n==1 else str(z+n)+(SYMBOLS[z] if z<len(SYMBOLS) else 'Z'+str(z)))))
    csvsave('tables/prism_default_extent_species.csv',extent)
    save('provenance/archive_probe_summary.json',dict(requested_tracers=sorted(requested),archived_controls=len(controls),template_extent_species=len(extent),warning='Archive prefixes and selected tracer inputs do not establish a complete archive integrity or historical-run census.'))
    log('Archive/input probes complete')
if __name__=='__main__':main()
