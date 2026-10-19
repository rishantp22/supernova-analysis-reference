#!/usr/bin/env python3
"""Stream every final PRISM EOC CSV from its compact archive and compare by ID,Z,A.

These are derivative x/y extracts, not substitutes for historical run controls.
Preserve explicit NaN and absent rows separately in the ragged NPZ database.
"""
import tarfile,hashlib,re,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from common import *

SOURCE=Path('/scratch/svierlin/prism-v02/EOC_csv.tar.gz')
def main():
    out=ROOT/'tables/eoc';out.mkdir(exist_ok=True)
    meta,a,rawids=table('new');cols=meta['columns'];nu,nrow=unique_rows(rawids);rowmap={int(x):int(y) for x,y in zip(nu,nrow)}
    isotopecols=[i for i,c in enumerate(cols) if species(c)[0]=='isotope'];codes=np.array([species(cols[i])[1]*512+species(cols[i])[2] for i in isotopecols]);code_to_column={int(c):i for i,c in enumerate(codes)}
    summary=Counter();per_species={int(c):Counter() for c in codes};all_species=defaultdict(Counter);all_ids=[];sums=[];yes=[];mismatches=[];missingcounts=[];statusrows=[];examples=[];members=[];blocks=[];block_ids=[];mass_yield=defaultdict(float);invalid=[]
    before=stamp(SOURCE);selected=set(pd.read_csv(ROOT/'tables/objective_representative_particles.csv').id.astype(int));selected.update([7753,157594,191136])
    block_number=0
    def flush():
        nonlocal block_number
        if not blocks:return
        offsets=np.r_[0,np.cumsum([len(v) for v in blocks])];v=np.concatenate(blocks)
        np.savez_compressed(out/f'block_{block_number:05d}.npz',id=np.array(block_ids,dtype='i8'),offset=offsets,Z=v[:,0].astype('i2'),A=v[:,1].astype('i2'),X=v[:,2],Y=v[:,3])
        block_number+=1
        blocks.clear();block_ids.clear()
    with tarfile.open(SOURCE,'r|gz') as tf:
        for member in tf:
            if not member.isfile() or not re.fullmatch(r'(?:.*/)?\d+\.csv',member.name):continue
            ident=int(Path(member.name).stem);payload=tf.extractfile(member).read()
            if not payload.startswith(b'Z,A,x,y\n'):raise ValueError('Unexpected EOC header '+member.name)
            data=np.fromstring(payload.split(b'\n',1)[1].replace(b',',b' '),sep=' ').reshape(-1,4)
            z=data[:,0].astype('i4');A=data[:,1].astype('i4');code=z*512+A;x=data[:,2];y=data[:,3]
            if len(np.unique(code))!=len(code):raise ValueError('Duplicate species '+member.name)
            if (x[np.isfinite(x)]<0).any():invalid.append(dict(id=ident,issue='negative_x'))
            valid=np.isfinite(x);xy=valid&np.isfinite(y)
            sumx=float(x[valid].sum());ye=float(np.sum(x[valid]*z[valid]/A[valid]));all_ids.append(ident);sums.append(sumx);yes.append(ye);missingcounts.append(int((~valid).sum()))
            summary['files']+=1;summary['species_rows']+=len(x);summary['finite_x_rows']+=int(valid.sum());summary['xy_pairs']+=int(xy.sum());summary['xy_mismatches_rtol1e12']+=int((abs(x[xy]-A[xy]*y[xy])>1e-12*abs(x[xy])).sum())
            mass=float(a[rowmap[ident],cols.index('mass')]/MSUN) if ident in rowmap else None
            for c,v in zip(code,x):
                q=all_species[int(c)];q['present_rows']+=1;q['finite_x']+=int(np.isfinite(v));q['positive_x']+=int(v>0)
                if mass is not None and np.isfinite(v):mass_yield[int(c)]+=mass*v
            count=0;row_status='absent_from_valid_new_ids'
            if ident in rowmap:
                plot=np.asarray(a[rowmap[ident],isotopecols]);lookup={int(c):float(v) for c,v in zip(code,x)}
                ref=np.array([lookup.get(int(c),np.nan) for c in codes]);good=np.isfinite(ref)&np.isfinite(plot);bad=good&(abs(plot-ref)>5.1e-6*abs(ref));count=int(bad.sum());row_status='compared'
                summary['compared_particles']+=1;summary['finite_pairs']+=int(good.sum());summary['mismatch_cells']+=count;summary['mismatch_particles']+=int(count>0)
                censored=(~np.isfinite(ref))&(plot==0);summary['zero_plot_missing_final_x']+=int(censored.sum());summary['plot_nonzero_missing_final_x']+=int(((~np.isfinite(ref))&(plot>0)).sum())
                for j,c in enumerate(codes):
                    q=per_species[int(c)];q['finite_pairs']+=int(good[j]);q['mismatch']+=int(bad[j]);q['zero_plot_missing_final_x']+=int(censored[j]);q['plot_nonzero_missing_final_x']+=int(not np.isfinite(ref[j]) and plot[j]>0)
                if count and len(examples)<40:
                    j=int(np.nanargmax(np.where(bad,abs(plot-ref),np.nan)));examples.append(dict(id=ident,species=cols[isotopecols[j]],final_x=ref[j],plotting=plot[j],final_sum=sumx))
            mismatches.append(count)
            if ident in selected:
                target=ROOT/'provenance/representative_eoc';target.mkdir(exist_ok=True);(target/f'{ident}.csv').write_bytes(payload)
            members.append(dict(id=ident,member=member.name,bytes=member.size,sha256=hashlib.sha256(payload).hexdigest()))
            blocks.append(data);block_ids.append(ident)
            if len(blocks)>=5000:flush()
            if summary['files']%25000==0:log(f'EOC: {summary["files"]:,} particles; {summary["mismatch_cells"]:,} mismatched cells')
    flush();assert stamp(SOURCE)==before
    u,c=np.unique(all_ids,return_counts=True);summary['duplicate_ids']=int((c>1).sum());summary['unique_ids']=len(u);summary['new_ids_without_eoc']=len(np.setdiff1d(nu,u));summary['eoc_ids_without_new']=len(np.setdiff1d(u,nu))
    np.savez_compressed(ROOT/'tables/eoc_particle_validation.npz',id=np.array(all_ids),sumX=np.array(sums),Ye=np.array(yes),missingX=np.array(missingcounts),plot_mismatches=np.array(mismatches))
    csvsave('tables/eoc_species_inventory.csv',[dict(code=c,Z=c//512,A=c%512,species=('n' if c==1 else str(c%512)+(SYMBOLS[c//512] if c//512<len(SYMBOLS) else 'Z'+str(c//512))),yield_msun_matched=mass_yield.get(c),**v) for c,v in sorted(all_species.items())])
    csvsave('tables/eoc_vs_new_by_species.csv',[dict(Z=c//512,A=c%512,species=cols[isotopecols[code_to_column[c]]],**v) for c,v in per_species.items()]);csvsave('provenance/eoc_member_hashes.csv',members);csvsave('tables/eoc_duplicate_ids.csv',{'id':u[c>1],'count':c[c>1]});csvsave('tables/eoc_invalid.csv',invalid)
    save('provenance/eoc_summary.json',dict(source=before,archive_sha256=sha(SOURCE),counts=dict(summary),sumX=percentiles(sums),Ye=percentiles(yes),sum_over_unity_1e5=int((np.array(sums)>1.00001).sum()),sum_under_unity_1e5=int((np.array(sums)<.99999).sum()),examples=examples,comparison_rtol=5.1e-6,atol=0))
    log('Full EOC archive audit complete')
if __name__=='__main__':main()
