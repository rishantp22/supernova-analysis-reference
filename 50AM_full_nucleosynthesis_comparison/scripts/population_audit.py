#!/usr/bin/env python3
"""Full-population audit; column-oriented memmaps and one-species-at-a-time statistics.

Input files are scanned, hashed and parsed without repairing invalid records.
NPZ files retain every matched ID, all numerical differences and undefined values.
No positive floor is used. Zero denominators and nonpositive log arguments are NA.
"""
import argparse, hashlib, json, platform, re, resource, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
from common import *

INPUTS={'old':Path('/data/pyoung3/supernova_data/50Am/analysis/50Am_plotting.out'),
        'new':Path('/scratch/rpraka15/50Am_plotting_all_isotopes.out')}

def ingest(label,path):
    before=stamp(path);cache=ROOT/f'provenance/{label}_scan.json'
    if cache.exists():
        m=json.loads(cache.read_text())
        if m['source']==before:log(f'Reusing validated {label} parsed cache');return
        raise RuntimeError('Source changed since cached analysis: '+str(path))
    h=hashlib.sha256();rows=0;width_errors=[];blank=[]
    log(f'Hashing and checking record widths: {label}')
    with path.open('rb') as f:
        header=f.readline();h.update(header)
        rawcols=[s.strip() for s in header.decode().strip().split(',')];cols=[canon(s) for s in rawcols]
        assert len(cols)==len(set(cols)) and cols[0]=='id'
        for lineno,line in enumerate(f,2):
            h.update(line)
            if not line.strip():blank.append(lineno);continue
            if len(line.rstrip(b'\r\n').split(b','))!=len(cols):width_errors.append(lineno)
            rows+=1
    assert not width_errors, f'Wrong-width rows require explicit parser review: {width_errors[:10]}'
    a=np.lib.format.open_memmap(ROOT/f'tables/{label}_values.npy',mode='w+',dtype='f8',shape=(rows,len(cols)),fortran_order=True)
    ids=np.lib.format.open_memmap(ROOT/f'tables/{label}_ids.npy',mode='w+',dtype='i8',shape=(rows,))
    invalid=[];counts={c:dict(NA=0,invalid=0,negative=0,zero=0,positive=0) for c in cols};offset=0
    for frame in pd.read_csv(path,skipinitialspace=True,dtype=str,keep_default_na=False,chunksize=8192):
        n=len(frame)
        for j,c in enumerate(cols):
            tokens=frame.iloc[:,j].str.strip()
            missing=tokens.eq('NA') | tokens.eq('')
            val=pd.to_numeric(tokens.mask(missing),errors='coerce').to_numpy(dtype='f8',copy=True)
            bad=(~missing.to_numpy())&(~np.isfinite(val))
            if j==0:
                valid=tokens.str.fullmatch(r'\d+').to_numpy() & np.isfinite(val) & (val < 2**53)
                bad=~valid
                ids[offset:offset+n]=np.where(valid,val,-1).astype('i8')
            for k in np.flatnonzero(bad):invalid.append(dict(row_index=offset+int(k),line=offset+int(k)+2,column=c,token=tokens.iloc[k]))
            val[bad]=np.nan;a[offset:offset+n,j]=val
            counts[c]['NA']+=int(missing.sum());counts[c]['invalid']+=int(bad.sum())
            for key,mask in [('negative',val<0),('zero',val==0),('positive',val>0)]:counts[c][key]+=int(mask.sum())
        offset+=n
        if offset//100000!=(offset-n)//100000:log(f'{label}: parsed {offset:,}/{rows:,}')
    assert offset==rows and before==stamp(path)
    a.flush();ids.flush()
    valid_ids=np.asarray(ids)[ids>=0];u,c=np.unique(valid_ids,return_counts=True)
    csvsave(f'tables/{label}_duplicate_ids.csv',[dict(id=int(x),count=int(y)) for x,y in zip(u[c>1],c[c>1])])
    csvsave(f'tables/{label}_invalid_cells.csv',invalid)
    save(f'provenance/{label}_scan.json',dict(source=before,sha256=h.hexdigest(),rows=rows,raw_columns=rawcols,columns=cols,unique_valid_ids=len(u),duplicate_ids=int((c>1).sum()),invalid_ids=int((ids<0).sum()),id_inversions=int((np.diff(valid_ids)<0).sum()),blank_lines=blank,wrong_width_lines=width_errors,column_counts=counts,invalid_cells=invalid))

def stats(old,new,ids,mass_old,mass_new,rtol):
    finite=np.isfinite(old)&np.isfinite(new);o=old[finite];n=new[finite];ii=ids[finite]
    d=n-o;exact=n==o;approx=np.abs(d)<=rtol*np.abs(o)
    ratio=np.full(len(o),np.nan);rel=ratio.copy();lg=ratio.copy()
    with np.errstate(divide='ignore',invalid='ignore',over='ignore',under='ignore'):
        np.divide(n,o,out=ratio,where=o!=0);np.divide(d,o,out=rel,where=o!=0)
        pos=(n>0)&(o>0);lg[pos]=np.log10(n[pos])-np.log10(o[pos])
    wm=np.isfinite(mass_old[finite])&np.isfinite(mass_new[finite])&(mass_old[finite]>0)&(mass_new[finite]>0)
    mo=mass_old[finite][wm]/MSUN;mn=mass_new[finite][wm]/MSUN
    yo=np.sum(mo*o[wm]);yn=np.sum(mn*n[wm]);dis=~approx[wm]
    row=dict(matched_ids=len(ids),finite_pairs=len(o),both_missing=int((~np.isfinite(old)&~np.isfinite(new)).sum()),missing_old_only=int((~np.isfinite(old)&np.isfinite(new)).sum()),missing_new_only=int((np.isfinite(old)&~np.isfinite(new)).sum()),exact=int(exact.sum()),approximate_not_exact=int((approx&~exact).sum()),significant=int((~approx).sum()),increasing=int((d>0).sum()),decreasing=int((d<0).sum()),newly_nonzero=int(((o==0)&(n>0)).sum()),newly_zero=int(((o>0)&(n==0)).sum()),both_zero=int(((o==0)&(n==0)).sum()),negative_old=int((o<0).sum()),negative_new=int((n<0).sum()),ratio_defined=int((o!=0).sum()),ratio_overflow=int(np.isinf(ratio).sum()),log_defined=int(pos.sum()),rtol=rtol,atol=0.,old_yield_shared_msun=yo,new_yield_shared_msun=yn,yield_ratio_shared=yn/yo if yo!=0 else None,absolute_particle_yield_change_msun=float(np.sum(np.abs(mn*n[wm]-mo*o[wm]))),signed_yield_change_msun=yn-yo,discrepant_old_yield_msun=np.sum(mo[dis]*o[wm][dis]),discrepant_new_yield_msun=np.sum(mn[dis]*n[wm][dis]))
    if not wm.any():
        for k in list(row):
            if 'yield' in k:row[k]=None
    for key,v in [('signed_delta',d),('absolute_delta',abs(d)),('relative_delta',rel),('ratio',ratio),('log10_ratio',lg)]:
        for q,x in percentiles(v).items():row[key+'_'+q]=x
    if len(o):
        k=int(np.argmax(abs(d)));row.update(max_abs_id=int(ii[k]),max_abs_old=float(o[k]),max_abs_new=float(n[k]))
        if np.any(pos):k=int(np.nanargmax(abs(lg)));row.update(max_log_id=int(ii[k]),max_abs_log10_ratio=float(abs(lg[k])))
    full={k:np.full(len(ids),np.nan) for k in ['delta','relative_delta','ratio','log10_ratio']}
    for k,v in zip(full,[d,rel,ratio,lg]):full[k][finite]=v
    full.update(old=old,new=new,valid_pair=finite,significant=np.where(finite,abs(new-old)>rtol*abs(old),False))
    return row,full

def run_stats():
    om,oa,oi=table('old');nm,na,ni=table('new');ids,ro,rn=matched()
    ou,orow=unique_rows(oi);nu,nrow=unique_rows(ni)
    csvsave('tables/ids_only_old.csv',{'id':np.setdiff1d(ou,nu)})
    csvsave('tables/ids_only_new.csv',{'id':np.setdiff1d(nu,ou)})
    np.savez_compressed(ROOT/'tables/alignment.npz',id=ids,old_row=ro,new_row=rn)
    mo=oa[ro,om['columns'].index('mass')];mn=na[rn,nm['columns'].index('mass')]
    inventory=[];summaries=[];hydro=[];all_yields=[];bins=[];examples=[]
    mapp=lambda label,c: (oa[:,om['columns'].index(c)] if label=='old' else na[:,nm['columns'].index(c)])
    for c in sorted(set(om['columns'])|set(nm['columns']),key=lambda c:(species(c)[0],species(c)[1] or -1,species(c)[2] or -1,c)):
        kind,z,a=species(c);inold=c in om['columns'];innew=c in nm['columns']
        inventory.append(dict(name=c,kind=kind,Z=z,A=a,in_old=inold,in_new=innew,old_name=om['raw_columns'][om['columns'].index(c)] if inold else '',new_name=nm['raw_columns'][nm['columns'].index(c)] if innew else '',old_finite=om['rows']-om['column_counts'][c]['NA']-om['column_counts'][c]['invalid'] if inold else 0,new_finite=nm['rows']-nm['column_counts'][c]['NA']-nm['column_counts'][c]['invalid'] if innew else 0))
        if kind!='physical':
            for label,meta,data,rows in [('old',om,oa,orow),('new',nm,na,nrow)]:
                if c not in meta['columns']:continue
                x=data[rows,meta['columns'].index(c)];mass=data[rows,meta['columns'].index('mass')]
                valid=np.isfinite(x)&np.isfinite(mass)&(mass>0)
                all_yields.append(dict(dataset=label,species=c,kind=kind,unique_particle_rows=len(rows),available_values=int(valid.sum()),zero=int((x==0).sum()),positive=int((x>0).sum()),negative=int((x<0).sum()),yield_msun=float(np.sum(x[valid]*mass[valid]/MSUN)) if valid.any() else None,**{k:v for k,v in percentiles(x).items()}))
        if not (inold and innew) or c=='id':continue
        old=np.asarray(mapp('old',c)[ro]);new=np.asarray(mapp('new',c)[rn]);rtol=1e-9 if kind=='physical' else 1e-5
        row,full=stats(old,new,ids,mo,mn,rtol);row.update(species=c,kind=kind,Z=z,A=a)
        (hydro if kind=='physical' else summaries).append(row)
        if kind=='isotope':
            out=ROOT/'tables/species';out.mkdir(exist_ok=True)
            np.savez_compressed(out/(c+'.npz'),**full)
            delta=np.abs(new-old);ix=np.argsort(np.where(np.isfinite(delta),delta,-1))[-3:][::-1]
            for k in ix:
                if np.isfinite(delta[k]):examples.append(dict(species=c,id=int(ids[k]),old=old[k],new=new[k],abs_delta=delta[k],criterion='top3_absolute_mass_fraction_difference'))
            for field,edges in [('peak temp',np.array([0,1e9,2e9,3e9,4e9,5e9,7e9,10e9,np.inf])),('peak density',np.r_[0,np.logspace(0,10,11),np.inf]),('Y_{e}',np.r_[0,.4,.45,.48,.49,.495,.5,.505,.51,.55,1,np.inf]),('mass',np.r_[0,np.logspace(20,32,13),np.inf])]:
                env=oa[ro,om['columns'].index(field)]
                for lo,hi in zip(edges[:-1],edges[1:]):
                    ok=(env>=lo)&(env<hi)&np.isfinite(old)&np.isfinite(new)
                    bins.append(dict(species=c,field=field,lower=lo,upper=hi,particles=int(ok.sum()),significant=int((ok&full['significant']).sum()),old_yield_msun=np.sum(mo[ok]*old[ok]/MSUN),new_yield_msun=np.sum(mn[ok]*new[ok]/MSUN),mean_abs_delta=float(np.mean(delta[ok])) if ok.any() else None))
        log(f'Compared {c}: {row["finite_pairs"]:,} pairs; {row["significant"]:,} significant')
    csvsave('tables/species_inventory.csv',inventory);csvsave('tables/species_statistics.csv',summaries);csvsave('tables/physical_statistics.csv',hydro);csvsave('tables/yields_all_available_particles.csv',all_yields);csvsave('tables/thermodynamic_binned_statistics.csv',bins);csvsave('tables/objective_representative_particles.csv',examples)
    save('provenance/population_summary.json',dict(shared_unique_valid_ids=len(ids),old_only_ids=len(np.setdiff1d(ou,nu)),new_only_ids=len(np.setdiff1d(nu,ou)),old_unique_valid_ids=len(ou),new_unique_valid_ids=len(nu),shared_particle_mass_old_msun=np.sum(mo)/MSUN,shared_particle_mass_new_msun=np.sum(mn)/MSUN,shared_isotopes=sum(r['kind']=='isotope' and r['in_old'] and r['in_new'] for r in inventory),new_isotopes=sum(r['kind']=='isotope' and r['in_new'] for r in inventory),old_isotopes=sum(r['kind']=='isotope' and r['in_old'] for r in inventory),tolerance={'physical':{'rtol':1e-9,'atol':0},'composition':{'rtol':1e-5,'atol':0}},maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

def consistency():
    ids,ro,rn=matched();closure=[];element_rows=[];ratio_rows=[];decay_rows=[]
    ratio_dir=ROOT/'tables/ratios';ratio_dir.mkdir(exist_ok=True)
    datasets={}
    for label in ['old','new']:
        meta,a,rawids=table(label);u,rows=unique_rows(rawids);cols=meta['columns']
        total=np.zeros(len(u));ye=np.zeros(len(u));missing=np.zeros(len(u),dtype='i4');elements={}
        for c in cols:
            kind,z,A=species(c)
            if kind!='isotope':continue
            x=np.asarray(a[rows,cols.index(c)]);finite=np.isfinite(x);total+=np.where(finite,x,0);ye+=np.where(finite,x*z/A,0);missing+=~finite
            if z not in elements:elements[z]=np.zeros(len(u))
            elements[z]+=np.where(finite,x,0)
        # Sum available isotope values only. Missing columns remain recorded and
        # these are explicitly partial sums, not normalization/Ye assertions.
        reported=a[rows,cols.index('Y_{e}')]
        np.savez_compressed(ROOT/f'tables/{label}_closure.npz',id=u,selected_isotope_sum=total,unaccounted_mass=1-total,partial_composition_ye=ye,reported_progenitor_ye=reported,missing_isotope_columns=missing)
        closure.append(dict(dataset=label,particles=len(u),isotope_columns=sum(species(c)[0]=='isotope' for c in cols),rows_missing=int((missing>0).sum()),over_unity_1e5=int((total>1.00001).sum()),under_unity_1e5=int((total<.99999).sum()),within_1e5=int((abs(total-1)<=1e-5).sum()),max_sum_id=int(u[np.nanargmax(total)]),**{'sum_'+k:v for k,v in percentiles(total).items()},**{'partial_ye_'+k:v for k,v in percentiles(ye).items()}))
        for z,subtotal in elements.items():
            el=SYMBOLS[z]
            if z==0:continue
            e=a[rows,cols.index(el)] if el in cols else np.full(len(u),np.nan)
            finite=np.isfinite(e)
            element_rows.append(dict(dataset=label,element=el,finite_element_totals=int(finite.sum()),selected_isotopes=sum(species(c)[0]=='isotope' and species(c)[1]==z for c in cols),subtotal_exceeds_reported=int((finite&(subtotal>e+1e-5*abs(e))).sum()),matches_selected_subtotal=int((finite&(abs(subtotal-e)<=1e-5*abs(e))).sum()),**{'subtotal_minus_total_'+k:v for k,v in percentiles(subtotal[finite]-e[finite]).items()}))
        datasets[label]=(meta,a,ro if label=='old' else rn)
    csvsave('tables/closure_statistics.csv',closure);csvsave('tables/elemental_total_checks.csv',element_rows)
    # Every isotope versus the highest population-integrated mass-fraction reference
    # in its element; one fixed reference for both datasets. Explicit interest pairs added.
    om,oa,_=datasets['old'];nm,na,_=datasets['new'];union=set(om['columns'])|set(nm['columns'])
    pairs=set();yieldtable=pd.read_csv(ROOT/'tables/yields_all_available_particles.csv')
    for z in range(1,len(SYMBOLS)):
        isotopes=[c for c in union if species(c)[0]=='isotope' and species(c)[1]==z]
        if len(isotopes)<2:continue
        candidates=yieldtable[(yieldtable.dataset=='new')&yieldtable.species.isin(isotopes)]
        if candidates.empty:candidates=yieldtable[yieldtable.species.isin(isotopes)]
        ref=str(candidates.loc[candidates.yield_msun.idxmax(),'species'])
        pairs.update((c,ref) for c in isotopes if c!=ref)
    pairs.update([('26Al','27Al'),('25Mg','24Mg'),('26Mg','24Mg'),('48Ca','40Ca'),('46Ca','40Ca')]+[(c,'48Ti') for c in ['46Ti','47Ti','49Ti','50Ti']]+[(c,'52Cr') for c in ['50Cr','53Cr','54Cr']]+[(c,'56Fe') for c in ['54Fe','57Fe','58Fe']])
    for num,den in sorted(pairs):
        vals={};status={}
        for label,(meta,a,rows) in datasets.items():
            if num not in meta['columns'] or den not in meta['columns']:vals[label]=np.full(len(ids),np.nan);status[label]='missing_column';continue
            x=a[rows,meta['columns'].index(num)];y=a[rows,meta['columns'].index(den)];ok=np.isfinite(x)&np.isfinite(y)&(y>0)&(x>=0)
            v=np.full(len(ids),np.nan)
            with np.errstate(over='ignore',under='ignore',invalid='ignore'):v[ok]=(x[ok]/y[ok])*(species(den)[2]/species(num)[2])
            vals[label]=v;status[label]=dict(valid=int(ok.sum()),zero_denominator=int((y==0).sum()),missing=int((~np.isfinite(x)|~np.isfinite(y)).sum()),overflow=int(np.isinf(v).sum()))
        row,_=stats(vals['old'],vals['new'],ids,np.ones(len(ids)),np.ones(len(ids)),2e-5)
        row={k:v for k,v in row.items() if 'yield' not in k};row.update(numerator=num,denominator=den,old_status=json.dumps(status.get('old')),new_status=json.dumps(status.get('new')))
        ratio_rows.append(row);np.savez_compressed(ratio_dir/f'{num}_over_{den}.npz',old=vals['old'],new=vals['new'])
    csvsave('tables/number_ratio_statistics.csv',ratio_rows)
    # A simple closed-chain mass-conservation diagnostic, never a decay correction.
    for group in [('26Al','26Mg'),('41Ca','41K'),('44Ti','44Sc','44Ca'),('53Mn','53Cr'),('56Ni','56Co','56Fe'),('57Ni','57Co','57Fe')]:
        record={'chain':' -> '.join(group)}
        sums={}
        for label,(meta,a,rows) in datasets.items():
            available=[c for c in group if c in meta['columns']];record[label+'_available']=';'.join(available);record[label+'_missing']=';'.join(c for c in group if c not in available)
            sums[label]=sum((a[rows,meta['columns'].index(c)] for c in available),np.zeros(len(ids))) if len(available)==len(group) else np.full(len(ids),np.nan)
        finite=np.isfinite(sums['old'])&np.isfinite(sums['new']);d=sums['new']-sums['old'];record.update(valid_closed_chain_pairs=int(finite.sum()),inconsistent_with_isolated_decay_1e5=int((finite&(abs(d)>1e-5*abs(sums['old']))).sum()),**{'chain_delta_'+k:v for k,v in percentiles(d).items()});decay_rows.append(record)
    csvsave('tables/decay_chain_diagnostics.csv',decay_rows)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['all','ingest','statistics','consistency'],default='all');args=parser.parse_args()
    save('provenance/environment.json',dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,platform=platform.platform(),slurm_job_id=__import__('os').environ.get('SLURM_JOB_ID'),scripts={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}))
    if args.phase in ['all','ingest']:
        for label,path in INPUTS.items():ingest(label,path)
    if args.phase in ['all','statistics']:run_stats()
    if args.phase in ['all','consistency']:consistency()
    for label,path in INPUTS.items():assert stamp(path)==json.loads((ROOT/f'provenance/{label}_scan.json').read_text())['source']
    log('Population audit phase complete')
if __name__=='__main__':main()

