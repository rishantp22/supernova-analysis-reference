#!/usr/bin/env python3
"""Compare every archived PRISM seed with the native initial SNSPH composition."""
import re,zipfile,sys
from common import *
from population_audit import stats
sys.path.insert(0,'/scratch/rpraka15/supernova-grain-analysis/code')
from readers import time_blocks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def main():
    h=np.load(ROOT/'tables/hydro_initial_composition.npz');hi=h['id'];hr={int(x):i for i,x in enumerate(hi)};hx=h['X'];hz=h['Z'];ha=h['A'];hye=h['composition_Ye'];assert (hz==hz[0]).all() and (ha==ha[0]).all()
    inv=pd.read_csv(ROOT/'tables/initial_species_inventory.csv');codes=(inv.Z*512+inv.A).to_numpy();nm,nt,ni=table('new');initial=pd.read_csv(ROOT/'tables/initial_composition_population.csv.gz');ids=initial.id.to_numpy();rows={int(v):i for i,v in enumerate(ids)}
    arr=np.lib.format.open_memmap(ROOT/'tables/prism_initial_X.npy',mode='w+',dtype='f8',shape=(len(ids),len(codes)));arr[:]=np.nan;np.save(ROOT/'tables/prism_initial_ids.npy',ids)
    p=Path('/scratch/svierlin/prism-v02/input/conditions/50Am_abund.zip');before=stamp(p);count=0
    with zipfile.ZipFile(p) as z:
        for member in z.infolist():
            if '__MACOSX' in member.filename:continue
            m=re.fullmatch(r'(\d+)_abund\.txt',Path(member.filename).name)
            if not m:continue
            ident=int(m[1]);v=np.fromstring(z.read(member).replace(b'D',b'E'),sep=' ').reshape(-1,3);code=v[:,0].astype('i8')*512+v[:,1].astype('i8');j=np.searchsorted(codes,code);assert np.array_equal(codes[j],code);arr[rows[ident],j]=v[:,2];count+=1
            if count%100000==0:log(f'Compared seed inventory: {count:,}')
    arr.flush();assert stamp(p)==before and count==len(ids)
    common,ii,jj=np.intersect1d(ids,hi,return_indices=True);total=initial.sumX.to_numpy()[ii];seed_ye=initial.Ye_after_PRISM_initial_normalization.to_numpy()[ii];diff=seed_ye-hye[jj];out=[];mass=np.ones(len(common));mismatch=np.zeros(len(common),dtype='i4');missing=np.zeros(len(common),dtype='i4')
    for k,(Z,A) in enumerate(zip(hz[0],ha[0])):
        code=int(Z)*512+int(A);name='n' if code==1 else str(A)+SYMBOLS[Z];x=hx[jj,k];y=np.full(len(common),np.nan)
        if code in codes:y=np.asarray(arr[ii,np.searchsorted(codes,code)])/total
        s,_=stats(x,y,common,mass,mass,1e-5);s={q:v for q,v in s.items() if 'yield' not in q};s.update(species=name,Z=int(Z),A=int(A));out.append(s);mismatch+=np.isfinite(y)&(abs(y-x)>1e-5*abs(x));missing+=~np.isfinite(y)
        if name in ['1H','4He','12C','16O','24Mg','28Si']:
            ok=(x>0)&(y>0)&np.isfinite(x)&np.isfinite(y);fig,ax=plt.subplots(figsize=(5.5,4.5))
            if ok.any():ax.hexbin(np.log10(x[ok]),np.log10(y[ok]),gridsize=90,bins='log',mincnt=1)
            ax.set(xlabel='Native initial log10 X('+name+')',ylabel='Normalized PRISM seed log10 X('+name+')',title=f'{name}: {len(common):,} ID-matched inputs; positive pairs shown')
            for ext in ['png','pdf']:fig.savefig(ROOT/f'figures/validation/seed_vs_native_{name}.{ext}',dpi=240,bbox_inches='tight')
            plt.close(fig)
    csvsave('tables/prism_seed_vs_native_species.csv',out);csvsave('tables/prism_seed_vs_native_particles.csv.gz',dict(id=common,mismatched_shared_species=mismatch,missing_shared_species=missing,native_composition_Ye=hye[jj],normalized_seed_Ye=seed_ye))
    examples=[]
    for p in sorted((ROOT/'provenance/raw_archive_prefixes').glob('*/*/xtime.txt')):
        ident=int(p.parent.name)
        if ident not in rows or ident not in hr:continue
        with p.open() as f:header,first=next(time_blocks(f))
        seed=arr[rows[ident]];sx=np.nansum(seed);native_index=hr[ident]
        for code,name in [(513,'1H'),(1028,'4He'),(3084,'12C'),(4112,'16O')]:
            Z,A=code//512,code%512;k=np.flatnonzero((hz[0]==Z)&(ha[0]==A));j=np.searchsorted(codes,code);normalized=seed[j]/sx if j<len(codes) and codes[j]==code else np.nan
            examples.append(dict(id=ident,species=name,native_initial_X=hx[native_index,k[0]] if len(k) else np.nan,seed_normalized_X=normalized,first_runtime_X=first.get((Z,A),np.nan),runtime_time=header['time_s'],runtime_T_GK=header['temp_GK']))
    csvsave('tables/seed_runtime_native_examples.csv',examples)
    save('provenance/seed_attribution_summary.json',dict(ids=len(common),input_species=len(codes),native_species=len(hz[0]),particles_any_shared_isotope_difference=int((mismatch>0).sum()),seed_minus_native_composition_Ye=percentiles(diff),composition_Ye_difference_gt_1e5=int((abs(diff)>1e-5).sum()),seed_vs_native_approximation_caveat='Native SNSPH has an approximate 20-species inventory. Small species differences need not be errors; large H/He-to-O inventory differences require seed provenance and cannot be explained by exporting more species.'))
    log('Complete seed-to-native initial composition comparison finished')
if __name__=='__main__':main()
