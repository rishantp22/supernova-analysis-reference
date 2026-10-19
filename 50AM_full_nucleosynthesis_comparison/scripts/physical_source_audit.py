#!/usr/bin/env python3
"""Independent complete-population physical-field and initial-composition check."""
import sys
from pathlib import Path
from common import *
sys.path.insert(0,'/scratch/rpraka15/supernova-grain-analysis/code')
from readers import open_sdf
BASE=Path('/data/pyoung3/supernova_data/50Am/run3g_50Am6')
def main():
    om,ot,oi=table('old');ou,orr=unique_rows(oi);rows=[];summary={}
    initial=pd.read_csv(ROOT/'tables/initial_composition_population.csv.gz').set_index('id')
    for label,suffix in [('initial','0000'),('final','2000')]:
        source=BASE/('run3g_50Am6_sph.'+suffix);before=stamp(source);h,a=open_sdf(source);ids=np.array(a['ident'],dtype='i8');u,rr=unique_rows(ids);common,ri,pi=np.intersect1d(u,ou,return_indices=True);ar=rr[ri];pr=orr[pi]
        factors={'x':6.955e10,'y':6.955e10,'z':6.955e10,'vx':6.955e8,'vy':6.955e8,'vz':6.955e8,'ax':6.955e6,'ay':6.955e6,'az':6.955e6,'mass':1.9889e27,'h':6.955e10,'density':1.9889e27/6.955e10**3}
        for name,cf in (factors.items() if label=='final' else [('Y_{e}',1)]):
            field='Y_el' if name=='Y_{e}' else ('rho' if name=='density' else name);x=np.asarray(a[field][ar],dtype='f8')*cf;y=ot[pr,om['columns'].index(name)];good=np.isfinite(x)&np.isfinite(y);match=good&(abs(x-y)<=1e-5*abs(x));rows.append(dict(dataset=label,field=name,pairs=len(x),matches_export_precision=int(match.sum()),mismatches=int((good&~match).sum()),source_zero_plot_nonzero=int(((x==0)&(y!=0)).sum()),source_nonzero_plot_zero=int(((x!=0)&(y==0)).sum()),max_abs_delta=float(np.max(abs(x[good]-y[good])))))
        if label=='initial':
            # Native p fields are Z, m fields are N. Source unburned.c explicitly
            # excludes slots 21 and 22, which repeat free proton/neutron slots.
            Z=np.column_stack([a[f'p{k}'] for k in range(1,21)]);N=np.column_stack([a[f'm{k}'] for k in range(1,21)]);X=np.column_stack([a[f'f{k}'] for k in range(1,21)]).astype('f8');A=Z+N;total=X.sum(axis=1);ye=np.sum(X*Z/A,axis=1)
            np.savez_compressed(ROOT/'tables/hydro_initial_composition.npz',id=ids,Z=Z,A=A,X=X,sumX=total,composition_Ye=ye,reported_Ye=a['Y_el'])
            seed=initial.reindex(ids);sye=seed.Ye_after_PRISM_initial_normalization.to_numpy();diff=sye-np.asarray(a['Y_el']);summary.update(native_initial_rows=len(ids),native_initial_sumX=percentiles(total),native_composition_minus_Ye=percentiles(ye-a['Y_el']),archived_seed_minus_native_initial_Ye=percentiles(diff),archived_seed_comparable=int(np.isfinite(diff).sum()),archived_seed_matches_native_initial_Ye_1e5=int((np.isfinite(diff)&(abs(diff)<=1e-5)).sum()),seed_epoch_caveat='Different seed composition is verified; the intended seed epoch/mapping cannot be inferred from this difference alone')
        assert stamp(source)==before
        summary[label+'_source']=before
    csvsave('tables/plotting_vs_original_hydro.csv',rows);save('provenance/physical_source_audit.json',summary);log('Independent physical source audit complete')
if __name__=='__main__':main()
