#!/usr/bin/env python3
"""Full-population density figures and all-species environmental summaries.

Every histogram uses all finite eligible particles; no random sampling or floors.
Missing values/zero denominators are excluded explicitly and counted in tables.
"""
import json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,TwoSlopeNorm
from common import *
plt.rcParams.update({'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':240,'figure.constrained_layout.use':True})

def finish(fig,name):
    p=ROOT/'figures'/name;p.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(p.with_suffix('.png'));fig.savefig(p.with_suffix('.pdf'));plt.close(fig)
def density(ax,x,y,xlabel,ylabel,title='',logx=True,logy=True,diagonal=False):
    valid=np.isfinite(x)&np.isfinite(y)
    if logx:valid &= x>0
    if logy:valid &= y>0
    xx=x[valid];yy=y[valid]
    if logx:xx=np.log10(xx)
    if logy:yy=np.log10(yy)
    if len(xx):
        h,xe,ye=np.histogram2d(xx,yy,bins=90)
        im=ax.pcolormesh(xe,ye,np.ma.masked_equal(h.T,0),norm=LogNorm(vmin=1,vmax=max(2,h.max())),cmap='viridis',shading='flat')
        if diagonal:
            lo=max(xx.min(),yy.min());hi=min(xx.max(),yy.max())
            if hi>lo:ax.plot([lo,hi],[lo,hi],color='tab:red',lw=.8,ls='--')
    else:ax.text(.5,.5,'No eligible pairs',ha='center',va='center',transform=ax.transAxes)
    ax.set(xlabel=('log₁₀ '+xlabel if logx else xlabel),ylabel=('log₁₀ '+ylabel if logy else ylabel),title=title+f'\nN = {len(xx):,} / {len(x):,}')
    return valid

def all_isotopes():
    summary=pd.read_csv(ROOT/'tables/species_statistics.csv')
    for p in sorted((ROOT/'tables/species').glob('*.npz')):
        with np.load(p) as d:
            fig,ax=plt.subplots(figsize=(5.5,4.5));density(ax,d['old'],d['new'],'old X','new X',p.stem+' mass fraction',diagonal=True)
            row=summary[summary.species==p.stem].iloc[0]
            fig.suptitle(f'Zero→positive: {int(row.newly_nonzero):,}; positive→zero: {int(row.newly_zero):,}',fontsize=8)
            finish(fig,'all_shared_isotopes/'+p.stem)
    summary_figures()

def summary_figures():
    summary=pd.read_csv(ROOT/'tables/species_statistics.csv')
    iso=summary[(summary.kind=='isotope')&(summary.finite_pairs>0)].copy()
    for start in range(0,len(iso),32):
        block=iso.iloc[start:start+32];data=np.array([block.exact/block.finite_pairs,block.significant/block.finite_pairs,block.newly_nonzero/block.finite_pairs,block.newly_zero/block.finite_pairs]).T
        fig,ax=plt.subplots(figsize=(7,max(4,len(block)*.21)));im=ax.pcolormesh(np.arange(5)-.5,np.arange(len(block)+1)-.5,data,vmin=0,vmax=1,cmap='magma',shading='flat');ax.invert_yaxis();ax.set_yticks(range(len(block)),block.species);ax.set_xticks(range(4),['Exact','Beyond 10⁻⁵ relative','Zero → positive','Positive → zero'],rotation=20,ha='right');fig.colorbar(im,ax=ax,label='Fraction of all matched particles');finish(fig,f'summaries/isotope_discrepancy_{start//32+1}')
    fig,ax=plt.subplots(figsize=(6,5));density(ax,iso.old_yield_shared_msun.to_numpy(),iso.new_yield_shared_msun.to_numpy(),'old yield (M☉)','new yield (M☉)','Integrated yields on the same particle IDs',diagonal=True)
    for c in ['26Al','26Mg','50Ti','54Cr','56Ni']:
        r=iso[iso.species==c]
        if len(r) and r.iloc[0].old_yield_shared_msun>0 and r.iloc[0].new_yield_shared_msun>0:
            rr=r.iloc[0];ax.annotate(c,(np.log10(rr.old_yield_shared_msun),np.log10(rr.new_yield_shared_msun)),xytext=(3,3),textcoords='offset points',fontsize=8)
    finish(fig,'summaries/shared_integrated_yields')
    fig,ax=plt.subplots(figsize=(13,4));x=np.arange(len(iso));ratio=iso.yield_ratio_shared.to_numpy();valid=np.isfinite(ratio)&(ratio>0);ax.bar(x[valid],np.log10(ratio[valid]),width=.8);ax.axhline(0,c='black',lw=.5);ax.set_xticks(x,iso.species,rotation=90,fontsize=7);ax.set(ylabel='log₁₀(new / old yield)',title='Matched-particle yield ratios; undefined ratios have no bar');finish(fig,'summaries/all_shared_yield_ratios')
    # Element sums use the SAME shared isotope set in old/new, never full totals.
    elements=iso.groupby('Z')[['old_yield_shared_msun','new_yield_shared_msun']].sum();fig,ax=plt.subplots(figsize=(11,4));x=np.arange(len(elements));ax.bar(x-.2,elements.old_yield_shared_msun,.4,label='Old');ax.bar(x+.2,elements.new_yield_shared_msun,.4,label='New');ax.set_yscale('log');ax.set_xticks(x,[SYMBOLS[int(z)] for z in elements.index],rotation=60);ax.set(ylabel='Shared-isotope subtotal (M☉)',title='Element comparison restricted to identical exported isotopes');ax.legend();finish(fig,'summaries/element_shared_subtotals');elements.to_csv(ROOT/'tables/element_shared_isotope_yields.csv')

def environments():
    meta,a,rawids=table('new');u,rows=unique_rows(rawids);cols=meta['columns'];mass=a[rows,cols.index('mass')]/MSUN
    fields={'peak temp':a[rows,cols.index('peak temp')]/1e9,'peak density':a[rows,cols.index('peak density')],'Y_{e}':a[rows,cols.index('Y_{e}')],'radius':np.sqrt(sum(a[rows,cols.index(c)]**2 for c in ['x','y','z'])),'mass':a[rows,cols.index('mass')]}
    edges={'peak temp':np.r_[0,1,2,3,4,5,7,10,np.inf],'peak density':np.r_[0,np.logspace(-8,10,19),np.inf],'Y_{e}':np.r_[0,.4,.45,.48,.49,.495,.5,.505,.51,.55,1,np.inf],'radius':np.r_[0,np.logspace(9,16,15),np.inf],'mass':np.r_[0,np.logspace(20,32,13),np.inf]}
    summaries=[];selected={'26Al','24Mg','25Mg','26Mg','40Ca','41Ca','44Ca','48Ca','44Ti','50Ti','54Cr','56Fe','60Fe','56Ni','88Sr','138Ba'}
    for c in cols:
        if species(c)[0]!='isotope':continue
        x=a[rows,cols.index(c)]
        for field,env in fields.items():
            for lo,hi in zip(edges[field][:-1],edges[field][1:]):
                ok=(env>=lo)&(env<hi)&np.isfinite(x)
                summaries.append(dict(species=c,field=field,lower=lo,upper=hi,particles=int(ok.sum()),positive=int((ok&(x>0)).sum()),yield_msun=float(np.sum(mass[ok]*x[ok])) if ok.any() else None,mean_X=float(np.mean(x[ok])) if ok.any() else None))
        if c in selected:
            fig,axes=plt.subplots(1,3,figsize=(12,3.8))
            for ax,field,label,lx in zip(axes,['peak temp','peak density','Y_{e}'],['legacy peak T (GK)','density at legacy peak T (g cm⁻³)','progenitor Yₑ'],[False,True,False]):density(ax,fields[field],x,label,c+' X',c,logx=lx)
            finish(fig,'thermodynamics/'+c)
    csvsave('tables/all_new_species_environmental_bins.csv',summaries)

def ratios_and_correlations():
    groups={'Mg_Al':['26Al_over_27Al','25Mg_over_24Mg','26Mg_over_24Mg'],'Ca':['41Ca_over_40Ca','42Ca_over_40Ca','43Ca_over_40Ca','44Ca_over_40Ca','46Ca_over_40Ca','48Ca_over_40Ca'],'Ti':['46Ti_over_48Ti','47Ti_over_48Ti','49Ti_over_48Ti','50Ti_over_48Ti'],'Cr':['50Cr_over_52Cr','53Cr_over_52Cr','54Cr_over_52Cr'],'Fe':['54Fe_over_56Fe','57Fe_over_56Fe','58Fe_over_56Fe'],'Ni':['58Ni_over_56Ni','60Ni_over_56Ni','61Ni_over_56Ni','62Ni_over_56Ni','64Ni_over_56Ni']}
    for group,names in groups.items():
        names=[n for n in names if (ROOT/'tables/ratios'/f'{n}.npz').exists()]
        if not names:continue
        fig,axes=plt.subplots(math.ceil(len(names)/3),3,figsize=(12,math.ceil(len(names)/3)*3.6),squeeze=False)
        for ax,name in zip(axes.ravel(),names):
            with np.load(ROOT/'tables/ratios'/f'{name}.npz') as d:density(ax,d['old'],d['new'],'old number ratio','new number ratio',name.replace('_over_',' / '),diagonal=True)
        for ax in axes.ravel()[len(names):]:ax.set_visible(False)
        finish(fig,'isotope_ratios/'+group)
    pairs=[('48Ca_over_40Ca','50Ti_over_48Ti'),('50Ti_over_48Ti','54Cr_over_52Cr'),('48Ca_over_40Ca','54Cr_over_52Cr'),('50Cr_over_52Cr','54Cr_over_52Cr'),('53Cr_over_52Cr','54Cr_over_52Cr'),('26Al_over_27Al','26Mg_over_24Mg'),('25Mg_over_24Mg','26Mg_over_24Mg')]
    for xname,yname in pairs:
        px=ROOT/'tables/ratios'/f'{xname}.npz';py=ROOT/'tables/ratios'/f'{yname}.npz'
        if not px.exists() or not py.exists():continue
        with np.load(px) as dx,np.load(py) as dy:
            fig,axes=plt.subplots(1,2,figsize=(10,4))
            for ax,label in zip(axes,['old','new']):density(ax,dx[label],dy[label],xname.replace('_over_','/'),yname.replace('_over_','/'),label+' ejecta ratios')
            finish(fig,'correlations/'+xname+'__'+yname)

def closure_figures():
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for label,color in [('old','tab:blue'),('new','tab:orange')]:
        with np.load(ROOT/f'tables/{label}_closure.npz') as d:
            values=d['selected_isotope_sum'];axes[0].hist(values,bins=np.linspace(0,2.1,150),histtype='step',label=label,color=color);res=values-1;axes[1].hist(res,bins=np.linspace(-.02,.02,150),histtype='step',label=label,color=color)
            f,ax=plt.subplots(figsize=(5.5,4.5));density(ax,d['reported_progenitor_ye'],d['partial_composition_ye'],'reported progenitor Yₑ','available-isotope Σ(Z/A)X',label+'; partial inventory, differing epochs',False,False,True);finish(f,'validation/'+label+'_partial_ye')
    for ax in axes:ax.set_yscale('log');ax.set_ylabel('Particles');ax.legend()
    axes[0].set_xlabel('Available-isotope ΣX');axes[1].set_xlabel('ΣX − 1 (central range only)');finish(fig,'validation/mass_fraction_closure')
    p=ROOT/'tables/chunk_integrity.csv'
    if p.exists():
        c=pd.read_csv(p);good=c[c.status=='complete'];fig,axes=plt.subplots(2,1,figsize=(11,6),sharex=True);axes[0].scatter(good.chunk,good.last_sum_max,s=5);axes[0].axhline(1,c='black',lw=.5);axes[0].set_ylabel('Maximum last ΣX');axes[1].plot(good.chunk,good.plot_exact_mismatch_cells,'.',ms=2);axes[1].set(xlabel='Source chunk number',ylabel='Plotting mismatch cells');finish(fig,'validation/chunk_integrity')

def main():
    all_isotopes();log('All shared isotope figures complete');environments();ratios_and_correlations();closure_figures()
    save('provenance/figure_methods.json',dict(population='all valid unique matched IDs; all eligible positive pairs in logarithmic plots',density_bins=90,log_floor=None,zero_missing_policy='Not plotted in log coordinates; counts remain in statistics',png_dpi=240,pdf='vector pcolormesh and text',physical_fields='Copied legacy progenitor Ye and density at early-snapshot peak temperature; not independent validation of PRISM trajectories'))
    log('Figures complete')
if __name__=='__main__':main()
