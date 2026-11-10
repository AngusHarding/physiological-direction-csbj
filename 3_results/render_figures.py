# =============================================================================
# READ ME FIRST: render_figures.py
# -----------------------------------------------------------------------------
# WHAT IT DOES  Draws Figures 1 and 2 in their published layout from the
#               numbers the Results 1 and 2 scripts hand it.
#
# METHODS       None. Drawing only; no fitting, preprocessing or uncertainty
#               calculation happens here.
# =============================================================================
"""Approved Figure 1/2 drawing routines, 31 October 2026.

These functions only draw supplied values. Inputs may be the pinned display
records used for the accepted figures, or newly computed pipeline results.
No fitting, preprocessing, uncertainty calculation or Word editing occurs.

Figure 1: 244 curve/interval rows and ten forest points/intervals; Figure 2:
15,971 participant coordinates and 14 forest points/intervals in the accepted
version. Populations, endpoints and interval interpretation remain those in
the paper captions and source analysis notes. Rebuilding from saved records
is display reproduction, not a fresh analysis or independent validation.

Keep local sex legends inside a/b; forest-era legends belong beneath their
panel headings, above row labels. Use adverse distance with subscript A.
The selected output directory is the only write location. No archive script
or revision helper is imported. See README.md and INPUT_MANIFEST.json.
"""
from pathlib import Path
import csv,json,math,random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator,FuncFormatter,NullLocator

INK,MUTED,GRID='#222222','#666666','#E5E5E5'
MEN,WOMEN,GREEN,PURPLE='#0072B2','#E69F00','#009E73','#5B4B9A'
ERAS={'modern':('NHANES 1999–2018',11497,2494),'iii':('NHANES III',4474,3401)}
RECIPES=('SS','SC','PC')
LABELS={'SS':'Sex-specific centre\nand SD (SS)',
    'SC':'Sex-specific centre,\ncommon SD (SC)','PC':'Pooled centre\nand SD (PC)'}
LONG={'SS':'Sex-specific centre and SD (SS)',
    'SC':'Sex-specific centre, common SD (SC)','PC':'Pooled centre and SD (PC)'}
ERA_LABEL={'modern':ERAS['modern'][0],'iii':ERAS['iii'][0]}



def write_csv(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def render_figure1(curves, samples, summary, output, source_paths=None):
    """Draw supplied values verbatim and record each plotted source row."""
    source_paths = source_paths or {}
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
        'axes.titlesize':9,'xtick.labelsize':7.5,'ytick.labelsize':7.5,
        'axes.edgecolor':MUTED,'axes.linewidth':.65,'text.color':INK,
        'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,
        'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none',
        'pdf.fonttype':42})
    fig = plt.figure(figsize=(185/25.4,165/25.4), facecolor='white')
    lower_bottom = .16 - 2 / 165  # 2 mm lower, with plot sizes unchanged.
    axes = dict(a=fig.add_axes([.105,.57,.35,.31]), b=fig.add_axes([.635,.57,.35,.31]),
                c=fig.add_axes([.22,lower_bottom,.235,.235]),
                d=fig.add_axes([.75,lower_bottom,.235,.235]))
    evidence = []
    for panel, era, title in [('a','NHANES 1999-2018','NHANES 1999–2008'),('b','NHANES III','NHANES III')]:
        ax = axes[panel]
        for sex, colour in [('Men',MEN),('Women',WOMEN)]:
            rows = [r for r in curves if r['era']==era and r['sex']==sex]
            x = [float(r['age']) for r in rows]
            y, lo, hi = [[float(r[k]) for r in rows] for k in ('probability','ci_low','ci_high')]
            ax.fill_between(x, lo, hi, color=colour, alpha=.17, lw=0)
            ax.plot(x, y, color=colour, lw=1.8)
            for r in rows:
                evidence.append(dict(panel=panel,era=era,label=sex,age=r['age'],
                    estimate=r['probability'],low=r['ci_low'],high=r['ci_high'],
                    units='probability',source=source_paths.get('curves','revision/notes/S1_sex_gap_two_era/age_curves.csv'),
                    source_key=f"{era}|{sex}|{r['age']}"))
        ax.axvspan(52,79,color=INK,alpha=.045,lw=0)
        ax.set(xlim=(17,80),ylim=(.0008,1),yscale='log',xlabel='Age at examination (years)' if panel=='a' else 'Age at interview (years)')
        ax.yaxis.set_major_locator(FixedLocator([.001,.003,.01,.03,.1,.3,1]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{100*v:g}%'))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_xticks([20,40,60,80])
        ax.grid(axis='y',color=GRID,lw=.5)
        ax.set_axisbelow(True)
        ax.text(65.5,.025,'Analysis ages\n52–79',transform=ax.get_xaxis_transform(),
            ha='center',va='bottom',fontsize=6.4,color=MUTED)
        ax.set_ylabel('Ten-year mortality probability')
        counts={r['sex']:r for r in samples if r['era']==era}
        n=sum(int(counts[s]['n']) for s in counts)
        events=sum(int(counts[s]['deaths_within_10y']) for s in counts)
        ax.set_title(title,loc='left',pad=19,fontweight='bold')
        ages='20–79' if panel=='a' else '18–79'
        ax.text(0,1.015,f'Ages {ages} · n = {n:,} · {events:,} deaths',
            transform=ax.transAxes,fontsize=6.7,color=MUTED,va='bottom')
        ax.legend(handles=[Line2D([0],[0],color=MEN,lw=2,label='Men'),
                           Line2D([0],[0],color=WOMEN,lw=2,label='Women')],
            loc='upper left',bbox_to_anchor=(.025,.975),borderaxespad=0,
            frameon=False,ncol=1,fontsize=7.6,handlelength=2,labelspacing=.35)

    def point(panel, era, y, row, colour, metric):
        ax=axes[panel]
        estimate, low, high = (float(row[k]) for k in
            (('observed','ci_low','ci_high') if metric=='contrast' else ('hr','hr_low','hr_high')))
        marker='o' if era=='modern' else 'D'
        ax.hlines(y,low,high,color=colour,lw=1.35)
        ax.plot(estimate,y,marker=marker,markersize=5,markeredgewidth=1.2,
            markeredgecolor=colour,markerfacecolor=colour if era=='modern' else 'white')
        evidence.append(dict(panel=panel,era=ERA_LABEL[era],label=row['measure'],age='',
            estimate=repr(estimate),low=repr(low),high=repr(high),
            units='cohort SD' if metric=='contrast' else 'hazard ratio',
            source=source_paths.get('summary','revision/notes/S1b_magnitude_two_era/bootstrap_summary.csv'),
            source_key=f"{era}|{row['measure']}|{row['quantity']}"))
    for measure,y,colour in [('FI_lin',1,GREEN),('D',0,PURPLE)]:
        for era, off in [('modern',.13),('iii',-.13)]:
            row=next(r for r in summary if r['era']==era and r['measure']==measure and r['quantity']=='contrast')
            point('c',era,y+off,row,colour,'contrast')
    ax=axes['c'];ax.set(yticks=[1,0],yticklabels=['Adverse distance\n($D_A$)','Mahalanobis\ndistance ($D_M$)'],
        xlim=(-.18,.28),ylim=(-.55,1.55),xlabel='Women minus men\n(cohort standard deviations)')
    ax.set_xticks([-.1,0,.1,.2]);ax.axvline(0,color=MUTED,lw=.7,ls=(0,(3,2)))
    for measure,y,colour in [('none',2,MUTED),('FI_lin',1,GREEN),('D',0,PURPLE)]:
        for era, off in [('modern',.14),('iii',-.14)]:
            quantity='b0_sex_only' if measure=='none' else 'b1'
            row=next(r for r in summary if r['era']==era and r['measure']==measure and r['quantity']==quantity)
            point('d',era,y+off,row,colour,'HR')
    ax=axes['d'];ax.set(yticks=[2,1,0],yticklabels=['Sex only','+ adverse distance\n($D_A$)','+ Mahalanobis distance\n($D_M$)'],
        xscale='log',xlim=(.47,1.035),ylim=(-.55,2.55),xlabel='Female-to-male hazard ratio\n(log scale)')
    ax.xaxis.set_major_locator(FixedLocator([.5,.6,.7,.8,1]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{v:g}'))
    ax.xaxis.set_minor_locator(NullLocator());ax.axvline(1,color=MUTED,lw=.7,ls=(0,(3,2)))
    for key in ('c','d'):
        axes[key].grid(axis='x',color=GRID,lw=.5);axes[key].set_axisbelow(True)
        axes[key].tick_params(axis='y',length=0,pad=6)
    for key,x,y,title in [('a',.025,.955,''),('b',.555,.955,''),
            ('c',.025,.445,'Distance from young-healthy reference'),('d',.555,.445,'Mortality after adjustment')]:
        fig.text(x,y,key,fontweight='bold',fontsize=10,va='center')
        if title:fig.text(x+.033,y,title,fontsize=8.6,fontweight='bold',va='center')
    handles=[Line2D([0],[0],marker='o',lw=0,color=INK,ms=5,label='NHANES 1999–2018'),
        Line2D([0],[0],marker='D',lw=0,color=INK,markerfacecolor='white',ms=5,label='NHANES III')]
    for key in ('c','d'):
        axes[key].legend(handles=handles,loc='upper left',
            bbox_to_anchor=(.058 if key=='c' else .588,.423),bbox_transform=fig.transFigure,
            borderaxespad=0,borderpad=0,ncol=1,frameon=False,fontsize=7.2,
            handlelength=1.2,handletextpad=.6,labelspacing=.35)
    # Record actual rendered legend extents, so the verifier can check that
    # each key belongs to its own panel and clears the data/title regions.
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    layout={}
    for key,ax in axes.items():
        legend=ax.get_legend()
        extent=legend.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        layout[key]=dict(axes=list(ax.get_position().bounds),legend=list(extent.bounds),
            labels=[t.get_text() for t in legend.get_texts()],
            legend_location='inside upper left' if key in ('a','b') else 'beneath heading above row labels')
    (output/'LAYOUT_CHECKS.json').write_text(json.dumps(layout,indent=2)+'\n',encoding='utf-8')
    fig.savefig(output/'figure1.svg',facecolor='white',metadata={'Date':None})
    # Matplotlib emits trailing spaces inside multiline SVG path attributes;
    # newline remains valid path whitespace, so strip those spaces for clean diffs.
    svg=output/'figure1.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    fig.savefig(output/'figure1.png',dpi=600,facecolor='white')
    fig.savefig(output/'figure1_preview.png',dpi=220,facecolor='white')
    plt.close(fig)
    write_csv(output/'figure1_source_data.csv',evidence)
    return dict(curve_rows=sum(r['panel'] in ['a','b'] for r in evidence),
        forest_points=sum(r['panel'] in ['c','d'] for r in evidence),
        unchanged_original_panels={'a':'a','b':'b','f':'c','g':'d'},
        omitted_original_panels=['c: causes','d: condition differences','e: subgroups'])


def render_figure2(participants, ledger, baseline, output, source_paths=None):
    source_paths = source_paths or {}
    points=[]
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
        'axes.titlesize':9,'xtick.labelsize':7.3,'ytick.labelsize':7.2,
        'axes.edgecolor':MUTED,'axes.linewidth':.65,'text.color':INK,
        'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,
        'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(185/25.4,185/25.4),facecolor='white')
    axes=dict(a=fig.add_axes([.105,.60,.35,.29]),b=fig.add_axes([.635,.60,.35,.29]),
        c=fig.add_axes([.235,.19,.22,.225]),d=fig.add_axes([.765,.19,.22,.225]))
    for panel,era in [('a','modern'),('b','iii')]:
        ax=axes[panel];rows=[r for r in participants if r['cohort']==era]
        random.Random(20261024).shuffle(rows)  # all points, mixed drawing order only
        x=[r['continuous_departure_index'] for r in rows];y=[r['q'] for r in rows]
        assert min(x)>.25 and max(x)<2.8 and min(y)>-2.8 and max(y)<10.5
        ax.scatter(x,y,s=6,c=[MEN if r['sex']==1 else WOMEN for r in rows],
            alpha=.28,linewidths=0,marker='.')
        ax.set(xlim=(.25,2.8),ylim=(-2.8,10.5),
            xlabel='Adverse distance, $D_A$\n(era-specific units)',
            ylabel='Mortality score q\n(era-specific units)')
        ax.set_xticks([.5,1,1.5,2,2.5]);ax.set_yticks([-2,0,2,4,6,8,10])
        ax.grid(color=GRID,lw=.45);ax.set_axisbelow(True)
        ax.set_title(ERAS[era][0],loc='left',fontweight='bold',pad=19)
        n,events=ERAS[era][1:]
        ax.text(0,1.015,f'Ages 52–79 · n = {n:,} · {events:,} deaths',
            transform=ax.transAxes,fontsize=6.7,color=MUTED,va='bottom')
        # Repeat the stacked sex legend inside each scatter panel, in the
        # sparse upper-left area rather than sharing a remote figure legend.
        ax.legend(handles=[Line2D([0],[0],marker='o',lw=0,color=MEN,label='Men',ms=4),
            Line2D([0],[0],marker='o',lw=0,color=WOMEN,label='Women',ms=4)],
            loc='upper left',bbox_to_anchor=(.015,.985),ncol=1,fontsize=8,
            frameon=True,facecolor='white',edgecolor='none',framealpha=.9,
            borderpad=.25,labelspacing=.3,handletextpad=.45)

    def point(panel,era,recipe,y,metric):
        r=ledger[(era,recipe,metric)]
        transform=math.exp if metric in ['b0','b1_q'] else float
        estimate,low,high=[transform(float(r[k])) for k in ['estimate','low','high']]
        colour=MUTED if metric=='b0' else PURPLE
        ax=axes[panel];marker='o' if era=='modern' else 'D'
        ax.hlines(y,low,high,color=colour,lw=1.35)
        ax.plot(estimate,y,marker=marker,ms=5,mew=1.2,mec=colour,
            mfc=colour if era=='modern' else 'white')
        points.append(dict(panel=panel,cohort=era,recipe=recipe,metric=metric,
            estimate=estimate,low=low,high=high,
            units='hazard ratio' if metric in ['b0','b1_q'] else 'cohort SD',
            source_estimate=r['estimate'],source_low=r['low'],source_high=r['high'],
            transform='exp' if metric in ['b0','b1_q'] else 'identity',
            source=source_paths.get('ledger','revision/text/common_scaling_20261022/RESULT_LEDGER.csv')))
    for recipe,y in [('SS',2),('SC',1),('PC',0)]:
        for era,offset in [('modern',.115),('iii',-.115)]:
            point('c',era,recipe,y+offset,'contrast_q')
            point('d',era,recipe,y+offset,'b1_q')
    for era,offset in [('modern',.115),('iii',-.115)]:
        r=ledger[(era,'SS','b0')]
        original=next(x for x in baseline if x['era']==era and x['quantity']=='b0_sex_only')
        assert abs(math.exp(float(r['estimate']))-float(original['hr']))<1e-6
        point('d',era,'SS',3+offset,'b0')
    ax=axes['c'];ax.set(yticks=[2,1,0],yticklabels=[LABELS[x] for x in RECIPES],
        xlim=(-.85,.15),ylim=(-.5,2.5),xlabel='Women minus men in q\n(cohort standard deviations)')
    ax.set_xticks([-.8,-.4,0]);ax.axvline(0,color=MUTED,lw=.7,ls=(0,(3,2)))
    ax=axes['d'];ax.set(yticks=[3,2,1,0],
        yticklabels=['Before q\nadjustment']+[LABELS[x] for x in RECIPES],
        xscale='log',xlim=(.46,1.10),ylim=(-.55,3.55),
        xlabel='Female-to-male hazard ratio\n(log scale)')
    ax.xaxis.set_major_locator(FixedLocator([.5,.6,.8,1]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{v:g}'))
    ax.xaxis.set_minor_locator(NullLocator());ax.axvline(1,color=MUTED,lw=.7,ls=(0,(3,2)))
    for key in ['c','d']:
        axes[key].grid(axis='x',color=GRID,lw=.5);axes[key].set_axisbelow(True)
        axes[key].tick_params(axis='y',length=0,pad=6)
    for key,x,y,title in [('a',.025,.955,''),('b',.555,.955,''),
        ('c',.025,.49,'Reference choice and q'),('d',.555,.49,'Mortality after q adjustment')]:
        fig.text(x,y,key,fontweight='bold',fontsize=10,va='center')
        if title:fig.text(x+.033,y,title,fontsize=8.5,fontweight='bold',va='center')
    fig.text(.5,.095,'SS → SC changes the scale; SC → PC changes the centre.',
        ha='center',fontsize=7.1,color=MUTED)
    # Place a local era legend above each forest's left-hand row labels.
    # The reserved strip separates it from the title and first data row.
    for x in (.058,.588):
        fig.legend(handles=[Line2D([0],[0],marker='o',lw=0,color=INK,label='NHANES 1999–2018',ms=5),
            Line2D([0],[0],marker='D',lw=0,color=INK,mfc='white',label='NHANES III',ms=5)],
            loc='upper left',bbox_to_anchor=(x,.468),frameon=False,ncol=1,
            fontsize=7.4,borderaxespad=0,labelspacing=.35,handletextpad=.55)
    fig.savefig(output/'figure2.svg',facecolor='white',metadata={'Date':None})
    path=output/'figure2.svg';path.write_text('\n'.join(x.rstrip() for x in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    fig.savefig(output/'figure2.png',dpi=600,facecolor='white')
    fig.savefig(output/'figure2_preview.png',dpi=220,facecolor='white');plt.close(fig)
    write_csv(output/'figure2_forest_source_data.csv',points)
    assert len(points)==14
