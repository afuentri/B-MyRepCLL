import os
import sys
import shutil
import matplotlib
import pandas as pd
import xlsxwriter
import glob
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from chord import Chord
import coverage_IGHs
sns.set(style="whitegrid")

# define tables to use (given as argument)
input_table = sys.argv[1] ## homology_table.csv (output IGH pipeline)
## include 4 digits in run ID (example: RUN0193)
run =  sys.argv[2] ## example --> RUN0327
## out dir
out = sys.argv[3]
reads = sys.argv[4]
resume_Vmapping = sys.argv[5]
mincov = sys.argv[6]

def clone_plotting(df):

    s = df['sample_name'].unique().tolist()
    for e in range(0, len(s), 20):
        nsamples = e + 19
        name = 'clones_{}-{}'.format(e+1, nsamples + 1)
        plot_clones(df[df['sample_name'].isin(s[e:nsamples])], name)


def plot_clones(df_clonal, n):
    sns.set(font_scale=4)
    sns.set_style("white")
    sns.despine()
    fig = plt.figure(figsize=(250,75))
    pal1=['#3fd4ca', '#f0eb71']
    pal2=['#191a19','#c90a10']

    ax2 = sns.barplot(x="rearrangement", y="invert_diff", hue="max_diff", dodge=False,
                      data=df_clonal.sort_values(by=['sample_name', 
                                                                'percent_reads_mapped_Vgene']), 
                     palette=pal2)
    #for p in ax2.patches:## uncomment for maxdiff value annotation
    #    ax2.annotate(format(p.get_height(), '.1f'), 
    #                   (p.get_x() + p.get_width() / 2., p.get_height() -3), 
    #                   ha = 'center', va = 'center', 
    #                   xytext = (0, 9), size=22,
    #                   textcoords = 'offset points')

    ax2.get_legend().remove()  
    ax = sns.barplot(x="rearrangement", y="percent_reads_mapped_Vgene",data=df_clonal.sort_values(by=['sample_name', 
                                                                    'percent_reads_mapped_Vgene']))
    # specify the lines and labels of the first legend
    #ax.get_legend().remove()

    # Create the second legend and add the artist manually.
    from matplotlib.legend import Legend
    handles, labels = ax.get_legend_handles_labels()

    leg2 = Legend(ax, handles[:2],labels=['normal diff', 'max diff'],
                 loc='upper right', frameon=False)

    ax.add_artist(leg2)

    name = '{}.png'.format(n)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(name)
    plt.gcf().clear()


def predicted_status(df):

    samples = df['sample_name'].unique().tolist()

    for e in samples:
        
        df2 = df[df['sample_name'] == e].sort_values(by='percent_reads_mapped_Vgene', 
                                                                     ascending=False)
        rear = df2['rearrangement'].tolist()
    
        for i in range(0, len(rear)):
            
            if i == (len(rear)-1):
                t = df[df['rearrangement'] == rear[i]].index[0]
                n = float(df.at[t, 'percent_reads_mapped_Vgene']/1) 
            else:
                follow = rear[i + 1]
                t = df[df['rearrangement'] == rear[i]].index[0]
                t2 = df[df['rearrangement'] == rear[i + 1]].index[0]
                n = float(df.at[t, 'percent_reads_mapped_Vgene'])/float(df.at[t2, 'percent_reads_mapped_Vgene'])
            
            
            df.at[t, 'diff'] = n

    df['invert_diff'] = -(df['diff'])

    maxdiff = df.loc[df.groupby(['sample_name'])['diff'].idxmax()]
    list_max_diff = maxdiff['rearrangement'].unique().tolist()

    df['max_diff'] = np.where(df['rearrangement'].isin(list_max_diff), 'YES', 'NO')

    df['rearrangement'] = df['rearrangement'].str.split('_').str[1]
    df['predicted_status'] = 'polyclonal'
    df['clone_status'] = 'SUBCLONAL'
    s = df['sample_name'].unique().tolist()
    for e in s:
        
        n = df[df['sample_name'] == e]
        if float(n[n['max_diff'] =='YES']['diff'].to_string().split()[1]) > 5:
        
            r = float(n[n['max_diff'] =='YES']['rearrangement'].to_string().split()[1])
            
            df['clone_status'] = np.where((df['sample_name'] == e) & (df['rearrangement'].astype('float') <= r),
                                          'CLONAL', df['clone_status'])
            
            df['predicted_status'] = np.where(df['sample_name'] == e, 
                                                      str(int(r)) + 'CLONE', df['predicted_status']) 
            
    return df


def cdr3chord(df):

    """"""
    dict_cdr3 = {}
    cdr3s = df['CDR3'].unique().tolist()
    sa = df['sample_name'].unique().tolist()
    for t in sa:
    
        dict_cdr3[t] = df[df['sample_name'] == t]['CDR3'].unique().tolist()
        if 'None' in dict_cdr3[t]:
            dict_cdr3[t].remove('None')
    
    asmatrix = pd.DataFrame(index=sa, columns=sa)
    
    for e in dict_cdr3:
        for f in dict_cdr3:
            if e == f:
                inters = 0
            else:
                inters = len(set(dict_cdr3[e]).intersection(set(dict_cdr3[f])))
            asmatrix.loc[e][e] = inters
            asmatrix.loc[f][e] = inters
    
    asmatrix = asmatrix[asmatrix.sum() > 0]
    
    asmatrix_small = asmatrix[asmatrix.index.tolist()]
       
    names = [n.replace('IGH-','') for n in asmatrix_small.columns.tolist()]
    #names = [r for r in range(0, len(asmatrix_small.columns.tolist()))]
    Chord(asmatrix_small.as_matrix().tolist(), names,
          font_size_large="7px").to_html("chord.html")
    
        
def homology_resume(hom, rescued, outtable, outcoverage, min_cov, 
                    cov_table, folder_consensus,
                    cdr3='default'):


    """"""
    first_rearrangement = outtable.replace('.xlsx', '_principal-rearrangement.csv')
    top_rearrangement = outtable.replace('.xlsx', '_top3-rearrangement.csv')
    clonal_rearrangement = outtable.replace('.xlsx', '_clonal-rearrangements.csv')

    # define fields for Vgene and Vgroup from the column Vregion
    hom.loc[:, 'Vregion'] = hom['Vregion'].str.replace('D-','-D')
    hom.loc[:, 'Vgene'] = hom['Vregion'].str.split('-').str[:2].str.join('-')

    ## changed -1 for 2, to include only the first two elements of the notation in V alleles
    #hom.loc[:, 'Vgroup'] = hom['Vregion'].str.split('-').str[0]

    ## replace whole CDR3 if it starts by certain patterns that indicate false productivity
    if cdr3 == 'major':
        CDR3 = 'major_CDR3'
    else:
        CDR3 = 'CDR3'
    
    hom[CDR3] = hom[CDR3].replace(regex=[r'\b(^CITV\w*)\b',
                                               r'\b(^CE\w*)\b', r'\b(^C(V*)LL\w*)\b'], value='')
        
    # get only needed columns
    homf = hom[['sample_name', 'Vregion', 'reads_mapped', 'Vgene', 'Vgroup', CDR3]]
    homf.columns = ['sample_name', 'Vregion', 'reads_mapped', 'Vgene', 'Vgroup', 'CDR3']
    
    ## get info per gene
    hom_new = homf.merge(homf.groupby(['sample_name','Vgene'],
                                      as_index=False).sum(), on=['sample_name','Vgene'])
    
    hom_new = hom_new.merge(hom_new.groupby(['sample_name'],
                                            as_index=False).sum(), on=['sample_name'])

    del(hom_new['reads_mapped_y_y'])
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    
    ## percentage fields
    hom_new['percent_reads_mapped_Vallele'] = hom_new['reads_mapped_x_x']/hom_new['reads_mapped_x_y'] * 100
    hom_new['percent_reads_mapped_Vgene'] = hom_new['reads_mapped_y_x']/hom_new['reads_mapped_x_y'] * 100

    ## filter out rescued rearrangements from FR3 filter if %reads_mapped_Vgene < 90%
    if not rescued.empty:
        rescued.loc[:, 'Vregion'] = rescued['Vregion'].str.replace('D-','-D')
        for e in rescued['sample_name'].unique().tolist():
            mrescued = rescued[rescued['sample_name'] == e]
            for g in mrescued['Vregion'].unique().tolist():
                r = hom_new[(hom_new['sample_name'] == e) & (hom_new['Vregion'] == g)]
                if float(r['percent_reads_mapped_Vgene'].to_string().split()[1]) < 90:
                    hom_new = hom_new[~((hom_new['sample_name'] == e) & (hom_new['Vregion'] == g))]
            
    
    ## merge grouped table with original
    hom_plots = hom_new.merge(hom[['sample_name', 'Vregion', 'J_assigned']], on=['sample_name', 'Vregion'])
    
    # select allele with productive CDR3 and more reads
    ## Iterate per sample
    samples = hom_new['sample_name'].unique().tolist()
    
    for s in samples:

        ## calculate percentages
        n1 = hom_new[hom_new['sample_name'] == s]
        
	# For each Vgene, we will evaluate CDR3 on each allele
        list_Vgene = n1['Vgene'].unique().tolist()
        for l in list_Vgene:
	 
            nm = n1[n1['Vgene'] == l]
            productive_data = nm[nm['CDR3'].str.endswith('W', na=False) & nm['CDR3'].str.startswith('C', na=False)]
                                
            if productive_data['Vregion'].tolist() != []:
                unproductive_data = nm[~nm['Vregion'].isin(productive_data['Vregion'].tolist())]
                alleles_to_remove = unproductive_data['Vregion'].tolist()
            else:
                # if there is no productive data I am not removing any alleles
                alleles_to_remove = []
	 
                for a in alleles_to_remove:
                    hom_new = hom_new[~((hom_new['sample_name'] == s) & (hom_new['Vregion'] == a))]
	
	
    #retrieve the allele with max number of mapped reads within Vgene group
    hom_max = hom_new.loc[hom_new.groupby(['sample_name', 'Vgene'])['reads_mapped_x_x'].idxmax()]

    # merge this with the original dataframe to retrieve the max allele fo each V gene and the total counts per V gene
    hom_to_save1 = hom_max[['sample_name', 'Vregion','reads_mapped_x_x', 'Vgene',
	                        'reads_mapped_y_x', 'percent_reads_mapped_Vallele', 'percent_reads_mapped_Vgene']]
    
    hom_to_save = hom.merge(hom_to_save1, on=['sample_name','Vregion'])
    hom_to_save.loc[:, 'nrearrangement'] = hom_to_save.groupby('sample_name')['reads_mapped_y_x'].rank(ascending=False, method='first')
    ## delete undesired columns
    del hom_to_save['reads_mapped_x_x']
    del hom_to_save['Vgene_y']
    del hom_to_save['Vgroup']
    
    hom_to_save['percent_clonal_Vgene'] = 0
    ## calculate percentages
    samples = hom_to_save['sample_name'].unique().tolist()
    
    for s in samples:
           
        n1 = hom_to_save[hom_to_save['sample_name'] == s]
        rear = sorted(n1['nrearrangement'].tolist())
           
        for e in range(0, len(rear)):
    
            t = n1[n1['nrearrangement'] == rear[e]].index[0]               
            f = e
            mr = float(n1[n1['nrearrangement'] == rear[e] ]['reads_mapped_y_x'].to_string().split()[1])
            
            if e == (len(rear) -1):
                ms = 0
            else:
                ms = float(n1[n1['nrearrangement'] == rear[e+1]]['reads_mapped_y_x'].to_string().split()[1])
            
                while ((mr/2) < ms):
                    if (f+2 <= (len(rear) - 1 )):
                        
                        f += 1
                        ms = float(n1[n1['nrearrangement'] == rear[f+1]]['reads_mapped_y_x'].to_string().split()[1])
                    else:
                        ms = 0
        
            hom_to_save.at[t, 'percent_clonal_Vgene'] = (mr - ms)/float(n1['count_reads_mapped'].unique().tolist()[0]) * 100
            hom_to_save['percent_clonal_Vgene'] = hom_to_save['percent_clonal_Vgene'].astype(float)
            
    ## define columns
    hom_to_savef = hom_to_save[['sample_name', 'Vregion', 'combined_alleles', 'reads_mapped', 'reads_mapped_leader',
	                        'reads_mapped_fr1', 'reads_mapped_fr2', 'reads_mapped_fr3',
	                        'region_length', 'nreads', 'count_reads_mapped', 'nrearrangement',
	                        'Vgene_x', 'reads_mapped_y_x', 'percent_reads_mapped_Vgene',
	                        'percent_clonal_Vgene','homology-IGHV_noFR3',
	                        'mutational_status_noFR3', 'J_assigned', 'IGHV-J',
	                        'consensus_length', CDR3, 'CDR3_length',  'IGHD', 'IGHD_emboss',
	                        'insertions','deletions','ORF disruption', 'majorproductive_seq',
                                'nreads_majorproductive_seq']]

    hom_to_savef.columns = ['sample_name', 'Vallele', 'combined_alleles', 'reads_mapped_allele', 'reads_mapped_leader',
	                    'reads_mapped_fr1', 'reads_mapped_fr2', 'reads_mapped_fr3',
	                    'ref_length', 'reads_trimming', 'reads_mappedV', 'nrearrangement', 'Vgene',
	                    'reads_mapped_gene', 'percent_reads_mapped_Vgene', 'percent_clonal_Vgene',
                            'homology-IGHV_noFR3', 'mutational_status_noFR3',
	                    'J_assigned', 'IGHV-J', 'consensus_length', 'CDR3',
	                    'CDR3_length', 'IGHD','IGHD_emboss', 'insertions','deletions',
	                    'ORF disruption','majorunique_productiveSeq','nreads_majorunique_productiveSeq']
	
    ## select the max entry for each V gene
    ## save tables
    hom_to_savef = hom_to_savef.copy()
    ## chord plot
    cdr3chord(hom_to_savef)

    ## predict status
    hom_to_savef['rearrangement'] = hom_to_savef['sample_name'] + '_' + hom_to_savef['nrearrangement'].astype(str)
    
    hom_to_savef = predicted_status(hom_to_savef)
    clone_plotting(hom_to_savef)
    
    ## colour max value and save table
    l = list(hom_to_savef.groupby('sample_name', as_index=False)['reads_mapped_gene'].idxmax())
    hom_to_savef.style.apply(lambda x: ['background-color: #00ff8c' if x.name in l else '' for i in x],
                             axis=1).to_excel(outtable, engine='openpyxl', index=False)
    
    hom_to_savef.to_csv(outtable.replace('.xlsx', '.csv'), sep=',')
    
    ## save first rearrangement
    hom_to_savef.iloc[l]
    hom_to_savef.iloc[l].to_csv(first_rearrangement, sep=',')
	
    ## retrieve top 3 clones
    l2 = []
    list_samples = hom_to_savef['sample_name'].unique().tolist()
	
    for i in list_samples:
        n = hom_to_savef[hom_to_savef['sample_name'] == i].nlargest(3,
                                                                    'reads_mapped_gene')
        l2 += n.index.tolist()
    
    ## save the first three rearrangements for each sample
    hom_to_savef.iloc[l2]
    hom_to_savef.iloc[l2].to_csv(top_rearrangement, sep=',')

    ## calculate coverage of the clonal rearrangements and include coverage info in homology_resume
    hom_to_savef[hom_to_savef['clone_status'] == 'CLONAL'].to_csv(clonal_rearrangement, sep=',', index=False)
    coverage_IGHs.calculate_coverage(outcoverage, out, clonal_rearrangement, 10, False, min_cov)
    ## annotate coverage information
    covt = pd.read_csv(cov_table, sep=',')
    colcov = '%{}X'.format(min_cov)
    covt['col'] = covt['sample_name'].str.replace('-sorted', '')
    covtab = covt[['col', colcov]]
    hom_to_savef['col'] = hom_to_savef['sample_name'] + '_' + hom_to_savef['IGHV-J']
    hom_to_saveff = hom_to_savef.merge(covtab, how='left', on='col')

    ## make accessible consensus sequences from clonal rearrangements
    list_clonal = hom_to_saveff[hom_to_savef['clone_status'] == 'CLONAL']['col'].unique().tolist()
    
    if not os.path.isdir(folder_consensus):
        os.mkdir(folder_consensus)
    for ffasta in list_clonal:
        fasta = glob.glob(os.path.join(out, 'results/consensus_complete/{}*.fa'.format(ffasta)))
        [shutil.copy(f, folder_consensus) for f in fasta]
        
        
    hom_to_saveff = hom_to_saveff.drop(columns=['col'])

    ## colour max value and save table                                                                                                                                                                     
    l = list(hom_to_saveff.groupby('sample_name', as_index=False)['reads_mapped_gene'].idxmax())
    hom_to_saveff.style.apply(lambda x: ['background-color: #00ff8c' if x.name in l else '' for i in x],
                             axis=1).to_excel(outtable, engine='openpyxl', index=False)

    hom_to_saveff.to_csv(outtable.replace('.xlsx', '.csv'), sep=',')
    ## save first rearrangement                                                                                                                                                                          
    hom_to_saveff.iloc[l]
    hom_to_saveff.iloc[l].to_csv(first_rearrangement, sep=',')

    hom_to_saveff.style.apply(lambda x: ['background-color: #00ff8c' if x.name in l else '' for i in x],
                             axis=1)

    return hom_to_saveff


def filterFR3(filtered, rest, CDR3='major_CDR3'):

    """"""
    rescued = pd.DataFrame()
    for e in rest['sample_name'].unique().tolist():
    
        s = rest[rest['sample_name'] == e ]
        s.loc[:, 'Vgene'] = s['Vregion'].str.split('-').str[:2].str.join('-')
        df = filtered[filtered['sample_name'] == e]
        for f in s['Vregion'].unique().tolist():

            r = s[s['Vregion'] == f ]
            family = r['Vgroup'].to_string().split()[1]
            cdr3 = r[CDR3].to_string().split()[1]
            JH = r['J_assigned'].to_string().split()[1]
            total_reads_sample = sum(df['reads_mapped']) + sum(s['reads_mapped'])
            
            if int(r['reads_mapped'].to_string().split()[1]) >= (total_reads_sample * 0.8):
                rescued = pd.concat([rescued, r]).sort_values(['sample_name'])
                filtered = pd.concat([filtered, r]).sort_values(['sample_name'])
            
            else:
                
                if cdr3 != 'None':
                    subdf = df[(df['Vgroup'] == family) & (df[CDR3] == cdr3) & (df['J_assigned'] == JH)]
                    #print(subdf.to_string())
                    if not subdf.empty:
                        ## if passing filter alleles are more than one
                        ## with same VH and CDR3, pick the most dominant
                        if len(subdf) > 1:
                    
                            max_allele = subdf.loc[subdf['reads_mapped'].idxmax()]['Vregion']
                            subdf_max = subdf[subdf['Vregion'] == max_allele]
                                                
                        else:
                            subdf_max = subdf
                            max_allele = subdf['Vregion'].to_string().split()[1]
                             
                        i = filtered[(filtered['sample_name'] == e) & (filtered['Vregion'] == max_allele)].index[0]
                        
                        filtered.at[i, 'reads_mapped'] = filtered.at[i, 'reads_mapped'] + r['reads_mapped']
                        filtered.at[i, 'reads_mapped_leader'] = filtered.at[i, 'reads_mapped_leader'] + r['reads_mapped_leader']
                        filtered.at[i, 'reads_mapped_fr1'] = filtered.at[i, 'reads_mapped_fr1'] + r['reads_mapped_fr1']
                        filtered.at[i, 'reads_mapped_fr2'] = filtered.at[i, 'reads_mapped_fr2'] + r['reads_mapped_fr2']
                        filtered.at[i, 'reads_mapped_fr3'] = filtered.at[i, 'reads_mapped_fr3'] + r['reads_mapped_fr3']

    #print(filtered[filtered['sample_name'] == e].to_string())
    return filtered, rescued
                

def usage_plots(hom_to_savef2, homt, folder_plots, naleles, tag=''):

    """"""

    hom_to_savef2 = hom_to_savef2.sort_values(['Vgene', 'Vallele']).reset_index(drop=True)

    # plot general allele and IGHV gene counts
    # size a4 paper
    plt.figure(figsize=(naleles/1.3, naleles/3))
    sns.set(font_scale=naleles/40, style='white')

    ## allele
    az = sns.barplot(x='Vregion', y='reads_mapped',
                     data=homt.sort_values('Vregion'))
    az.set_xticklabels(az.get_xticklabels(), rotation=90)

    plt.tight_layout()
    name = os.path.join(folder_plots, 'counts-Valleles_{}{}.png'.format(tag, str(run)))
    plt.savefig(name)
    plt.gcf().clear()
        
    ## percent

    # size a4 paper
    plt.figure(figsize=(50,30))
    sns.set(font_scale=3, style='white')
    allelesJ = sorted(list(hom_to_savef2['J_assigned'].dropna().unique()))
    ## gene usage run
    pivot_df = hom_to_savef2.pivot_table(index='Vgene', columns='J_assigned',
                                        values='reads_mapped_gene',
                                        aggfunc=sum).fillna(0)
    plotname = os.path.join(folder_plots, 'counts-genes-IGHVJ_{}{}.png'.format(tag, str(run)))
    # counts
    ## size a4 paper
    plt.figure(figsize=(15,8))
    with sns.color_palette("Paired", 15):
        az = pivot_df.plot.bar(stacked=True, figsize=(30,15))
    az.set_xticklabels(az.get_xticklabels(), rotation=90)

    plt.tight_layout()
    plt.savefig(plotname)
    plt.gcf().clear()

    # percent
    pivot_df = hom_to_savef2.pivot_table(index='Vgene', columns='J_assigned',
                            values='percent_reads_mapped_Vgene',
                                        aggfunc=sum).fillna(0)
    plotname = os.path.join(folder_plots, 'percent-genes-IGHVJ_{}{}.png'.format(tag, str(run)))
    ## size a4 paper
    plt.figure(figsize=(15,8))
    with sns.color_palette("Paired", 15):
        az = pivot_df.plot.bar(stacked=True, figsize=(30,15))
    az.set_xticklabels(az.get_xticklabels(), rotation=90)

    plt.tight_layout()
    plt.savefig(plotname)
    plt.gcf().clear()


def main():

    ## HOMOLOGY RESUMEN TABLE ##
    ### total reads annotation
    re = pd.read_excel(reads)
    
    re['sample_name'] = re['Muestra'].str.split('_').str[0]
    re['nreads'] = re['numero lecturas trim']*2
    re1 = re[['nreads', 'sample_name']].drop_duplicates()
    
    ### reads mapped against V annotation
    rm = pd.read_csv(resume_Vmapping, sep=',')
    rm.columns = ['sample_name', 'count_reads_mapped',
                  'percent_reads_mapped']
    re2 = re1.merge(rm)
    
    outtable = os.path.join(out, 'homology_resume_{}.xlsx'.format(run))
    outtable2 = os.path.join(out, 'homology_resume_filterFR3_{}.xlsx'.format(run))
    outcoverage = os.path.join(out, 'coverage')

    cov_table = os.path.join(outcoverage, 'samplecovstats.csv')
    filter_FR3 = os.path.join(out, 'homology_table_filterFR3.csv')
    folder_plots = os.path.join(out, 'usage_plots')
    folder_consensus = os.path.join(out, 'consensus_sequences_clonal')    
    # read homology_table.csv
    homt = pd.read_csv(input_table, sep=',')

    homt.loc[:, 'Vgroup'] = homt['Vregion'].str.split('-').str[0]
    homt2 = homt.merge(re2)
    
    ## filter rearrangements supported by one fragment only
   
    hom_filt = homt2[((homt2['reads_mapped'] - homt2['reads_mapped_fr1']) > 2)
                    & ((homt2['reads_mapped'] - homt2['reads_mapped_fr2']) > 2)
                    & ((homt2['reads_mapped'] - homt2['reads_mapped_fr3']) > 2)]

    hom_drop = homt2[((homt2['reads_mapped'] - homt2['reads_mapped_fr1']) <= 2)
                    | ((homt2['reads_mapped'] - homt2['reads_mapped_fr2']) <= 2)
                    | ((homt2['reads_mapped'] - homt2['reads_mapped_fr3']) <= 2)]

    ## filter rearrangement if 90% of reads are FR3 supported
    ## add read % and add or percent_reads_rearrangement > 92
    hom_filterFR3 = hom_filt[(hom_filt['reads_mapped_fr3']/hom_filt['reads_mapped'] < 0.92)] ## or percent reads >90
    hom_drop2 = hom_filt[hom_filt['reads_mapped_fr3']/hom_filt['reads_mapped'] >= 0.92] ## and percent reads < 90
    result = pd.concat([hom_drop, hom_drop2], ignore_index=True)
    
    ## iterate and for each of them, sum the reads totally and by fragments
    ### reads will be added to the major rearrangement of the same VH family if CDR3 is equal    
    hom_filterFR3defmajor, rescued = filterFR3(hom_filterFR3.reset_index().copy(), result, CDR3='major_CDR3')
    hom_filterFR3.to_csv(filter_FR3, sep=',')
             
    hom_to_savef2 = homology_resume(hom_filterFR3defmajor.reset_index(), rescued,
                                    outtable2.replace('.xlsx', '-majorcdr3.xlsx'),
                                    outcoverage, mincov, cov_table, folder_consensus,
                                    cdr3='major')
               
    
    ## DATABASE UPDATE ##
    
    alleles = hom_to_savef2['Vallele'].nunique()
    genes = hom_to_savef2['Vgene'].nunique()
    naleles = len(homt['Vregion'].unique().tolist())
    print('n unique IGHV alleles: {}'.format(alleles))
    print('n unique IGHV genes: {}'.format(genes))
    
    ## plots
    if not os.path.isdir(folder_plots):
        os.mkdir(folder_plots)
    
    usage_plots(hom_to_savef2, homt, folder_plots, naleles)
    samples = hom_to_savef2['sample_name'].unique().tolist()
    for s in samples:
        ds = hom_to_savef2[hom_to_savef2['sample_name'] == s]
        usage_plots(ds, homt, folder_plots, naleles, tag='{}-'.format(s))
        

if __name__ == "__main__":
    main()
