import streamlit as st
import pandas as pd
import numpy as np
import base64
from scipy import stats
import csv
from ete3 import NCBITaxa
import datetime
import altair as alt
import plotly.express as px
import re
from datetime import date
import seaborn as sns
import webbrowser
from streamlit.script_runner import RerunException
import matplotlib.pyplot as plt


## load data from NCBI ftp site
@st.cache(allow_output_mutation=True) 
def load_data(DATA_URL):
    '''
    function to load prokaryote genome metadata from the NCBI ftp site
    :param DATA_URL: NCBI ftp site
    :return: genome_df: a dataframe of the genome metadata
    '''

    data = pd.read_table(DATA_URL, usecols = [0,1,6,7,8,11,12,13,15,20])  # the column 15: Status
    complete_genomes = data.loc[data['Status'] == "Complete Genome"]    # Select only complete genomes
    genome_df = complete_genomes.reset_index(drop = True) # reset index from 0
    genome_df = genome_df.drop(['Status'], axis = 1) # remove the column "Status"
    genome_df = genome_df.rename(columns={"#Organism/Name": "Genome Name"}) # change the first column name
    return genome_df

def rerun():
    raise st.script_runner.RerunException(st.script_request_queue.RerunData(None))


def change_ftp(ftp_path):
    '''
    function to change from ftp: to https: of the downloadable table
    :param ftp_path: a string of the ftp path
    :return:
    '''
    new_ftp_path = 'https' + ftp_path[3:]  
    return new_ftp_path

# def path_hyperlink(ftp_path):
#     webbrowser.open_new_tab(ftp_path)

# def make_clickable(val):
#     # target _blank to open new window
#     return '<a target="_blank" href="{}">{}</a>'.format(val, val)


@st.cache(allow_output_mutation=True) 
def get_desired_ranks(taxid, desired_ranks):
    '''
    function to get taxonomical information from taxID
    :param taxid: e.g., 257
    :param desired_ranks: taxonimical ranks
    :return:
    '''
    ncbi = NCBITaxa()
    count = 0   
    lineage_list = [taxid]
    try:
        lineage = ncbi.get_lineage(taxid)
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        rank_dict = {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}
                
        for k, v in rank_dict.items():
            if v != "<not present>":
                try:
                    new_v = ncbi.get_taxid_translator([v])  ## new_v is dictionary key=pylum_id, value=xxxxxx
                    lineage_list.append(new_v[v])
                except:
                    pass
            else: 
                lineage_list.append("<not present>")
            
    except:
        count = 1
    return (lineage_list, count)

@st.cache(allow_output_mutation=True)
def taxID_lineage_df(taxids, desired_ranks):
    '''
    Function to make a taxonomical information dtatframe of taxids
    :param taxids:
    :param desired_ranks: superkingdom, phylum, class, order, family, genus, species
    :return: lineage_df: a dataframe of taxids and their taxonomical information
    '''

    missed_taxID = 0
    lineage_dict = {}   ## key=taxid, value=[lineage]
    for i, taxid in enumerate(taxids):
        (lineage_list, count) = get_desired_ranks(taxid, desired_ranks)
        lineage_dict[i] = lineage_list
        if count == 1:
            missed_taxID +=1
    
    desired_ranks.insert(0, 'TaxID')
    lineage_df = pd.DataFrame.from_dict(lineage_dict, orient='index', columns=desired_ranks)
    return (lineage_df, missed_taxID)

## download data
@st.cache(allow_output_mutation=True)
def download_link(object_to_download, download_filename, download_link_text):
    """
    Generates a link to download the given object_to_download.

    object_to_download (str, pd.DataFrame):  The object to be downloaded.
    download_filename (str): filename and extension of file. e.g. mydata.csv, some_txt_output.txt
    download_link_text (str): Text to display for download link.

    Examples:
    download_link(YOUR_DF, 'YOUR_DF.csv', 'Click here to download data!')
    download_link(YOUR_STRING, 'YOUR_STRING.txt', 'Click here to download your text!')

    """
    if isinstance(object_to_download,pd.DataFrame):
        object_to_download = object_to_download.to_csv(index=False)

    # some strings <-> bytes conversions necessary here
    b64 = base64.b64encode(object_to_download.encode()).decode()

    return f'<a href="data:file/txt;base64,{b64}" download="{download_filename}">{download_link_text}</a>'


@st.cache(allow_output_mutation=True)
def making_all_select_genome_df(final_genome_df, items_selected_text, subitems_selected):
    '''
    Function to select a certain group of genome metadata
    :param final_genome_df:
    :param items_selected_text:
    :param subitems_selected:
    :return: 3 information
    '''


    # print(type(subitems_selected), subitems_selected)
    col_names = final_genome_df.columns.values.tolist()
    final_all_selected_genome = pd.DataFrame(columns=col_names)
    
    # print(len(subitems_selected))
    # print('items_selected_text: ', items_selected_text)
    selected_genomes = []
    
    if 1 <= len(subitems_selected) <= 10:
        for t in subitems_selected:
            selected_by_subitem = final_genome_df.loc[final_genome_df[items_selected_text] == t]
            selected_genomes.append(selected_by_subitem)
        
        all_selected_genomes = pd.concat(selected_genomes)  ## Concatenating dataframes
        final_all_selected_genomes = pd.concat((final_all_selected_genome, all_selected_genomes))
    
        final_all_selected_genome_shape = final_all_selected_genomes.shape
        if final_all_selected_genome_shape[0] >= 1:
            final_all_selected_genome = final_all_selected_genome.reset_index(drop = True)
        
            # final_genome_df['GC%'] = pd.to_numeric(final_all_selected_genome['GC%'], errors='coerce') 
            # final_genome_df['Genes'] = pd.to_numeric(final_all_selected_genome['Genes'], errors='coerce') 
            # final_genome_df['Proteins'] = pd.to_numeric(final_all_selected_genome['Proteins'], errors='coerce')
            final_genome_df['Chromosome'] = pd.to_numeric(final_all_selected_genome['Chromosome'], errors='coerce')
            final_genome_df['Plasmid'] = pd.to_numeric(final_all_selected_genome['Plasmid'], errors='coerce')
        
        selected_genomes = final_genome_df.shape
        return (final_all_selected_genomes, selected_genomes[0], final_all_selected_genome_shape[0])

    elif len(subitems_selected) > 10:
        message = 'Please select no more than 10 items'
        return message
    
    elif len(subitems_selected) == 0:
        message = 'Please select item(s)'
        return message



    

@st.cache(allow_output_mutation=True)    
def barchart_maker(dataframe, column_name):
    '''
    Function to make a input dataframe for barchart plotting of genomic data
    argument 1: dataframe
    argument 2: column_name
    '''
    if column_name == 'Release Date':
        ## No. of genome per year
        dataframe['Release Date'] = pd.to_datetime(dataframe['Release Date'])
        dataframe['year'] = dataframe['Release Date'].dt.year   # extract year from date datatype
        year_series = dataframe['year'].value_counts()    # count year values and save it to a Series dattype
        year_df = pd.DataFrame({'Year': year_series.index, 'Count': year_series.values})
        year_index_df = year_df.set_index('Year')
        # st.bar_chart(year_index_df)
        return year_index_df
        
    else:
        column_name_series = dataframe[column_name].value_counts()    # count year values and save it to a Series dattype
        column_name_df = pd.DataFrame({column_name: column_name_series.index, 'Count': column_name_series.values})
        column_index_df = column_name_df.set_index(column_name)
        return column_index_df

@st.cache(allow_output_mutation=True) 
def boxplot_maker(graph_choice, genome_dataframe):
    
    '''
    Function to draw a boxplot of the genomic features selected
    argument 1: graph_choice
    argument 2: genome dataframe
    '''
    fig = px.box(genome_dataframe, y=graph_choice, points='all')
    return fig

@st.cache(allow_output_mutation=True) 
def scatterplot_maker(scatter_items_selected, genome_dataframe):
    '''
    Function to draw a boxplot of the genomic features selected
    argument 1: scatter_items_selected (datatype: list)
    argument 2: genome dataframe
    '''
    
    if len(scatter_items_selected) == 2:
        fig = px.scatter(genome_dataframe, x=scatter_items_selected[0], y=scatter_items_selected[1], trendline='ols', trendline_color_override="red", color='superkingdom')
        # fig = px.scatter(genome_dataframe, x=scatter_items_selected[0], y=scatter_items_selected[1])
        # fig = sns.jointplot(data=genome_dataframe, x = scatter_items_selected[0], y = scatter_items_selected[1], kind='reg')
        return fig

    elif len(scatter_items_selected) == 3:
        fig = px.scatter_3d(genome_dataframe, x=scatter_items_selected[0], y=scatter_items_selected[1], z=scatter_items_selected[2], opacity=0.5)
        return fig
    
@st.cache(allow_output_mutation=True) 
def scatterplot_maker_3(scatter_items_selected, genome_dataframe):
    '''
    Function to draw a boxplot of the genomic features selected
    argument 1: scatter_items_selected (datatype: list)
    argument 2: genome dataframe
    '''
    
    if len(scatter_items_selected) == 2:
        # st.subheader(f'2d scatter plot of {scatter_items_selected[0]} and {scatter_items_selected[1]}')
        # fig = px.scatter(genome_dataframe, x=scatter_items_selected[0], y=scatter_items_selected[1], trendline='ols', trendline_color_override="red", color='superkingdom')
        # fig = px.scatter(genome_dataframe, x=scatter_items_selected[0], y=scatter_items_selected[1])

        classes = genome_dataframe['superkingdom']
        x = scatter_items_selected[0]
        y = scatter_items_selected[1]
        ## calculation of R2 value
        slope, intercept, r_value, p_value, std_err = stats.linregress(genome_dataframe[x], genome_dataframe[y])
        R2 = r_value ** 2
        title = "y={0:.1f}x+{1:.1f} (R^2 value: {2:.2})".format(slope, intercept, R2)

        # fig = sns.jointplot(data=genome_dataframe, x = scatter_items_selected[0], y = scatter_items_selected[1], kind='reg', ratio=7, height=7, scatter=False)
        fig = sns.jointplot(data=genome_dataframe, x = scatter_items_selected[0], y = scatter_items_selected[1], kind='reg')
        fig.ax_joint.scatter(x, y, c=classes)
        fig.fig.suptitle(title)
        fig.ax_joint.collections[0].set_alpha(1)
        fig.fig.tight_layout()
        fig.fig.subplots_adjust(top=0.95)
        
        return fig

    elif len(scatter_items_selected) == 3:
        # st.subheader(f'3d scatter plot of {scatter_items_selected[0]}, {scatter_items_selected[1]}, and {scatter_items_selected[2]}')
        fig = px.scatter_3d(genome_dataframe, x=scatter_items_selected[0], y=scatter_items_selected[1], z=scatter_items_selected[2], opacity=0.5)
        # st.plotly_chart(fig)
        return fig

@st.cache(allow_output_mutation=True)
def genome_selector_features(final_genome_df, selected_features):
    '''
    Function to select a certain group of genome data based on the selected genomeic features
    :param final_genome_df:
    :param selected_features:
    :return: all_selected_genome, all_genomes[0], selected_genomes[0]
    '''


    if len(selected_features) == 0:
        message = 'Select genomical feature(s)'
        return message

    else:
        count = 1
        for k, v in selected_features.items():
            if count == 1:
                selected_df = final_genome_df[final_genome_df[k].between(v[0], v[1])]
                count += 1
            else:
                selected_df = selected_df[selected_df[k].between(v[0], v[1])]
                count += 1
        
        all_selected_genome = selected_df.reset_index(drop = True)
        # st.write(all_selected_genome)
        all_genomes = final_genome_df.shape
        selected_genomes = all_selected_genome.shape
        
    return all_selected_genome, all_genomes[0], selected_genomes[0]
    

def main():
    
    st.image('./picture7.png', use_column_width=True)
    
    pd.set_option('mode.chained_assignment', None)
    ## load genome data from NCBI 
    DATA_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
    genome_df = load_data(DATA_URL)
            
    ## Count the numbers of chromosome and plasmid of each genomes and 
    @st.cache(allow_output_mutation=True)
    def count_chro_plas(genome_df):
        '''
        Function to update the final genome_df by adding two columns, i.e., No. of chromosome and No. of plasmid
        :param genome_df:
        :return: the updated genome_df
        '''

        replicons_list = list(genome_df['Replicons'])
        chromosome_list = []
        plasmid_list = []
        word_chromosome = ['chromosome']
        word_plasmid = ['plasmid']
        exactMatch_chromosome = re.compile(r'\b%s\b' % '\\b|\\b'.join(word_chromosome), flags=re.IGNORECASE)
        exactMatch_plasmid = re.compile(r'\b%s\b' % '\\b|\\b'.join(word_plasmid), flags=re.IGNORECASE)
        for r in replicons_list:
            no_chromosome = len(exactMatch_chromosome.findall(r))
            chromosome_list.append(no_chromosome)
            no_plasmid = len(exactMatch_plasmid.findall(r))
            plasmid_list.append(no_plasmid)
        
        chr_plas_dic = {'Chromosome': chromosome_list, 'Plasmid': plasmid_list}
        chr_plas_df = pd.DataFrame(chr_plas_dic)
        genome_df = pd.concat([genome_df, chr_plas_df], axis=1)
        
        return genome_df
    
    
    genome_df = count_chro_plas(genome_df)

    ## Reorder of columns of the genome_df
    @st.cache(allow_output_mutation=True)
    def making_final_df (genome_df):
        '''
        Function to make a final genome_df from the genome_df by reordering of columns, changing datatypes, etc...
        :param genome_df:
        :return: final_genome_df
        '''

        columnsTitles = ['Genome Name', 'TaxID', 'Size (Mb)', 'GC%', 'Replicons', 'Chromosome', 'Plasmid', 'Genes', 'Proteins', 'Release Date', 'FTP Path']
        genome_df = genome_df.reindex(columns=columnsTitles)
        genome_df=genome_df.rename(columns = {'Size (Mb)':'Genome size (Mb)', 'FTP Path':'Genome download (FTP Path)'})
        
        taxID_list = list(set(genome_df['TaxID']))
        ncbi = NCBITaxa()
        desired_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        lineage_df, missed_taxID = taxID_lineage_df(taxID_list, desired_ranks)
        
        ## Join of two dataframes
        final_genome_df = pd.merge(genome_df, lineage_df, on='TaxID')
        
        ## ftp to https
        final_genome_df['Genome download (FTP Path)'] = final_genome_df['Genome download (FTP Path)'].apply(change_ftp)  
        
        # OK_df = final_genome_df.style.format({'Genome download (FTP Path)': make_clickable})
        
        ## str to numeric data of the columns or number to str
        final_genome_df['GC%'] = pd.to_numeric(final_genome_df['GC%'], errors='coerce') 
        final_genome_df['Genes'] = pd.to_numeric(final_genome_df['Genes'], errors='coerce') 
        final_genome_df['Proteins'] = pd.to_numeric(final_genome_df['Proteins'], errors='coerce')
        final_genome_df['Chromosome'] = pd.to_numeric(final_genome_df['Chromosome'], errors='coerce')
        final_genome_df['Plasmid'] = pd.to_numeric(final_genome_df['Plasmid'], errors='coerce')
        final_genome_df['TaxID'] = final_genome_df['TaxID'].astype(str)

        return final_genome_df
    
    final_genome_df = making_final_df (genome_df)
    ## main page
    # st.title('Complete Prokaryote Genome Metadata Explorer')
    # c1, c2 = st.beta_columns(2)
    # with c2:
    #     st.write('**Created by Microbiology Team**') 
    
    ## Introduction; reference; data sources; etc
    with st.beta_expander("See notes"):

        st.markdown("""
        This app allows the user to easily access and explore the metadata of complete prokaryote genomes to support education and genome-based research.  
        
        * **All complete prokaryte genomes**
        * **Genome selection by taxonomic lineage**
        * **Genome selection by numerical genomic features, including genome sieze, GC%, etc.**
        
        """)
        
        
    ## Data 1: Date and the no. of genomes
    today = date.today()
    d1 = today.strftime("%B %d, %Y")
    
    ## No. of genomes
    (no_genomes, no_columns) = final_genome_df.shape
    st.write(f'**Data sourece:** GenBank prokaryotes.txt file downloaded **_{d1}_** **(a total of {no_genomes} genomes) (https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/)**')
    
    if st.checkbox('Show raw prokaryote genome table (5 genomes only) and download all genome table'):
        st.dataframe(final_genome_df.head())
                
        if st.button('Download complete genome data as CSV', key = 1):
            tmp_download_link = download_link(final_genome_df, 'Genome_infor.csv', 'Click here to download your data!')
            st.markdown(tmp_download_link, unsafe_allow_html=True)
    
    # st.write('Description of the raw data')
    # st.write(final_genome_df.isnull().sum())    # missing values 
    final_genome_df = final_genome_df.dropna()  # remove missing values
    # st.write(final_genome_df.describe(include=['O']))   # description of object columns
    # st.write(final_genome_df.describe())    # description of numerical data columns
    # st.write(final_genome_df.isnull().sum())

    
    ## Sidebar
    st.sidebar.image('./Picture10.png')
    # st.sidebar.image('./bacteria_all.png', use_column_width=True)
    # st.sidebar.image('./bacteria_256_test.png')
    # st.sidebar.success('**Complete Prokaryote Genome Metadata Explorer**')    
    menu = ['1. All prokaryote genomes', '2. Selected by taxonomic lineage', '3. Selected by numerical feature(s)']
    choice = st.sidebar.selectbox('Select mode', menu, key = 1)
    
    if choice == '1. All prokaryote genomes':
        st.write('')
        st.write('')
        st.success('**All complete prokaryote genomes**')
        
        ## Data 1:
        st.markdown('**1. Number of bacterial and archaeal genomes sequenced each year and submitted to NCBI.**')
        year_index_df = barchart_maker(final_genome_df, 'Release Date')
        
        st.bar_chart(year_index_df)
                
        ## Data 2: Descriptive statistics of genomes
        st.markdown('**2. Descriptive statistics of complete genomes (boxplot)**')
        graph_menu = ['Genome size (Mb)', 'Chromosome', 'Plasmid', 'GC%', 'Genes', 'Proteins']
        graph_choice = st.selectbox('Boxplot data of genome size, No. of chromosome, plasmid, gene, and protein, and GC%', graph_menu, key = 1, index = 0)
        
        if graph_choice == 'Select one':
            st.write('Please select a boxplot')
        
        elif graph_choice == 'GC%':
            fig = boxplot_maker(graph_choice, final_genome_df)
            st.plotly_chart(fig)
        
        elif graph_choice == 'Genome size (Mb)':
            fig = boxplot_maker(graph_choice, final_genome_df)
            st.plotly_chart(fig)
            
        elif graph_choice == 'Chromosome':
            fig = boxplot_maker(graph_choice, final_genome_df)
            st.plotly_chart(fig)
            
        elif graph_choice == 'Plasmid':
            fig = boxplot_maker(graph_choice, final_genome_df)
            st.plotly_chart(fig)
                
        elif graph_choice == 'Genes':
            fig = boxplot_maker(graph_choice, final_genome_df)
            st.plotly_chart(fig)
            
        elif graph_choice == 'Proteins':
            fig = boxplot_maker(graph_choice, final_genome_df)
            st.plotly_chart(fig)
        
        ## Data 3: Descriptive statistics of genomes (scatter plot)
        st.markdown('**3. Descriptive statistics of complete genomes (Pearson correlation heatmap and scatter plot [2d or 3d])**')
        # st.write('<For a 3d scatter plot, please select three values>')
        st.write('Pearson correlation heatmap')
        fig, ax = plt.subplots()
        corrMatrix = final_genome_df.iloc[:, 0:-1].corr()
        sns.heatmap(corrMatrix, annot=True, cmap='BrBG')
        st.pyplot(fig)
        
        scatter_items = ['Genome size (Mb)', 'Chromosome','Plasmid','GC%', 'Genes', 'Proteins']
        scatter_items_selected = st.multiselect('Select 2 or 3 items', scatter_items, key = 1, default=['Genome size (Mb)', 'GC%'])
        # print(type(scatter_items_selected))
        
        if 0 <= len(scatter_items_selected) < 2:
            st.write('Select two or three items') 
        
        elif 2 <= len(scatter_items_selected) <= 3:
            fig = scatterplot_maker(scatter_items_selected, final_genome_df)
            st.plotly_chart(fig)
                
        else:
            st.write('Select two or three items')
        
        ## Data 4: Descriptive statistics of genomes (table data)
        st.markdown('**4. Descriptive statistics of complete genomes (Count values of taxonomic groups)**')
        # st.write('<Table data>')        
        df_menu = ['TaxID', 'Superkingdom','Phylum', 'Class', 'Order', 'Family', 'Genus','Species']
        df_choice = st.selectbox('Select one item', df_menu, key = 5, index = 7)
        
        try:
            if df_choice == 'Select data':
                st.write('Please select table data') 
            
            elif df_choice == 'TaxID':
                TaxID_series = final_genome_df['TaxID'].value_counts()    # count year values and save it to a Series dattype
                TaxID_df = pd.DataFrame({'TaxID': TaxID_series.index, 'Count': TaxID_series.values})
                st.dataframe(TaxID_df.style.highlight_max(axis=0))
                
                # ## download table data
                if st.button('Download the selected complete genome data as CSV', key = 23):
                    tmp_download_link = download_link(TaxID_df, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                    st.markdown(tmp_download_link, unsafe_allow_html=True)
                
            
            elif df_choice == 'Superkingdom':
                Superkingdom_series = final_genome_df['superkingdom'].value_counts()    # count year values and save it to a Series dattype
                Superkingdom_df = pd.DataFrame({'Superkingdom': Superkingdom_series.index, 'Count': Superkingdom_series.values})
                st.dataframe(Superkingdom_df.style.highlight_max(axis=0))
                
                # ## download table data
                if st.button('Download the selected complete genome data as CSV', key = 23):
                    tmp_download_link = download_link(Superkingdom_df, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                    st.markdown(tmp_download_link, unsafe_allow_html=True)
            
            elif df_choice == 'Phylum':
                Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_genome_df)
                st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                
                # ## download table data
                if st.button('Download the selected complete genome data as CSV', key = 23):
                    tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                    st.markdown(tmp_download_link, unsafe_allow_html=True)        
            
            elif df_choice == 'Class':
                Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_genome_df)
                st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                
                # ## download table data
                if st.button('Download the selected complete genome data as CSV', key = 23):
                    tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                    st.markdown(tmp_download_link, unsafe_allow_html=True)
                
            elif df_choice == 'Order':
                Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_genome_df)
                st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                
                # ## download table data
                if st.button('Download the selected complete genome data as CSV', key = 23):
                    tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                    st.markdown(tmp_download_link, unsafe_allow_html=True)
                
            elif df_choice == 'Family':
                Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_genome_df)
                st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                
                # ## download table data
                if st.button('Download the selected complete genome data as CSV', key = 23):
                    tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                    st.markdown(tmp_download_link, unsafe_allow_html=True)
                
            elif df_choice == 'Genus':
                Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_genome_df)
                st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                
                # ## download table data
                if st.button('Download the selected complete genome data as CSV', key = 23):
                    tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                    st.markdown(tmp_download_link, unsafe_allow_html=True)
                
            elif df_choice == 'Species':
                Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_genome_df) 
                st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                
                # ## download table data
                if st.button('Download the selected complete genome data as CSV', key = 23):
                    tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                    st.markdown(tmp_download_link, unsafe_allow_html=True)
        except:
            st.write('Please select genome(s)')   
        
        
    elif choice == '2. Selected by taxonomic lineage':
        st.write('')
        st.write('')
        st.success('**Complete genomes selected by taxonomic group(s)**')
        st.subheader('**Step 1**: Select complete genome(s) by taxonomic group(s)')
        
        col1,col2 = st.beta_columns(2)
    
        with col1:
            items_text = ['TaxID', 'superkingdom', 'phylum', 'class', 'order', 'family','genus', 'species']
            items_selected_text = st.selectbox('TaxID and taxonomic ranks', items_text, key = 6, index= 1)
        
        with col2:
            if items_selected_text == 'TaxID':
            
                taxID_list = list(set(final_genome_df['TaxID']))
                taxID_list.sort()
                subitems_selected = st.multiselect('Select a taxID or multiple taxIDs from ' + str(len(taxID_list)) + ' (maximum 10)', taxID_list, key = 2)
                                
            elif items_selected_text == 'superkingdom':
                    
                superkingdom_list = list(set(final_genome_df['superkingdom']))
                superkingdom_list = list(filter(None, superkingdom_list))
                superkingdom_list.sort()
                subitems_selected = st.multiselect('Select a superkingdom (from ' + str(len(superkingdom_list)) + ' superkingdoms)', superkingdom_list, key = 7)
        
            elif items_selected_text == 'phylum':
                
                phylum_list = list(set(final_genome_df['phylum']))
                phylum_list = list(filter(None, phylum_list))
                phylum_list.sort()
                subitems_selected = st.multiselect('Select a phylum or multiple phyla from ' + str(len(phylum_list)) + ' phyla (maximum 10)', phylum_list, key = 3)
            
            elif items_selected_text == 'class':
                class_list = list(set(final_genome_df['class']))
                class_list = list(filter(None, class_list))
                class_list.sort()
                subitems_selected = st.multiselect('Select a class or multiple class names from ' + str(len(class_list)) + ' classes (maximum 10)', class_list, key =4)
            
            elif items_selected_text == 'order':
                order_list = list(set(final_genome_df['order']))
                order_list = list(filter(None, order_list))
                order_list.sort()
                subitems_selected = st.multiselect('Select a order or multiple order names from ' + str(len(order_list)) + ' orders (maximum 10)', order_list, key = 5)
            
            elif items_selected_text == 'family':
                family_list = list(set(final_genome_df['family']))
                family_list = list(filter(None, family_list))
                family_list.sort()
                subitems_selected = st.multiselect('Select a family or multiple family names from ' + str(len(family_list)) + ' families (maximum 10)', family_list, key = 6)
            
            elif items_selected_text == 'genus':
                genus_list = list(set(final_genome_df['genus']))
                genus_list = list(filter(None, genus_list))
                genus_list.sort()
                subitems_selected = st.multiselect('Select a genus or multiple genus names from ' + str(len(genus_list)) + ' genera (maximum 10)', genus_list, key = 7)
            
            elif items_selected_text == 'species':
                species_list = list(set(final_genome_df['species']))
                species_list = list(filter(None, species_list))
                species_list.sort()
                subitems_selected = st.multiselect('Select a species or multiple species names from ' + str(len(species_list)) + ' species (maximum 10)', species_list, key = 8)
            
        # final_all_selected_genome, total_no_genomes, selected_no_genomes = making_all_select_genome_df(final_genome_df, items_selected_text, subitems_selected)
        output = making_all_select_genome_df(final_genome_df, items_selected_text, subitems_selected)
        if len(output) == 3:
            final_all_selected_genome = output[0]
            st.write(output[0])
            st.write('The selected genomes: ', output[2], '/', output[1])
            
            if st.button('Download the selected complete genome data as CSV', key = 23):
                tmp_download_link = download_link(final_all_selected_genome, 'selectGenome_{}_table.csv'.format(items_selected_text), 'Click here to download your data!')
                st.markdown(tmp_download_link, unsafe_allow_html=True)
            
        else:
            st.write(output)
        
        try:                
            # all_selected_genome = method_2_descript(final_genome_df)
            
            if final_all_selected_genome.shape[0] != 0:        
                
                ## Step 2: 
                st.subheader('**Step 2**: Overview of the selected complete genomes')
                # st.write('Short Introducton')
                
                ## Data 1:
                st.markdown('**1. Number of the selected bacterial and/or archaeal genomes sequenced each year and submitted to NCBI.**')
                
                try:
                    year_index_df = barchart_maker(final_all_selected_genome, 'Release Date')
                    st.bar_chart(year_index_df)        
                except:
                    pass
                
                ## Data 2: Descriptive statistics of genomes
                st.markdown('**2. Descriptive statistics of complete genomes (boxplot)**')
                graph_menu = ['Genome size (Mb)', 'Chromosome', 'Plasmid', 'GC%', 'Genes', 'Proteins']
                graph_choice = st.selectbox('Boxplot data of genome size, No. of chromosome, plasmid, gene, and protein, and GC%', graph_menu, key = 8, index = 0)
                
                try:
                    if graph_choice == 'Select one':
                        st.write('Please select a boxplot')
                    
                    elif graph_choice == 'GC%':
                        fig = boxplot_maker(graph_choice, final_all_selected_genome)
                        st.plotly_chart(fig)
                        
                    elif graph_choice == 'Genome size (Mb)':
                        fig = boxplot_maker(graph_choice, final_all_selected_genome)
                        st.plotly_chart(fig)
                        
                    elif graph_choice == 'Chromosome':
                        fig = boxplot_maker(graph_choice, final_all_selected_genome)
                        st.plotly_chart(fig)
                        
                    elif graph_choice == 'Plasmid':
                        fig = boxplot_maker(graph_choice, final_all_selected_genome)
                        st.plotly_chart(fig)
                    
                    elif graph_choice == 'Genes':
                        fig = boxplot_maker(graph_choice, final_all_selected_genome)
                        st.plotly_chart(fig)
                        
                    elif graph_choice == 'Proteins':
                        fig = boxplot_maker(graph_choice, final_all_selected_genome)
                        st.plotly_chart(fig)
                        
                except:
                    pass
                
                ## Data 3: Descriptive statistics of genomes (scatter plot)
                st.markdown('**3. Descriptive statistics of complete genomes (Pearson correlation heatmap and scatter plot [2d or 3d])**')
                # st.write('<For a 3d scatter plot, please select three values>')
                st.write('Pearson correlation heatmap')
                
                # st.write('<For a 3d scatter plot, please select three values>')
                
                # fig, ax = plt.subplots()
                # corrMatrix = final_all_selected_genome.iloc[:, 0:-1].corr()
                # sns.heatmap(corrMatrix, annot=True, cmap='BrBG')
                # st.pyplot(fig)
                
                
                # fig, ax = plt.subplots()
                
                data_df = final_all_selected_genome.iloc[:, 0:-1]
                data_df['Chromosome'] = pd.to_numeric(final_all_selected_genome['Chromosome'], errors='coerce')
                data_df['Plasmid'] = pd.to_numeric(final_all_selected_genome['Plasmid'], errors='coerce')
                
                
                fig, ax = plt.subplots()
                mask = np.triu(np.ones_like(data_df.corr()))
                sns.heatmap(data_df.corr(), mask=mask, annot=True, cmap='BrBG')
                # sns.heatmap(data_df.corr(), mask=mask, annot=True)
                st.pyplot(fig)
                
                
                scatter_items = ['Genome size (Mb)', 'Chromosome','Plasmid','GC%', 'Genes', 'Proteins']
                scatter_items_selected = st.multiselect('Select 2 or 3 items', scatter_items, key = 9, default=['Genome size (Mb)', 'GC%'])
                # print(type(scatter_items_selected))
                
                ## plotting test
                
                # st.write('test')
                # x = final_all_selected_genome[scatter_items_selected[0]]
                # y = final_all_selected_genome[scatter_items_selected[1]]
                # st.write(y)
                # st.write(x)
                # fig, ax = plt.subplots(figsize=(8, 4))
                # # ax.plot(X, response, color='k', label='Regression model')
                # ax.scatter(x, y, edgecolor='k', facecolor='grey', alpha=0.7, label='Sample data')
                # # ols = linear_model.LinearRegression()
                # # model = ols.fit(x, y)
                # # r2 = model.score(x, y)
                # st.pyplot(fig)
                # # st.write(r2)
                # st.write('test')
                # plt.style.use('default')
                # plt.style.use('ggplot')

                # fig, ax = plt.subplots(figsize=(8, 4))

                # ax.plot(X, response, color='k', label='Regression model')
                # ax.scatter(X, y, edgecolor='k', facecolor='grey', alpha=0.7, label='Sample data')
                # ax.set_ylabel('Gas production (Mcf/day)', fontsize=14)
                # ax.set_xlabel('Porosity (%)', fontsize=14)
                # ax.text(0.8, 0.1, 'aegis4048.github.io', fontsize=13, ha='center', va='center',
                #         transform=ax.transAxes, color='grey', alpha=0.5)
                # ax.legend(facecolor='white', fontsize=11)
                # ax.set_title('$R^2= %.2f$' % r2, fontsize=18)
                
                if 0 <= len(scatter_items_selected) < 2:
                    st.write('Select two or three items') 
        
                elif 2 <= len(scatter_items_selected) <= 3:
                    fig = scatterplot_maker(scatter_items_selected, final_all_selected_genome)
                    st.plotly_chart(fig)
                else:
                    st.write('Select two or three items')
                
                ## Data 4: Descriptive statistics of genomes (table data)
                st.markdown('**4. Descriptive statistics of complete genomes (Count values of taxonomic groups)**')
                # st.write('<Table data>')        
                df_menu = ['TaxID', 'Superkingdom','Phylum', 'Class', 'Order', 'Family', 'Genus','Species']
                df_choice = st.selectbox('Select one', df_menu, key = 14, index = 7)
                
                try:
                    if df_choice == 'Select data':
                        st.write('Please select table data') 
                    
                    elif df_choice == 'TaxID':
                        TaxID_series = final_all_selected_genome['TaxID'].value_counts()    # count year values and save it to a Series dattype
                        TaxID_df = pd.DataFrame({'TaxID': TaxID_series.index, 'Count': TaxID_series.values})
                        st.dataframe(TaxID_df.style.highlight_max(axis=0))
                        
                        ## download table data
                        if st.button('Download the table data as CSV', key = 23):
                            tmp_download_link = download_link(TaxID_df, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                            st.markdown(tmp_download_link, unsafe_allow_html=True)
                    
                    elif df_choice == 'Superkingdom':
                        Superkingdom_series = final_all_selected_genome['superkingdom'].value_counts()    # count year values and save it to a Series dattype
                        Superkingdom_df = pd.DataFrame({'Superkingdom': Superkingdom_series.index, 'Count': Superkingdom_series.values})
                        st.dataframe(Superkingdom_df.style.highlight_max(axis=0))
                        
                        ## download table data
                        if st.button('Download the selected complete genome data as CSV', key = 24):
                            tmp_download_link = download_link(Superkingdom_df, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                            st.markdown(tmp_download_link, unsafe_allow_html=True)
                        
                            
                    elif df_choice == 'Phylum':
                        Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                        st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                        
                        ## download table data
                        if st.button('Download the selected complete genome data as CSV', key = 25):
                            tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                            st.markdown(tmp_download_link, unsafe_allow_html=True)
                        
                        
                    elif df_choice == 'Class':
                        Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                        st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                        
                        ## download table data
                        if st.button('Download the selected complete genome data as CSV', key = 26):
                            tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                            st.markdown(tmp_download_link, unsafe_allow_html=True)
                        
                    elif df_choice == 'Order':
                        Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                        st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                        
                        ## download table data
                        if st.button('Download the selected complete genome data as CSV', key = 27):
                            tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                            st.markdown(tmp_download_link, unsafe_allow_html=True)
                        
                    elif df_choice == 'Family':
                        Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                        st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                        
                        ## download table data
                        if st.button('Download the selected complete genome data as CSV', key = 28):
                            tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                            st.markdown(tmp_download_link, unsafe_allow_html=True)
                        
                    elif df_choice == 'Genus':
                        Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                        st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                        
                        ## download table data
                        if st.button('Download the selected complete genome data as CSV', key = 29):
                            tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                            st.markdown(tmp_download_link, unsafe_allow_html=True)
                        
                    elif df_choice == 'Species':
                        Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                        st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                        
                        ## download table data
                        if st.button('Download the selected complete genome data as CSV', key = 30):
                            tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                            st.markdown(tmp_download_link, unsafe_allow_html=True)
                except:
                    pass

        except:
            st.write("Plese select genome(s) for further analyses in the Step 2")

    
    else:
        st.write('')
        st.write('')
        st.success('**Complete genome(s) selected by genomic feature(s)**')
        st.subheader('**Step 1**: Select complete genomes by genomic feature(s)')
        
        try:
            
            ## multiselect
            items = ['Genome Size', 'No. of Chromosome','No. of Plasmid','GC%', 'No. of Genes', 'No. of Proteins']
            items_selected = st.multiselect('Select Genome features', items, key = 11)
            
            
            selected_features = {}
                    
            if 'Genome Size' in items_selected:
                # Add a slider to the sidebar:
                ## min and max values of "Genome Size" from the reset_data
                genomeSize_min = float(final_genome_df['Genome size (Mb)'].min())
                genomeSize_max = float(final_genome_df['Genome size (Mb)'].max())
                # print(genomeSize_max)
                add_slider_size = st.slider('Select a range of Genome Size (Mb)', genomeSize_min, genomeSize_max, (genomeSize_min, genomeSize_max), key = 1)
                selected_features['Genome size (Mb)'] = add_slider_size
            
            if 'No. of Chromosome' in items_selected:
                chromosome_max = int(final_genome_df['Chromosome'].max())
                chromosome_min = int(final_genome_df['Chromosome'].min())
                add_slider_chromosome = st.slider('Select a range of No. of Chromosome', chromosome_min, chromosome_max, (chromosome_min, chromosome_max), key = 2)
                selected_features['Chromosome'] = add_slider_chromosome
            
            if 'No. of Plasmid' in items_selected:
                plasmid_max = int(final_genome_df['Plasmid'].max())
                plasmid_min = int(final_genome_df['Plasmid'].min())
                add_slider_plasmid = st.slider('Select a range of No. of Plasmid', plasmid_min, plasmid_max, (plasmid_min, plasmid_max), key = 3)
                selected_features['Plasmid'] = add_slider_plasmid
                    
            if 'GC%' in items_selected:
                ## min and max values of "GC%" from the reset_data
                final_genome_df['GC%'] = final_genome_df['GC%'].astype(float) # object to float
                gc_max = float(final_genome_df['GC%'].max())
                gc_min = float(final_genome_df['GC%'].min())
                
                add_slider_gc = st.slider('Select a range of GC%', gc_min, gc_max, (gc_min, gc_max), key = 4)
                selected_features['GC%'] = add_slider_gc
            
            if 'No. of Genes' in items_selected:
                # final_genome_df['Genes'] = final_genome_df['Genes'].astype(int) # object to float
                genes_max = int(final_genome_df['Genes'].max())
                genes_min = int(final_genome_df['Genes'].min())
                
                add_slider_genes = st.slider('Select a range of No. of Genes', genes_min, genes_max, (genes_min, genes_max), key = 5)
                selected_features['Genes'] = add_slider_genes
            
            if 'No. of Proteins' in items_selected:
                proteins_max = int(final_genome_df['Proteins'].max())
                proteins_min = int(final_genome_df['Proteins'].min())
                add_slider_proteins = st.slider('Select a range of No. of Proteins', proteins_min, proteins_max, (proteins_min, proteins_max), key = 6)
                selected_features['Proteins'] = add_slider_proteins
            
            
            output_features = genome_selector_features(final_genome_df, selected_features)

            if len(output_features) != 3:
                st.write(output_features) 

            else:
                final_all_selected_genome = output_features[0]
                st.write(final_all_selected_genome)
                
                st.write('The selected genomes: ', output_features[2], '/', output_features[1])
            
                if st.button('Download the selected complete genome data as CSV', key = 23):
                    tmp_download_link = download_link(final_all_selected_genome, 'selectGenome_genomicFeatures_table.csv', 'Click here to download your data!')
                    st.markdown(tmp_download_link, unsafe_allow_html=True)
            
            
                if final_all_selected_genome.shape[0] != 0:  
                                
                    ## Step 2: 
                    st.subheader('**Step 2**: Overview of the selected complete genomes')
                    # st.write('Short Introducton')
                    
                    ## Data 1: barchart
                    st.markdown('**1. Number of the selected bacterial and/or archaeal genomes sequenced each year and submitted to NCBI.**')
                    
                    try:
                        year_index_df = barchart_maker(final_all_selected_genome, 'Release Date')
                        st.bar_chart(year_index_df) 
                    except:
                        pass
                    
                    ## Data 2: Descriptive statistics of genomes
                    st.markdown('**2. Descriptive statistics of complete genomes (boxplot)**')
                    graph_menu = ['Genome size (Mb)', 'Chromosome', 'Plasmid', 'GC%', 'Genes', 'Proteins']
                    graph_choice = st.selectbox('Boxplot data of genome size, No. of chromosome, plasmid, gene, and protein, and GC%', graph_menu, key = 12, index = 0)
                    
                    try:
                        if graph_choice == 'Select one':
                            st.write('Please select a boxplot')
                        
                        elif graph_choice == 'GC%':
                            fig = boxplot_maker(graph_choice, final_all_selected_genome)
                            st.plotly_chart(fig)
                        
                        elif graph_choice == 'Genome size (Mb)':
                            fig = boxplot_maker(graph_choice, final_all_selected_genome)
                            st.plotly_chart(fig)
                        
                        elif graph_choice == 'Chromosome':
                            fig = boxplot_maker(graph_choice, final_all_selected_genome)
                            st.plotly_chart(fig)
                            
                        elif graph_choice == 'Plasmid':
                            fig = boxplot_maker(graph_choice, final_all_selected_genome)
                            st.plotly_chart(fig)
                            
                        elif graph_choice == 'Genes':
                            fig = boxplot_maker(graph_choice, final_all_selected_genome)
                            st.plotly_chart(fig)
                            
                        elif graph_choice == 'Proteins':
                            fig = boxplot_maker(graph_choice, final_all_selected_genome)
                            st.plotly_chart(fig)
                    except:
                        pass
                    
                    ## Data 3: Descriptive statistics of genomes (scatter plot)
                    
                    st.markdown('**3. Descriptive statistics of complete genomes (Pearson correlation heatmap and scatter plot [2d or 3d])**')
                    # st.write('<For a 3d scatter plot, please select three values>')
                    st.write('Pearson correlation heatmap')
                    
                    data_df = final_all_selected_genome.iloc[:, 0:-1]
                    data_df['Chromosome'] = pd.to_numeric(final_all_selected_genome['Chromosome'], errors='coerce')
                    data_df['Plasmid'] = pd.to_numeric(final_all_selected_genome['Plasmid'], errors='coerce')
                    
                    
                    fig, ax = plt.subplots()
                    mask = np.triu(np.ones_like(data_df.corr()))
                    sns.heatmap(data_df.corr(), mask=mask, annot=True, cmap='BrBG')
                    # sns.heatmap(data_df.corr(), mask=mask, annot=True)
                    st.pyplot(fig)
                    
                    
                    scatter_items = ['Genome size (Mb)', 'Chromosome','Plasmid','GC%', 'Genes', 'Proteins']
                    scatter_items_selected = st.multiselect('Select 2 or 3 items', scatter_items, key = 12, default=['Genome size (Mb)', 'GC%'])
                    
                    # fig = scatterplot_maker_3(('Genome size (Mb)', 'GC%'), final_all_selected_genome)
                    # st.pyplot(fig)
                    
                    if 0 <= len(scatter_items_selected) < 2:
                        st.write('Select two or three items') 
                        # fig = scatterplot_maker(('Genome size (Mb)', 'GC%'), final_all_selected_genome)
                        # st.pyplot(fig)
        
                    elif len(scatter_items_selected) == 2:
                        fig = scatterplot_maker(scatter_items_selected, final_all_selected_genome)
                        # st.pyplot(fig)
                        
                        st.plotly_chart(fig)

                    elif len(scatter_items_selected) == 3:
                        fig = scatterplot_maker(scatter_items_selected, final_all_selected_genome)
                        # st.pyplot(fig)
                        st.plotly_chart(fig)

                    else:
                        st.write('Select two or three items')
                    
                                                
                    ## Data 4: Descriptive statistics of genomes (table data)
                    st.markdown('**4. Descriptive statistics of complete genomes (Count values of taxonomic groups)**')
                    # st.write('<Table data>')        
                    df_menu = ['TaxID', 'Superkingdom','Phylum', 'Class', 'Order', 'Family', 'Genus','Species']
                    df_choice = st.selectbox('Select one item', df_menu, key = 15, index= 7)
                    
                    try:
                        if df_choice == 'Select data':
                            st.write('Please select table data') 
                        
                        elif df_choice == 'TaxID':
                            TaxID_series = final_all_selected_genome['TaxID'].value_counts()    # count year values and save it to a Series dattype
                            TaxID_df = pd.DataFrame({'TaxID': TaxID_series.index, 'Count': TaxID_series.values})
                            st.dataframe(TaxID_df.style.highlight_max(axis=0))
                            
                            ## download table data
                            if st.button('Download the table data as CSV', key = 23):
                                tmp_download_link = download_link(TaxID_df, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                                st.markdown(tmp_download_link, unsafe_allow_html=True)
                        
                        elif df_choice == 'Superkingdom':
                            Superkingdom_series = final_all_selected_genome['superkingdom'].value_counts()    # count year values and save it to a Series dattype
                            Superkingdom_df = pd.DataFrame({'Superkingdom': Superkingdom_series.index, 'Count': Superkingdom_series.values})
                            st.dataframe(Superkingdom_df.style.highlight_max(axis=0))
                            
                            ## download table data
                            if st.button('Download the selected complete genome data as CSV', key = 24):
                                tmp_download_link = download_link(Superkingdom_df, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                                st.markdown(tmp_download_link, unsafe_allow_html=True)
                            
                                
                        elif df_choice == 'Phylum':
                            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                            
                            ## download table data
                            if st.button('Download the selected complete genome data as CSV', key = 25):
                                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                                st.markdown(tmp_download_link, unsafe_allow_html=True)
                            
                            
                        elif df_choice == 'Class':
                            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                            
                            ## download table data
                            if st.button('Download the selected complete genome data as CSV', key = 26):
                                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                                st.markdown(tmp_download_link, unsafe_allow_html=True)
                            
                        elif df_choice == 'Order':
                            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                            
                            ## download table data
                            if st.button('Download the selected complete genome data as CSV', key = 27):
                                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                                st.markdown(tmp_download_link, unsafe_allow_html=True)
                            
                        elif df_choice == 'Family':
                            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                            
                            ## download table data
                            if st.button('Download the selected complete genome data as CSV', key = 28):
                                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                                st.markdown(tmp_download_link, unsafe_allow_html=True)
                            
                        elif df_choice == 'Genus':
                            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                            
                            ## download table data
                            if st.button('Download the selected complete genome data as CSV', key = 29):
                                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                                st.markdown(tmp_download_link, unsafe_allow_html=True)
                            
                        elif df_choice == 'Species':
                            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, final_all_selected_genome)
                            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
                            
                            ## download table data
                            if st.button('Download the selected complete genome data as CSV', key = 30):
                                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                                st.markdown(tmp_download_link, unsafe_allow_html=True)
                    except:
                        pass
                    
                    
        except:
            pass
            # st.write("If you cann't see a table (with selected genomes or empty), please hit the 'Refresh Data' button")
            # if st.button('Refresh Data'):
            #     rerun()
                # raise RerunException(st.ScriptRequestQueue.RerunData(None))
                # raise RerunException
                # raise RerunException()
            
            # st.write("Plese select genome(s) for further analyses in the Step 2")             

if __name__ == '__main__':
    main()