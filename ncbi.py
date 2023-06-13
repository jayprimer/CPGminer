import streamlit as st
import pandas as pd
import pprint
from ete3 import NCBITaxa
import base64
import re
from datetime import date
import os
import numpy as np

def change_ftp(ftp_path):
    '''
    function to change from ftp: to https: of the downloadable table
    :param ftp_path: a string of the ftp path
    :return:
    '''
    new_ftp_path = 'https' + ftp_path[3:]  
    return new_ftp_path

## download data
@st.cache_data(persist="disk")
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

@st.cache_data(persist="disk")
def count_tableMaker_groupby (Taxlank, genome_df):

    '''
    Function to make a dataframe of the count values of taxonomic groups of the selected genomes
    argument: Taxlank
    argument: final_genome_df 
    '''
    
    if Taxlank == 'Phylum':
        groupby_list = ['superkingdom','phylum']
    elif Taxlank == 'Class':
        groupby_list = ['superkingdom','phylum', 'class']
    elif Taxlank == 'Order':
        groupby_list = ['superkingdom','phylum', 'class', 'order']    
    elif Taxlank == 'Family':
        groupby_list = ['superkingdom','phylum', 'class', 'order', 'family']      
    elif Taxlank == 'Genus':
        groupby_list = ['superkingdom','phylum', 'class', 'order', 'family', 'genus'] 
    elif Taxlank == 'Species':
        groupby_list = ['superkingdom','phylum', 'class', 'order', 'family', 'genus', 'species'] 
                    
    Taxlank_series = genome_df.groupby(groupby_list).size()
    Taxlank_df = Taxlank_series.to_frame()
    Taxlank_df_reindex = Taxlank_df.reset_index()
    Taxlank_df_final = Taxlank_df_reindex.rename(columns={0: 'Count'})
    Taxlank_df_dsend = Taxlank_df_final.sort_values(by='Count', ascending=False)
    Taxlank_df_dsend_reindex = Taxlank_df_dsend.reset_index(drop=True)
    return Taxlank_df_dsend_reindex

@st.cache_data(persist="disk")
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

@st.cache_data(persist="disk")
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



class NCBIdata:
    def __init__(self):
        self.url = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
        self.tax_item_texts = ['TaxID', 'superkingdom', 'phylum', 'class', 'order', 'family','genus', 'species']
        self.column_names = ['#Organism/Name', 'TaxID', 'Size (Mb)', 'GC%', 'Replicons', 'Genes', 'Proteins', 'Release Date', 'Status', 'FTP Path']
        
        self.tax_items = {}
        self.filters = {}

        self.genome_df = None
        self.filtered_df = None

        self.cache_path = './cache'

        self.size_menus = { 
            'Genome Size': {
                'menu': 'Genome Size',
                'title': 'Select a range of Genome Size (Mb)',
                'col_name': 'Genome size (Mb)',
                },
            'GC%': {
                'menu': 'GC%',
                'title': 'Select a range of GC%',
                'col_name': 'GC%',
                },
            'Number of Chromosomes': {
                'menu': 'Number of Chromosomes',
                'title': 'Select a range of number of Chromosomes',
                'col_name': 'Chromosome',
                },
            'Number of Plasmids': {
                'menu': 'Number of Plasmids',
                'title': 'Select a range of number of Plasmids',
                'col_name': 'Plasmid',
                },
            'Number of Genes': {
                'menu': 'Number of Genes',
                'title': 'Select a range of number of Genes',
                'col_name': 'Genes',
                },
            'Number of Proteins': {
                'menu': 'Number of Proteins',
                'title': 'Select a range of number of Proteins',
                'col_name': 'Proteins',
                },
       }

    def get_cache_filename(self):
        today = date.today()

        filename = 'genome_df' + today.strftime("-%Y-%m-%d") + ".feather"
        return os.path.join(self.cache_path, filename)

    def save_to_cache(self, overwrite=False): 
        cache_file = self.get_cache_filename()
        if (not os.path.exists(cache_file)) or (os.path.exists(cache_file) and overwrite):
            if not os.path.exists(self.cache_path):
                os.makedirs(self.cache_path)
            self.genome_df.to_feather(cache_file)

    def load_from_cache(self):
        cache_file = self.get_cache_filename()
        if os.path.exists(cache_file):
            print('load from cache...')
            genome_df = pd.read_feather(cache_file)
            genome_df.set_index('index', inplace=True)
            return genome_df
        else:
            return None
    
    def save_df(self, df, filename): 
        filepath = os.path.join(self.cache_path, filename + '.feather')
        df.to_feather(filepath)
    
    def load_df(self, filename):
        filepath = os.path.join(self.cache_path, filename + '.feather')
        df = pd.read_feather(filepath)
        try:
            df.set_index('index', inplace=True)
        except:
            pass

        return df


    def processing(self):

        data = self.load_df('original')
        complete_genomes = data.loc[data['Status'] == "Complete Genome"]    # Select only complete genomes
        genome_df = complete_genomes.reset_index(drop = True) # reset index from 0
        genome_df = genome_df.drop(['Status'], axis = 1) # remove the column "Status"
        genome_df = genome_df.rename(columns={"#Organism/Name": "Genome Name"}) # change the first column name
        self.save_df(genome_df.reset_index(), 'step1')

        self.genome_df = self.count_chro_plas(genome_df)
        self.save_df(self.genome_df.reset_index(), 'step2')
        
        self.genome_df = self.making_final_df(self.genome_df)
        self.save_df(self.genome_df.reset_index(), 'step3')

        self.genome_df = self.genome_df.dropna()  # remove missing values
        self.save_df(self.genome_df.reset_index(), 'step4')

        self.calc_tax_items()
        self.save_df(self.genome_df.reset_index(), 'step5')


    def load_from_ncbi(self):
        print('load from ncbi...')
        dtype_data = {
            '#Organism/Name':np.string_, 
            'TaxID':np.string_, 
            'Size (Mb)':np.string_, 
            'GC%':np.string_, 
            'Replicons':np.string_, 
            'Genes':np.string_,
            'Proteins':np.string_, 
            'Release Date':np.string_, 
            'Status':np.string_, 
            'FTP Path':np.string_
        }
        data = pd.read_table(self.url, usecols=self.column_names, dtype=dtype_data)
        # print(data.shape)
        # print(data.columns)
        # self.save_df(data, 'original')
        
        complete_genomes = data.loc[data['Status'] == "Complete Genome"]    # Select only complete genomes
        genome_df = complete_genomes.reset_index(drop = True) # reset index from 0
        genome_df = genome_df.drop(['Status'], axis = 1) # remove the column "Status"
        genome_df = genome_df.rename(columns={"#Organism/Name": "Genome Name"}) # change the first column name
 
        genome_df = self.count_chro_plas(genome_df)
        genome_df = self.making_final_df(genome_df)
        genome_df = genome_df.dropna()  # remove missing values

        genome_df = genome_df.reset_index()
        return genome_df

    def load(self):
        self.genome_df = self.load_from_cache()
        if self.genome_df is None:
            self.genome_df = self.load_from_ncbi()
            self.save_to_cache()
        
        # perform processing for this class
        self.calc_tax_items()
        self.calc_range_for_size_menus()
        return self.genome_df

    def calc_tax_items(self):
        for item_text in self.tax_item_texts:
            self.tax_items[item_text] = list(set(self.genome_df[item_text]))
            self.tax_items[item_text] = list(filter(None, self.tax_items[item_text]))
            self.tax_items[item_text].sort()
    
    def calc_range_for_size_menus(self):
        for title, menu in self.size_menus.items():
            (min_v, max_v) = self.get_range(self.genome_df, menu['col_name'])
            menu['range'] = (min_v, max_v)

    def get_range(self, df, col_name):
        min_value = float(df[col_name].min())
        max_value = float(df[col_name].max())
        return (min_value, max_value)

    def count_chro_plas(self, genome_df):
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
    
    def making_final_df(self, genome_df):
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
        final_genome_df['Genome size (Mb)'] = pd.to_numeric(final_genome_df['Genome size (Mb)'], errors='coerce') 
        final_genome_df['GC%'] = pd.to_numeric(final_genome_df['GC%'], errors='coerce') 
        final_genome_df['Genes'] = pd.to_numeric(final_genome_df['Genes'], errors='coerce') 
        final_genome_df['Proteins'] = pd.to_numeric(final_genome_df['Proteins'], errors='coerce')
        final_genome_df['Chromosome'] = pd.to_numeric(final_genome_df['Chromosome'], errors='coerce')
        final_genome_df['Plasmid'] = pd.to_numeric(final_genome_df['Plasmid'], errors='coerce')
        final_genome_df['TaxID'] = final_genome_df['TaxID'].astype(str)

        return final_genome_df

    def setFilter(self, key, value):
        self.filters[key] = value
        pprint.pprint(self.filters)
        self.apply_filter()

    def __repr__(self):
        # return pp.pformat(self.filters)
        return str(self.filters)


    def is_filter(self):
        """
        return true if filter is set
        """
        for title, filter in self.filters.items():
            if title == 'Taxonomic Ranks':
                if filter['values']: 
                    return True
            elif 'checked' in filter and filter['checked']:
                return True
                
    def apply_filter(self):
        # apply Taxonomic Ranks first
        print('applying filter...')
        if self.genome_df is not None:
            self.filtered_df = self.genome_df

            tax_filter = self.filters['Taxonomic Ranks']
            if tax_filter and ('values' in tax_filter): 
                print(tax_filter)
                if tax_filter['values']:
                    print('applying taxonomic ranks', tax_filter)

                    self.filtered_df = self.genome_df.loc[self.genome_df[tax_filter['menu']].isin(tax_filter['values'])]

                    print('filtered', self.filtered_df.shape)
                    print(self.filtered_df.tail())

            for title, filter in self.filters.items():
                if 'checked' in filter and filter['checked']:
                    col_name = filter['col_name']
                    filter_values = filter['values']
                    
                    self.filtered_df = self.filtered_df.loc[
                        self.filtered_df[col_name].between(filter_values[0], filter_values[1], inclusive='both') 
                    ]
            
            print('filtered', self.filtered_df.shape)
            print(self.filtered_df.tail())


    


        
        

