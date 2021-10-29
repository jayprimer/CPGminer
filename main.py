import streamlit as st
import pandas as pd
import pprintpp as pp
from ncbi import NCBIdata, download_link, count_tableMaker_groupby
from datetime import date
import base64
import matplotlib.pyplot as plt
import seaborn as sns


from charts import barchart_maker, boxplot_maker, scatterplot_maker, scatterplot_maker_3

def create_sidebar(ncbi_data):
    print('sidebar rendering started...')
    st.sidebar.image('./images/Picture10.png')

    # Taxonomic Ranks
    first_menu = 'Taxonomic Ranks'
    items_text = ['TaxID', 'superkingdom', 'phylum', 'class', 'order', 'family','genus', 'species']
    selected_text = st.sidebar.selectbox('TaxID and taxonomic ranks', items_text, index= 1)
    
    ncbi_data.setFilter(first_menu, selected_text)

    # Subitems for Taxonomic Ranks
    subitems_selected = st.sidebar.multiselect(f'Select a {selected_text}', ncbi_data.tax_items[selected_text])
                
        
    # Genome Features
    second_menu = 'Genome Features'
    items = ['Genome Size', 'No. of Chromosome','No. of Plasmid','GC%', 'No. of Genes', 'No. of Proteins']
    selected_features = st.sidebar.multiselect('Select Genome features', items, )
    
    ncbi_data.setFilter(second_menu, selected_features)

    ncbi_df = ncbi_data.genome_df

    menu_items = { 
        'Genome Size': {
            'menu': 'Genome Size',
            'title': 'Select a range of Genome Size (Mb)',
            'col_name': 'Genome size (Mb)',
            },
        'No. of Chromosome': {
            'menu': 'No. of Chromosome',
            'title': 'Select a range of No. of Chromosome',
            'col_name': 'Chromosome',
            },
        'No. of Plasmid': {
            'menu': 'No. of Plasmid',
            'title': 'Select a range of No. of Plasmid',
            'col_name': 'Plasmid',
            },
        'GC%': {
            'menu': 'GC%',
            'title': 'Select a range of GC%',
            'col_name': 'GC%',
            },
        'No. of Genes': {
            'menu': 'No. of Genes',
            'title': 'Select a range of No. of Genes',
            'col_name': 'Genes',
            },
        'No. of Proteins': {
            'menu': 'No. of Proteins',
            'title': 'Select a range of No. of Proteins',
            'col_name': 'Proteins',
            },
   

    }

    for menu in selected_features:
        menu_item = menu_items[menu]
        (min_v, max_v) = ncbi_data.get_range(ncbi_df, menu_item['col_name'])
        add_slider_size = st.sidebar.slider(menu_item['title'], min_v, max_v, (min_v, max_v))
        # selected_features[menu_item['col_name']] = add_slider_size
    
    print('sidebar rendering complete.')
    

@st.cache(allow_output_mutation=True) 
def initialize_data():
    pd.set_option('mode.chained_assignment', None)

    ncbi_data = NCBIdata()
    ncbi_df = ncbi_data.load()

    return ncbi_data


def analysis_section1(ncbi_df):
    ## Data 1:
    st.markdown('**1. Number of bacterial and archaeal genomes sequenced each year and submitted to NCBI.**')
    year_index_df = barchart_maker(ncbi_df, 'Release Date')
    
    st.bar_chart(year_index_df)

def analysis_section2(ncbi_df):
    ## Data 2: Descriptive statistics of genomes
    st.markdown('**2. Descriptive statistics of complete genomes (boxplot)**')
    graph_menu = ['Genome size (Mb)', 'Chromosome', 'Plasmid', 'GC%', 'Genes', 'Proteins']
    graph_choice = st.selectbox('Boxplot data of genome size, No. of chromosome, plasmid, gene, and protein, and GC%', graph_menu, key = 1, index = 0)
    
    if graph_choice == 'Select one':
        st.write('Please select a boxplot')
    
    elif graph_choice == 'GC%':
        fig = boxplot_maker(graph_choice, ncbi_df)
        st.plotly_chart(fig)
    
    elif graph_choice == 'Genome size (Mb)':
        fig = boxplot_maker(graph_choice, ncbi_df)
        st.plotly_chart(fig)
        
    elif graph_choice == 'Chromosome':
        fig = boxplot_maker(graph_choice, ncbi_df)
        st.plotly_chart(fig)
        
    elif graph_choice == 'Plasmid':
        fig = boxplot_maker(graph_choice, ncbi_df)
        st.plotly_chart(fig)
            
    elif graph_choice == 'Genes':
        fig = boxplot_maker(graph_choice, ncbi_df)
        st.plotly_chart(fig)
        
    elif graph_choice == 'Proteins':
        fig = boxplot_maker(graph_choice, ncbi_df)
        st.plotly_chart(fig)

def analysis_section3(ncbi_df):
    ## Data 3: Descriptive statistics of genomes (scatter plot)
    st.markdown('**3. Descriptive statistics of complete genomes (Pearson correlation heatmap and scatter plot [2d or 3d])**')
    # st.write('<For a 3d scatter plot, please select three values>')
    st.write('Pearson correlation heatmap')
    fig, ax = plt.subplots()
    corrMatrix = ncbi_df.iloc[:, 0:-1].corr()
    sns.heatmap(corrMatrix, annot=True, cmap='BrBG')
    st.pyplot(fig)
    
    scatter_items = ['Genome size (Mb)', 'Chromosome','Plasmid','GC%', 'Genes', 'Proteins']
    scatter_items_selected = st.multiselect('Select 2 or 3 items', scatter_items, key = 1, default=['Genome size (Mb)', 'GC%'])
    # print(type(scatter_items_selected))
    
    if 0 <= len(scatter_items_selected) < 2:
        st.write('Select two or three items') 
    
    elif 2 <= len(scatter_items_selected) <= 3:
        fig = scatterplot_maker(scatter_items_selected, ncbi_df)
        st.plotly_chart(fig)
            
    else:
        st.write('Select two or three items')

def analysis_section4(ncbi_df):
    ## Data 4: Descriptive statistics of genomes (table data)
    st.markdown('**4. Descriptive statistics of complete genomes (Count values of taxonomic groups)**')
    # st.write('<Table data>')        
    df_menu = ['TaxID', 'Superkingdom','Phylum', 'Class', 'Order', 'Family', 'Genus','Species']
    df_choice = st.selectbox('Select one item', df_menu, key = 5, index = 7)
    
    try:
        if df_choice == 'Select data':
            st.write('Please select table data') 
        
        elif df_choice == 'TaxID':
            TaxID_series = ncbi_df['TaxID'].value_counts()    # count year values and save it to a Series dattype
            TaxID_df = pd.DataFrame({'TaxID': TaxID_series.index, 'Count': TaxID_series.values})
            st.dataframe(TaxID_df.style.highlight_max(axis=0))
            
            # ## download table data
            if st.button('Download the selected complete genome data as CSV', key = 23):
                tmp_download_link = download_link(TaxID_df, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                st.markdown(tmp_download_link, unsafe_allow_html=True)
            
        
        elif df_choice == 'Superkingdom':
            Superkingdom_series = ncbi_df['superkingdom'].value_counts()    # count year values and save it to a Series dattype
            Superkingdom_df = pd.DataFrame({'Superkingdom': Superkingdom_series.index, 'Count': Superkingdom_series.values})
            st.dataframe(Superkingdom_df.style.highlight_max(axis=0))
            
            # ## download table data
            if st.button('Download the selected complete genome data as CSV', key = 23):
                tmp_download_link = download_link(Superkingdom_df, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                st.markdown(tmp_download_link, unsafe_allow_html=True)
        
        elif df_choice == 'Phylum':
            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, ncbi_df)
            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
            
            # ## download table data
            if st.button('Download the selected complete genome data as CSV', key = 23):
                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                st.markdown(tmp_download_link, unsafe_allow_html=True)        
        
        elif df_choice == 'Class':
            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, ncbi_df)
            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
            
            # ## download table data
            if st.button('Download the selected complete genome data as CSV', key = 23):
                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                st.markdown(tmp_download_link, unsafe_allow_html=True)
            
        elif df_choice == 'Order':
            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, ncbi_df)
            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
            
            # ## download table data
            if st.button('Download the selected complete genome data as CSV', key = 23):
                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                st.markdown(tmp_download_link, unsafe_allow_html=True)
            
        elif df_choice == 'Family':
            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, ncbi_df)
            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
            
            # ## download table data
            if st.button('Download the selected complete genome data as CSV', key = 23):
                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                st.markdown(tmp_download_link, unsafe_allow_html=True)
            
        elif df_choice == 'Genus':
            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, ncbi_df)
            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
            
            # ## download table data
            if st.button('Download the selected complete genome data as CSV', key = 23):
                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                st.markdown(tmp_download_link, unsafe_allow_html=True)
            
        elif df_choice == 'Species':
            Taxlank_df_dsend_reindex = count_tableMaker_groupby (df_choice, ncbi_df) 
            st.write(Taxlank_df_dsend_reindex.style.highlight_max(axis=0))
            
            # ## download table data
            if st.button('Download the selected complete genome data as CSV', key = 23):
                tmp_download_link = download_link(Taxlank_df_dsend_reindex, 'selectGenome_{}_count.csv'.format(df_choice), 'Click here to download your data!')
                st.markdown(tmp_download_link, unsafe_allow_html=True)
    except:
        st.write('Please select genome(s)')  

def main():
    ncbi_data = initialize_data()
    ncbi_df = ncbi_data.genome_df

    st.image('./images/Picture7.png', use_column_width=True)
    
    create_sidebar(ncbi_data)

    print('page rendering started...')

    ## Introduction; reference; data sources; etc
    with st.expander("See notes"):

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
    (no_genomes, no_columns) = ncbi_df.shape
    st.write(f'**Data sourece:** GenBank prokaryotes.txt file downloaded **_{d1}_** **(a total of {no_genomes} genomes) (https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/)**')
    
    if st.checkbox('Show raw prokaryote genome table (5 genomes only) and download all genome table'):
        st.dataframe(ncbi_df.head())
                
        if st.button('Download complete genome data as CSV', key = 1):
            tmp_download_link = download_link(ncbi_df, 'Genome_infor.csv', 'Click here to download your data!')
            st.markdown(tmp_download_link, unsafe_allow_html=True)

    
    st.write('')
    st.write('')

    analysis_section1(ncbi_df)
    analysis_section2(ncbi_df)
    analysis_section3(ncbi_df)
    analysis_section4(ncbi_df)
    
    print('page rendering complete.')
    

    
    
     

if __name__ == '__main__':
    main()
