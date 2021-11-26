from numpy.lib.function_base import select
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
    
    ncbi_data.setFilter(first_menu, {'menu':selected_text} )

    # Subitems for Taxonomic Ranks
    
    subitems_selected = st.sidebar.multiselect(f'Select a {selected_text}', ncbi_data.tax_items[selected_text])
    ncbi_data.setFilter(first_menu, {'menu':selected_text, 'values': subitems_selected})
        
    # Genome Features
    # second_menu = 'Genome Features'
    # items = ['Genome Size', 'No. of Chromosome','No. of Plasmid','GC%', 'No. of Genes', 'No. of Proteins']
    
    for title, menu in ncbi_data.size_menus.items(): 
        checked = st.sidebar.checkbox(title)
        (min_v, max_v) = menu['range']
        add_slider_size = st.sidebar.slider(menu['title'], min_v, max_v, (min_v, max_v))
  
        filter_item = {k:v for k, v in menu.items()}
        filter_item['values'] = add_slider_size
        filter_item['checked'] = checked
        ncbi_data.setFilter(title, filter_item)
    
    
    # selected_features = st.sidebar.multiselect('Select Genome features', items, )
    
    # print(selected_features)

    # filters = []
    # for menu_title in selected_features:
    #     menu_item = ncbi_data.size_menus[menu_title]
    #     (min_v, max_v) = menu_item['range']
    #     add_slider_size = st.sidebar.slider(menu_item['title'], min_v, max_v, (min_v, max_v))
    #     filter_item = {k:v for k, v in menu_item.items()}
    #     filter_item['values'] = add_slider_size
    #     filters.append(filter_item)
    # ncbi_data.setFilter(second_menu, filters)
    
    # clicked = st.sidebar.button('Apply')

    print('sidebar rendering complete.')
    # return clicked
    

# @st.cache(allow_output_mutation=True) 
def initialize_data():
    pd.set_option('mode.chained_assignment', None)

    ncbi_data = NCBIdata()
    ncbi_df = ncbi_data.load()

    return ncbi_data

def display_filters(ncbi_data):
    st.header('Selected Complete Genomes')
    for title, filter in ncbi_data.filters.items():
        if title == 'Taxonomic Ranks':
            if filter['values']: 
                filter_menu = filter['menu']
                filter_values = str(', '.join(filter['values']))
                st.markdown(f'**{filter_menu}** : {filter_values}')
        elif 'checked' in filter and filter['checked']:
            filter_values = str(filter['values'][0]) + ' - ' + str(filter['values'][1])
            st.markdown(f"**{filter['menu']}**: {filter_values}")
    
    st.dataframe(data=ncbi_data.filtered_df)
        
    col1, col2 = st.columns(2)
    with col1:
        st.write(ncbi_data.filtered_df.shape[0], 'rows x', ncbi_data.filtered_df.shape[1], 'columns')
    with col2:
        tmp_download_link = download_link(ncbi_data.filtered_df, 'genome.csv', 'Download (csv)')
        st.markdown(tmp_download_link, unsafe_allow_html=True)
    

def analysis_num_submission(ncbi_df):
    ## Data 1:
    st.header('Number of Submissions to NCBI')
    year_index_df = barchart_maker(ncbi_df, 'Release Date')
    
    st.bar_chart(year_index_df)

def analysis_descriptive(ncbi_df):
    ## Data 2: Descriptive statistics of genomes
    st.header('Descriptive Statistics')
    graph_menu = ['Genome size (Mb)', 'Chromosome', 'Plasmid', 'GC%', 'Genes', 'Proteins']
    graph_choice = st.selectbox('Boxplot data of genome size, No. of chromosome, plasmid, gene, and protein, and GC%', graph_menu, key = 1, index = 0)
    
    with st.spinner('Please wait...'):
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

def analysis_heatmap(ncbi_df):
    ## Data 3: Descriptive statistics of genomes (scatter plot)
    st.header('Pearson Correlation Heatmap')
    # st.write('<For a 3d scatter plot, please select three values>')
    # st.write('Pearson correlation heatmap')
    fig, ax = plt.subplots()
    corrMatrix = ncbi_df.iloc[:, 0:-1].corr()
    sns.heatmap(corrMatrix, annot=True, cmap='BrBG')
    st.pyplot(fig)
    
def analysis_scatterplot(ncbi_df):
    st.header('Scatter plot [2d or 3d])')
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
    st.header('Distribution by Taxonomic Groups')
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
    st.image('./images/Picture7.png', use_column_width=True)

    with st.spinner('Downloading data from NCBI...'):
        ncbi_data = initialize_data()
        ncbi_df = ncbi_data.genome_df

    create_sidebar(ncbi_data)
    # print('apply_clicked = ', apply_clicked)

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
    st.write(f'**Data sourece:** GenBank prokaryotes.txt file downloaded **_{d1}_** **(a total of {no_genomes} completed genomes) (https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/)**')
    
    # if st.checkbox('Show raw prokaryote genome table (5 genomes only) and download all genome table'):
    #     st.dataframe(ncbi_df.head())
                
    #     if st.button('Download complete genome data as CSV', key = 1):
    #         tmp_download_link = download_link(ncbi_df, 'Genome_infor.csv', 'Click here to download your data!')
    #         st.markdown(tmp_download_link, unsafe_allow_html=True)

    
    st.write('')
    st.write('')

    display_filters(ncbi_data)

        

    analysis_num_submission(ncbi_data.filtered_df)
    analysis_descriptive(ncbi_data.filtered_df)
    analysis_heatmap(ncbi_data.filtered_df)
    analysis_scatterplot(ncbi_data.filtered_df)
    analysis_section4(ncbi_data.filtered_df)
    
    print('page rendering complete.')
    

    
    
     

if __name__ == '__main__':
    main()
