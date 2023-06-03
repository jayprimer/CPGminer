import streamlit as st
import pandas as pd
from scipy import stats
import plotly.express as px
import seaborn as sns


@st.cache_data(persist="disk")  
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

@st.cache_data(persist="disk")
def boxplot_maker(graph_choice, genome_dataframe):
    
    '''
    Function to draw a boxplot of the genomic features selected
    argument 1: graph_choice
    argument 2: genome dataframe
    '''
    fig = px.box(genome_dataframe, y=graph_choice, points='all')
    return fig

@st.cache_data(persist="disk")
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
    
@st.cache_data(persist="disk")
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