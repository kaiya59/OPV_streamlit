import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
import pickle
from rdkit.Chem import MolFromSmiles as s2m
from rdkit.Chem import MolToSmiles as m2s
# from rdkit.Chem import Draw

# Load data
df = pd.read_csv('predicted_pce1.csv')
with open('lookup_table.pkl', 'rb') as f:
    lookup_dict = pickle.load(f)

# Streamlit App
st.title('OPV Molecular Data Visualization')

# Sidebar for options
option = st.sidebar.selectbox(
    'Select an option:',
    ['Show Data', 'Show Statistics']
)

if option == 'Show Data':
    st.subheader('Data Preview')

    # Filter options
    filter_column = st.selectbox(
        'Select the column to filter by (or select "All" for no filter):',
        ['All', 'acceptor', 'tin_reagents', 'ring_b1', 'ring_b2', 'chain1', 'chain2']
    )
    ascending = st.checkbox('Sort PCE in ascending order', value=False)  # Default is descending (False)

    if filter_column != 'All':
        # Get unique values for selected column
        unique_values = df[filter_column].unique()
        selected_value = st.selectbox(f'Select a value for {filter_column}:', unique_values)
        
        # Filter and sort the data
        filtered_df = df[df[filter_column] == selected_value]
        sorted_df = filtered_df.sort_values(by='PCE', ascending=ascending)
        sorted_df.reset_index(inplace=True)
        sorted_df.drop(columns=['index'], inplace=True)
    else:
        sorted_df = df.sort_values(by='PCE', ascending=ascending)
        if ascending:
            sorted_df.reset_index(inplace=True)
            sorted_df.drop(columns=['index'], inplace=True)

    st.write(sorted_df)

    # Display table with an option to select a row
    selected_index = st.selectbox('Select a row to view molecule details:', sorted_df.index)

    # Extract the selected row data
    selected_row = sorted_df.iloc[selected_index]
    
    # Extract SMILES strings
    tin_reagent_smiles = selected_row['tin_reagents']
    ring_b1_smiles = selected_row['ring_b1']
    ring_b2_smiles = selected_row['ring_b2']
    chain1_smiles = selected_row['chain1']
    chain2_smiles = selected_row['chain2']
    
    # Display the molecules
    st.subheader('Molecule Structures:')
    st.write(f"PCE = {selected_row['PCE']}")
    
    st.write("Tin Reagent:")
    image_path = f"./img1/{lookup_dict[tin_reagent_smiles]}.png"
    st.image(image_path, use_column_width=True)
    
    st.write("Ring B1:")
    image_path = f"./img1/{lookup_dict[ring_b1_smiles]}.png"
    st.image(image_path, use_column_width=True)
    
    st.write("Ring B2:")
    image_path = f"./img1/{lookup_dict[ring_b2_smiles]}.png"
    st.image(image_path, use_column_width=True)

    if chain1_smiles != '-':
        st.write("Side Chain 1:")
        image_path = f"./img1/{lookup_dict[chain1_smiles]}.png"
        st.image(image_path, use_column_width=True)
    else:
        st.write("Side Chain 1: -")
    
    if chain2_smiles != '-':
        st.write("Side Chain 2:")
        image_path = f"./img/chain2_{lookup_dict[chain2_smiles]}.png"
        st.image(image_path, use_column_width=True)
    else:
        st.write("Side Chain 2: -")

elif option == 'Show Statistics':
    st.subheader('Statistics Options')

    histogram_type = st.selectbox(
        'Select the type of histogram:',
        ['All Data', 'By Acceptor']
    )

    if histogram_type == 'All Data':
        st.subheader('Histogram of PCE Values')
        fig, ax = plt.subplots()
        ax.hist(df['PCE'].values, bins=100, color='blue', alpha=0.5)
        # df['PCE'].plot.kde(ax=ax, secondary_y=True, color='blue', alpha=0.7)  # Density plot
        st.pyplot(fig)
        st.write(f"Maximum value of predicted PCE: {df['PCE'].max()}")
        st.write(f"Minimum value of predicted PCE: {df['PCE'].min()}")

    elif histogram_type == 'By Acceptor':
        acceptor = st.selectbox('Select an acceptor:', df['acceptor'].unique())
        filtered_df = df[df['acceptor'] == acceptor]

        if not filtered_df.empty:
            st.subheader(f'Histogram of PCE Values for Acceptor: {acceptor}')
            fig, ax = plt.subplots()
            ax.hist(filtered_df['PCE'].values, bins=100, color='green', alpha=0.5)
            # filtered_df['PCE'].plot.kde(ax=ax, secondary_y=True, color='green', alpha=0.7)  # Density plot
            st.pyplot(fig)
            st.write(f"Maximum value of predicted PCE: {filtered_df['PCE'].max()}")
            st.write(f"Minimum value of predicted PCE: {filtered_df['PCE'].min()}")
        else:
            st.warning(f'No data available for acceptor: {acceptor}')
