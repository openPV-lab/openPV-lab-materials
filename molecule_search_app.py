import streamlit as st
import json
import os
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

# Configure Streamlit page
st.set_page_config(page_title="Molecule Search", layout="wide")

@st.cache_data
def load_material_data():
    """Load all material data from JSON files"""
    curr_dir = os.getcwd()
    material_data_dir = os.path.join(curr_dir, 'material_data')
    chemical_structure_dir = os.path.join(curr_dir, 'chemical_structure')
    
    materials = []
    
    if not os.path.exists(material_data_dir):
        st.error(f"Material data directory not found: {material_data_dir}")
        return materials, chemical_structure_dir
    
    json_files = [f for f in os.listdir(material_data_dir) if f.endswith('.json')]
    
    for json_file in json_files:
        try:
            with open(os.path.join(material_data_dir, json_file), 'r') as f:
                data = json.load(f)
                data['filename'] = json_file
                materials.append(data)
        except Exception as e:
            st.warning(f"Error loading {json_file}: {e}")
    
    return materials, chemical_structure_dir

def get_all_searchable_terms(materials):
    """Get all searchable terms (names and synonyms) with their corresponding materials"""
    search_dict = {}
    
    for material in materials:
        name = material.get('name', 'Unknown')
        # Add main name
        search_dict[name] = material
        
        # Add synonyms
        synonyms = material.get('synonyms', [])
        if isinstance(synonyms, list):
            for synonym in synonyms:
                if synonym and synonym.strip():
                    search_dict[synonym] = material
    
    return search_dict

def search_materials(materials, query):
    """Search materials by name and synonyms"""
    if not query:
        return []
    
    query_lower = query.lower()
    results = []
    
    for material in materials:
        # Search in name
        name = material.get('name', '').lower()
        if query_lower in name:
            results.append(material)
            continue
        
        # Search in synonyms
        synonyms = material.get('synonyms', [])
        if isinstance(synonyms, list):
            for synonym in synonyms:
                if query_lower in synonym.lower():
                    results.append(material)
                    break
    
    return results

def get_filtered_options(search_dict, query):
    """Get filtered options based on search query"""
    if not query:
        return list(search_dict.keys())
    
    query_lower = query.lower()
    filtered = [term for term in search_dict.keys() if query_lower in term.lower()]
    return sorted(filtered)

def get_molecule_image(material, chemical_structure_dir):
    """Get molecule image from .mol file or SMILES"""
    # Try .mol file first
    mol_file = material.get('mol_file')
    if mol_file:
        mol_path = os.path.join(chemical_structure_dir, mol_file)
        if os.path.exists(mol_path):
            try:
                mol = Chem.MolFromMolFile(mol_path)
                if mol:
                    return Draw.MolToImage(mol, size=(400, 400))
            except:
                pass
    
    # Fallback to SMILES
    smiles = material.get('smiles')
    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Draw.MolToImage(mol, size=(400, 400))
        except:
            pass
    
    return None

def display_material_info(material):
    """Display material information in a formatted way"""
    st.subheader(f"📄 Material Information")
    
    # Create a formatted display of the JSON data
    display_data = {}
    
    for key, value in material.items():
        if key == 'filename':
            continue
        
        if isinstance(value, list):
            display_data[key.replace('_', ' ').title()] = ', '.join(map(str, value))
        elif isinstance(value, dict):
            display_data[key.replace('_', ' ').title()] = str(value)
        else:
            display_data[key.replace('_', ' ').title()] = value
    
    # Display as a nice table
    df = pd.DataFrame(list(display_data.items()), columns=['Property', 'Value'])
    st.dataframe(df, use_container_width=True, hide_index=True)
    
    # Show raw JSON in expander
    with st.expander("View Raw JSON Data"):
        st.json(material)

def main():
    st.title("🔍 Molecule Search App")
    st.markdown("Search for molecules by name or synonyms and view their structure and properties.")
    
    # Load data
    materials, chemical_structure_dir = load_material_data()
    
    if not materials:
        st.error("No material data found. Please check the material_data directory.")
        return
    
    # Get all searchable terms
    search_dict = get_all_searchable_terms(materials)
    all_options = sorted(search_dict.keys())
    
    # Search interface
    st.sidebar.header("Search")
    
    # Text input for search
    search_query = st.sidebar.text_input("Type to search:", placeholder="e.g. ITIC, BTP-4F")
    
    # Get filtered options based on search query
    filtered_options = get_filtered_options(search_dict, search_query)
    
    # Dropdown with filtered options
    if filtered_options:
        selected_option = st.sidebar.selectbox(
            "Select molecule:",
            options=[""] + filtered_options,
            index=0,
            help="Choose from matching results"
        )
    else:
        selected_option = st.sidebar.selectbox(
            "Select molecule:",
            options=["No matches found"],
            index=0,
            disabled=True
        )
    
    # Determine which material to display
    selected_material = None
    if selected_option and selected_option != "No matches found":
        selected_material = search_dict.get(selected_option)
    
    # Display results count
    if search_query:
        st.sidebar.markdown(f"**Found {len(filtered_options)} matches**")
    else:
        st.sidebar.markdown(f"**{len(all_options)} materials available**")
    
    # Main content area
    if not selected_material:
        if search_query and not filtered_options:
            st.warning(f"No materials found matching '{search_query}'")
        else:
            st.info("Enter a search term or select a material from the dropdown to get started.")
        
        # Show all available materials as a fallback
        st.subheader("Available Materials:")
        material_names = [mat.get('name', 'Unknown') for mat in materials]
        st.write(", ".join(sorted(material_names)))
    
    else:
        # Display selected material
        name = selected_material.get('name', 'Unknown')
        
        st.header(f"🧪 {name}")
        
        # Show if this was found via synonym
        if selected_option != name:
            st.info(f"Found via synonym: **{selected_option}**")
        
        # Create two columns for layout
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("🧬 Molecular Structure")
            
            # Get and display molecule image
            mol_image = get_molecule_image(selected_material, chemical_structure_dir)
            
            if mol_image:
                st.image(mol_image, caption=f"Structure of {name}")
                
                # Show source information
                mol_file = selected_material.get('mol_file')
                if mol_file and os.path.exists(os.path.join(chemical_structure_dir, mol_file)):
                    st.caption(f"📁 Source: {mol_file}")
                else:
                    st.caption("📝 Generated from SMILES")
            else:
                st.error("Could not generate molecular structure")
                st.info("No .mol file or valid SMILES found")
        
        with col2:
            display_material_info(selected_material)

if __name__ == "__main__":
    main()
