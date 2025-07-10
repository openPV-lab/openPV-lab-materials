import streamlit as st
import json
import os
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

# Configure Streamlit page
st.set_page_config(page_title="Molecule Search", layout="wide")

@st.cache_data
def load_chemical_data():
    """Load all material data from JSON files"""
    curr_dir = os.getcwd()
    chemical_data_dir = os.path.join(curr_dir, 'chemical_data')
    chemical_structure_dir = os.path.join(curr_dir, 'chemical_structure')
    
    materials = []
    
    if not os.path.exists(chemical_data_dir):
        st.error(f"Material data directory not found: {chemical_data_dir}")
        return materials, chemical_structure_dir
    
    json_files = [f for f in os.listdir(chemical_data_dir) if f.endswith('.json')]
    
    for json_file in json_files:
        try:
            with open(os.path.join(chemical_data_dir, json_file), 'r') as f:
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
        iupac_name = material.get('iupac_name', '')
        
        # Add main name
        search_dict[name] = material
        
        # Add IUPAC name if available and not empty
        if iupac_name and iupac_name.strip():
            search_dict[iupac_name] = material
        
        # Add synonyms
        synonyms = material.get('synonyms', [])
        if isinstance(synonyms, list):
            for synonym in synonyms:
                if synonym and synonym.strip():
                    # Don't duplicate if synonym is same as name or IUPAC name
                    if synonym != name and synonym != iupac_name:
                        search_dict[synonym] = material
    
    return search_dict

def search_materials(materials, query):
    """Search materials by name, IUPAC name, and synonyms"""
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
        
        # Search in IUPAC name
        iupac_name = material.get('iupac_name', '').lower()
        if iupac_name and query_lower in iupac_name:
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

def get_all_classifications(materials):
    """Get all unique classifications from materials"""
    classifications = set()
    for material in materials:
        material_classifications = material.get('classification', [])
        if isinstance(material_classifications, list):
            classifications.update(material_classifications)
    return sorted(list(classifications))

def filter_materials_by_classification(materials, selected_classifications):
    """Filter materials by classification"""
    if not selected_classifications or "All" in selected_classifications:
        return materials
    
    filtered_materials = []
    for material in materials:
        material_classifications = material.get('classification', [])
        if isinstance(material_classifications, list):
            # Check if any of the material's classifications match the selected ones
            if any(cls in selected_classifications for cls in material_classifications):
                filtered_materials.append(material)
    
    return filtered_materials

def get_filtered_options(search_dict, query, materials, selected_classifications):
    """Get filtered options based on search query and classification"""
    # First filter by classification
    filtered_materials = filter_materials_by_classification(materials, selected_classifications)
    
    # Create search dict from filtered materials
    filtered_search_dict = {}
    for material in filtered_materials:
        name = material.get('name', 'Unknown')
        iupac_name = material.get('iupac_name', '')
        
        # Add main name
        filtered_search_dict[name] = material
        
        # Add IUPAC name if available and not empty
        if iupac_name and iupac_name.strip():
            filtered_search_dict[iupac_name] = material
        
        # Add synonyms
        synonyms = material.get('synonyms', [])
        if isinstance(synonyms, list):
            for synonym in synonyms:
                if synonym and synonym.strip():
                    # Don't duplicate if synonym is same as name or IUPAC name
                    if synonym != name and synonym != iupac_name:
                        filtered_search_dict[synonym] = material

    # Then filter by search query
    if not query:
        return list(filtered_search_dict.keys())
    
    query_lower = query.lower()
    filtered = [term for term in filtered_search_dict.keys() if query_lower in term.lower()]
    return sorted(filtered)

def calculate_image_size(container_width_fraction=0.8, max_size=800, min_size=300):
    """Calculate optimal image size based on screen/container width"""
    # Streamlit's default column width is roughly 700px for half-width columns
    # We'll estimate based on typical screen sizes and column layout
    base_width = 600 if container_width_fraction > 0.6 else 400
    
    # Scale with container width but respect min/max bounds
    calculated_size = int(base_width * container_width_fraction)
    size = max(min_size, min(calculated_size, max_size))
    
    return (size, size)

def get_molecule_image(material, chemical_structure_dir, image_size=None):
    """Get molecule image from .mol file or SMILES with configurable size"""
    if image_size is None:
        image_size = calculate_image_size()
    
    # Try .mol file first
    mol_file = material.get('mol_file')
    if mol_file:
        mol_path = os.path.join(chemical_structure_dir, mol_file)
        if os.path.exists(mol_path):
            try:
                mol = Chem.MolFromMolFile(mol_path)
                if mol:
                    return Draw.MolToImage(mol, size=image_size)
            except:
                pass
    
    # Fallback to SMILES
    smiles = material.get('smiles')
    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Draw.MolToImage(mol, size=image_size)
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
    materials, chemical_structure_dir = load_chemical_data()
    
    if not materials:
        st.error("No material data found. Please check the chemical_data directory.")
        return
    
    # Get all searchable terms and classifications
    search_dict = get_all_searchable_terms(materials)
    all_options = sorted(search_dict.keys())
    all_classifications = get_all_classifications(materials)
    
    # Search interface
    st.sidebar.header("Search")
    
    # Text input for search
    search_query = st.sidebar.text_input("Type to search:", placeholder="e.g. ITIC, BTP-4F")
    
    # Classification filter
    st.sidebar.subheader("Filter by Classification")
    classification_options = ["All"] + all_classifications
    selected_classifications = st.sidebar.multiselect(
        "Select classifications:",
        options=classification_options,
        default=["All"],
        help="Filter materials by their classification type"
    )
    
    # Get filtered options based on search query and classification
    filtered_options = get_filtered_options(search_dict, search_query, materials, selected_classifications)
    
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
        # Get the material from the filtered results
        filtered_materials = filter_materials_by_classification(materials, selected_classifications)
        filtered_search_dict = get_all_searchable_terms(filtered_materials)
        selected_material = filtered_search_dict.get(selected_option)
    
    # Display results count
    filtered_materials_count = len(filter_materials_by_classification(materials, selected_classifications))
    if search_query:
        st.sidebar.markdown(f"**Found {len(filtered_options)} matches**")
    else:
        st.sidebar.markdown(f"**{len(filtered_options)} materials available**")
    
    if "All" not in selected_classifications:
        st.sidebar.markdown(f"*({filtered_materials_count} total after classification filter)*")
    
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
            
            # Calculate optimal image size for the column
            image_size = calculate_image_size(container_width_fraction=0.9, max_size=600, min_size=400)
            
            # Get and display molecule image
            mol_image = get_molecule_image(selected_material, chemical_structure_dir, image_size)
            
            if mol_image:
                st.image(mol_image, caption=f"Structure of {name}", use_container_width=True)
                
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
