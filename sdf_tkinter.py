import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
import requests
import io
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import tempfile
import logging

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)

class CompoundFetcherApp:
    def __init__(self, root):
        logger.info("Initializing CompoundFetcherApp")
        self.root = root
        self.root.title("Compound Data Fetcher")

        # Configure the root grid to allow resizing
        root.grid_rowconfigure(0, weight=0)  # Row 0 for the label
        root.grid_rowconfigure(1, weight=0)  # Row 1 for the button
        root.grid_rowconfigure(2, weight=0)  # Row 2 for the button
        root.grid_rowconfigure(3, weight=1)  # Row 3 for the scrollable frame

        root.grid_columnconfigure(0, weight=1)  # Column 0 for the label and buttons
        root.grid_columnconfigure(1, weight=5)  # Column 1 for the entry field

        # Label for the input field
        self.label = tk.Label(root, text="Enter CIDs (comma-separated):")
        self.label.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        # Entry field to input CIDs
        self.cid_entry = tk.Entry(root)
        self.cid_entry.grid(row=0, column=1, padx=10, pady=10, sticky="ew")

        # Fetch data button
        self.fetch_button = tk.Button(root, text="Fetch Compound Data", command=self.fetch_data)
        self.fetch_button.grid(row=1, column=0, columnspan=2, padx=10, pady=10)

        # Save to Excel button
        self.save_button = tk.Button(root, text="Save to Excel", command=self.save_to_excel)
        self.save_button.grid(row=2, column=0, columnspan=2, padx=10, pady=10)

        # Scrollable frame to hold the text box and scrollbar
        self.scroll_frame = tk.Frame(root)
        self.scroll_frame.grid(row=3, column=0, columnspan=2, padx=10, pady=10, sticky="nsew")

        self.text_box = tk.Text(self.scroll_frame, width=80, height=15, wrap=tk.WORD)
        self.text_box.grid(row=0, column=0, sticky="nsew")

        self.scrollbar = tk.Scrollbar(self.scroll_frame, command=self.text_box.yview)
        self.scrollbar.grid(row=0, column=1, sticky="ns")

        self.text_box.config(yscrollcommand=self.scrollbar.set)

        # Allow scroll_frame to resize with the window
        self.scroll_frame.grid_rowconfigure(0, weight=1)
        self.scroll_frame.grid_columnconfigure(0, weight=1)
        
    def fetch_data(self):
        logger.info("Fetching compound data...")
        # Get the CIDs from the entry field and process them
        cids = self.cid_entry.get().split(',')
        self.compound_data = []

        # Clear the text box before displaying new data
        self.text_box.delete(1.0, tk.END)

        for cid in cids:
            cid = cid.strip()  # Clean up the CID input
            logger.info(f"Fetching data for CID {cid}")
            if cid.isdigit():  # Check if CID is numeric
                compound_info = self.fetch_compound_by_cid(int(cid))
                logger.info(f"Compound data for CID {cid}: {compound_info}")
                if compound_info:
                    self.compound_data.extend(compound_info)

        if self.compound_data:
            self.display_data()
            logger.info("Compound data fetched successfully.")
            messagebox.showinfo("Success", "Compound data fetched successfully.")
        else:
            messagebox.showerror("Error", "No compound data found.")

    def fetch_compound_by_cid(self, cid):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
        response = requests.get(url)
        
        if response.status_code != 200:
            print(f"Failed to fetch data for CID {cid}")
            return None

        # Create a temporary file to save the SDF data
        with tempfile.NamedTemporaryFile(delete=False, suffix='.sdf') as tmpfile:
            tmpfile.write(response.content)
            tmpfile_path = tmpfile.name
        
        # Use SDMolSupplier with the temporary file path
        suppl = Chem.SDMolSupplier(tmpfile_path)
        
        compound_info = []
        for mol in suppl:
            if mol is None:
                continue  # Skip invalid molecules

            data = {
                "CID": cid,
                "OpenEye Name": mol.GetProp("PUBCHEM_IUPAC_OPENEYE_NAME") if mol.HasProp("PUBCHEM_IUPAC_OPENEYE_NAME") else "N/A",
                "InchI": Chem.MolToInchi(mol) if mol else "N/A",
                "InchI Key": mol.GetProp("PUBCHEM_IUPAC_INCHIKEY") if mol.HasProp("PUBCHEM_IUPAC_INCHIKEY") else "N/A",
                # "SMILES": Chem.MolToSmiles(mol) if mol else "N/A",
                "SMILES": mol.GetProp("PUBCHEM_SMILES") if mol.HasProp("PUBCHEM_SMILES") else "N/A",
                "OpenEye Can Smiles": mol.GetProp("PUBCHEM_OPENEYE_CAN_SMILES") if mol.HasProp("PUBCHEM_OPENEYE_CAN_SMILES") else "N/A",
                "OpenEye Iso Smiles": mol.GetProp("PUBCHEM_OPENEYE_ISO_SMILES") if mol.HasProp("PUBCHEM_OPENEYE_ISO_SMILES") else "N/A",
                "Molecular Formula": mol.GetProp("PUBCHEM_MOLECULAR_FORMULA") if mol.HasProp("PUBCHEM_MOLECULAR_FORMULA") else "N/A",
                # "CID": mol.GetProp("PUBCHEM_COMPOUND_CID") if mol.HasProp("PUBCHEM_COMPOUND_CID") else "N/A",
                "Canonicalized": mol.GetProp("PUBCHEM_COMPOUND_CANONICALIZED") if mol.HasProp("PUBCHEM_COMPOUND_CANONICALIZED") else "N/A",
                "Complexity": mol.GetProp("PUBCHEM_CACTVS_COMPLEXITY") if mol.HasProp("PUBCHEM_CACTVS_COMPLEXITY") else "N/A",
                "H-bond Acceptors": mol.GetProp("PUBCHEM_CACTVS_HBOND_ACCEPTOR") if mol.HasProp("PUBCHEM_CACTVS_HBOND_ACCEPTOR") else "N/A",
                "H-bond Donors": mol.GetProp("PUBCHEM_CACTVS_HBOND_DONOR") if mol.HasProp("PUBCHEM_CACTVS_HBOND_DONOR") else "N/A",
                "Rotatable Bonds": mol.GetProp("PUBCHEM_CACTVS_ROTATABLE_BOND") if mol.HasProp("PUBCHEM_CACTVS_ROTATABLE_BOND") else "N/A",
                # "Substructure Keys": mol.GetProp("PUBCHEM_CACTVS_SUBSKEYS") if mol.HasProp("PUBCHEM_CACTVS_SUBSKEYS") else "N/A",
                "Iupac Cas Name": mol.GetProp("PUBCHEM_IUPAC_CAS_NAME") if mol.HasProp("PUBCHEM_IUPAC_CAS_NAME") else "N/A",
                "Iupac Name Markup": mol.GetProp("PUBCHEM_IUPAC_NAME_MARKUP") if mol.HasProp("PUBCHEM_IUPAC_NAME_MARKUP") else "N/A",
                "Iupac Name": mol.GetProp("PUBCHEM_IUPAC_NAME") if mol.HasProp("PUBCHEM_IUPAC_NAME") else "N/A",
                "Iupac Systematic Name": mol.GetProp("PUBCHEM_IUPAC_SYSTEMATIC_NAME") if mol.HasProp("PUBCHEM_IUPAC_SYSTEMATIC_NAME") else "N/A",
                "Iupac Traditional Name": mol.GetProp("PUBCHEM_IUPAC_TRADITIONAL_NAME") if mol.HasProp("PUBCHEM_IUPAC_TRADITIONAL_NAME") else "N/A",
                "XLogP3": mol.GetProp("PUBCHEM_XLOGP3") if mol.HasProp("PUBCHEM_XLOGP3") else "N/A",
                "Exact Mass": mol.GetProp("PUBCHEM_EXACT_MASS") if mol.HasProp("PUBCHEM_EXACT_MASS") else "N/A",
                "Molecular Weight": mol.GetProp("PUBCHEM_MOLECULAR_WEIGHT") if mol.HasProp("PUBCHEM_MOLECULAR_WEIGHT") else "N/A",    
                "Monoisotopic Weight": mol.GetProp("PUBCHEM_MONOISOTOPIC_WEIGHT") if mol.HasProp("PUBCHEM_MONOISOTOPIC_WEIGHT") else "N/A",    
                "Total Charge": mol.GetProp("PUBCHEM_TOTAL_CHARGE") if mol.HasProp("PUBCHEM_TOTAL_CHARGE") else "N/A",
                "Heavy Atom Count": mol.GetProp("PUBCHEM_HEAVY_ATOM_COUNT") if mol.HasProp("PUBCHEM_HEAVY_ATOM_COUNT") else "N/A",
                "TPSA": mol.GetProp("PUBCHEM_CACTVS_TPSA") if mol.HasProp("PUBCHEM_CACTVS_TPSA") else "N/A",
                "Atom Def Stereo Count": mol.GetProp("PUBCHEM_ATOM_DEF_STEREO_COUNT") if mol.HasProp("PUBCHEM_ATOM_DEF_STEREO_COUNT") else "N/A",
                "Atom Udef Stereo Count": mol.GetProp("PUBCHEM_ATOM_UDEF_STEREO_COUNT") if mol.HasProp("PUBCHEM_ATOM_UDEF_STEREO_COUNT") else "N/A",
                "Bond Def Stereo Count": mol.GetProp("PUBCHEM_BOND_DEF_STEREO_COUNT") if mol.HasProp("PUBCHEM_BOND_DEF_STEREO_COUNT") else "N/A",
                "Bond Udef Stereo Count": mol.GetProp("PUBCHEM_BOND_UDEF_STEREO_COUNT") if mol.HasProp("PUBCHEM_BOND_UDEF_STEREO_COUNT") else "N/A",
                "Isotopic Atom Count": mol.GetProp("PUBCHEM_ISOTOPIC_ATOM_COUNT") if mol.HasProp("PUBCHEM_ISOTOPIC_ATOM_COUNT") else "N/A"
            }
            compound_info.append(data)
            # print(data)
        return compound_info


    def display_data(self):
        # Format the fetched compound data and insert into the text box
        display_text = ""
        for data in self.compound_data:
            display_text += f"CID: {data['CID']}\n"
            display_text += f"OpenEye Name: {data['OpenEye Name']}\n"
            display_text += f"InchI: {data['InchI']}\n"
            display_text += f"InChI Key: {data.get('InchI Key', 'N/A')}\n"
            display_text += f"SMILES: {data.get('SMILES', 'N/A')}\n" 
            display_text += f"OpenEye Can Smiles: {data['OpenEye Can Smiles']}\n"
            display_text += f"OpenEye Iso Smiles: {data['OpenEye Iso Smiles']}\n"
            display_text += f"Molecular Formula: {data['Molecular Formula']}\n"
            display_text += f"Canonicalized: {data['Canonicalized']}\n"
            display_text += f"Complexity: {data['Complexity']}\n"
            display_text += f"H-bond Acceptors: {data['H-bond Acceptors']}\n"
            display_text += f"H-bond Donors: {data['H-bond Donors']}\n"
            display_text += f"Rotatable Bonds: {data['Rotatable Bonds']}\n"            
            display_text += f"IUPAC Cas Name: {data.get('Iupac Cas Name', 'N/A')}\n"
            display_text += f"IUPAC Name Markup: {data.get('Iupac Name Markup', 'N/A')}\n"
            display_text += f"IUPAC Name: {data.get('Iupac Name', 'N/A')}\n"
            display_text += f"IUPAC Systematic Name: {data.get('Iupac Systematic Name', 'N/A')}\n"
            display_text += f"IUPAC Traditional Name: {data.get('Iupac Traditional Name', 'N/A')}\n"
            display_text += f"XLogP3: {data['XLogP3']}\n"
            display_text += f"Exact Mass: {data['Exact Mass']}\n"
            display_text += f"Molecular Weight: {data['Molecular Weight']}\n"
            display_text += f"Monoisotopic Weight: {data['Monoisotopic Weight']}\n"
            display_text += f"Total Charge: {data['Total Charge']}\n"
            display_text += f"Heavy Atom Count: {data['Heavy Atom Count']}\n"
            display_text += f"TPSA: {data['TPSA']}\n"
            display_text += f"Atom Def Stereo Count: {data['Atom Def Stereo Count']}\n"
            display_text += f"Atom Udef Stereo Count: {data['Atom Udef Stereo Count']}\n"
            display_text += f"Bond Def Stereo Count: {data['Bond Def Stereo Count']}\n"
            display_text += f"Bond Udef Stereo Count: {data['Bond Udef Stereo Count']}\n"
            display_text += f"Isotopic Atom Count: {data['Isotopic Atom Count']}\n"
            display_text += "-"*40 + "\n"  # Separator         

        self.text_box.insert(tk.END, display_text)

    def save_to_excel(self):
        if not self.compound_data:
            messagebox.showerror("Error", "No data to save.")
            return
        
        df = pd.DataFrame(self.compound_data)
        try:
            df.to_excel("compound_data.xlsx", index=False)
            messagebox.showinfo("Success", "Compound data has been saved to compound_data.xlsx")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save data: {e}")

# Set up the Tkinter root window
root = tk.Tk()
app = CompoundFetcherApp(root)

# Run the application
root.mainloop()
