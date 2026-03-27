import pandas as pd
import re
import scipy.io
import numpy as np

# ========== 1. Load iML1515 model ==========
mat = scipy.io.loadmat('../models/eciML1515-basic.mat')
model = mat['eciML1515_basic'].flatten()[0]

mets = np.array([m[0] if hasattr(m, '__len__') else m for m in model['mets'].flatten()])

# ========== 2. Define mapping from simplified names to model metabolites ==========
metabolite_mapping = {
    '*ATP': 'atp_c',
    'GLC': 'glc__D_c',
    'ATP': 'atp_c',
    'ADP': 'adp_c',
    'PI': 'pi_c',
    'H2O': 'h2o_c',
    'CO2': 'co2_c',
    'NADH': 'nadh_c',
    'NADPH': 'nadph_c',
    'NADP': 'nadp_c',
    'NAD': 'nad_c',
    'NH4': 'nh4_c',
    'NH3': 'nh3_c',
    'H': 'h_c',
    'H+': 'h_c',
    'PGA': '3pg_c',
    'PG': '3pg_c',
    'PEP': 'pep_c',
    'PYR': 'pyr_c',
    'AcCoA': 'accoa_c',
    'Acetate': 'ac_c',
    'G6P': 'g6p_c',
    'F6P': 'f6p_c',
    'FUM': 'fum_c',
    'MAL': 'mal__L_c',
    'OAA': 'oaa_c',
    'ICT': 'akg_c',
    'OGA': 'akg_c',
    'S7P': 's7p_c',
    'E4P': 'e4p_c',
    'P5P': 'r5p_c',
    'T2P': 'dhap_c',
    'T3P': 'g3p_c',
}

# ========== 3. Read the Excel file ==========
df = pd.read_excel('13CFluxdata.xlsx', header=None)

# ========== 4. Build output data ==========
output_rows = []

# Get number of conditions (columns 2 onwards)
n_conditions = 17  # Cond1 to Cond17

# --- Process constraints (rows 2-4) ---
for i in range(2, 5):
    rxn = df.iloc[i, 0]
    rxn_name = df.iloc[i, 1]
    if pd.notna(rxn):
        fluxes = df.iloc[i, 2:2+n_conditions].values
        # Filter out NaN
        valid_fluxes = [f if pd.notna(f) else '' for f in fluxes]

        output_rows.append({
            'Type': 'constraint',
            'RxnName': rxn_name,
            'ModelMetabolite': '{' + rxn_name + '}',
            'Coefficient': '{-1}',
            'Cond1': valid_fluxes[0],
            'Cond2': valid_fluxes[1],
            'Cond3': valid_fluxes[2],
            'Cond4': valid_fluxes[3],
            'Cond5': valid_fluxes[4],
            'Cond6': valid_fluxes[5],
            'Cond7': valid_fluxes[6],
            'Cond8': valid_fluxes[7],
            'Cond9': valid_fluxes[8],
            'Cond10': valid_fluxes[9],
            'Cond11': valid_fluxes[10],
            'Cond12': valid_fluxes[11],
            'Cond13': valid_fluxes[12],
            'Cond14': valid_fluxes[13],
            'Cond15': valid_fluxes[14],
            'Cond16': valid_fluxes[15],
            'Cond17': valid_fluxes[16],
        })

# --- Process flux reactions (row 8 onwards) ---
for i in range(8, len(df)):
    eq = df.iloc[i, 0]
    name = df.iloc[i, 1]
    if pd.notna(eq) and pd.notna(name):
        fluxes = df.iloc[i, 2:2+n_conditions].values

        # Parse equation to extract metabolites and coefficients
        parts = eq.split('->')
        if len(parts) != 2:
            continue

        reactants = [r.strip() for r in parts[0].split('+')]
        products = [p.strip() for p in parts[1].split('+')]

        # Build metabolite list and coefficient list
        metabolites = []
        coefficients = []

        # Reactants (positive coefficient)
        for r in reactants:
            match = re.match(r'^(\d+\.?\d*)\s*(.+)$', r)
            if match:
                coef = float(match.group(1))
                met = match.group(2)
            else:
                coef = 1.0
                met = r
            model_met = metabolite_mapping.get(met, met)
            metabolites.append(model_met)
            coefficients.append(int(coef) if coef == int(coef) else coef)

        # Products (negative coefficient)
        for p in products:
            match = re.match(r'^(\d+\.?\d*)\s*(.+)$', p)
            if match:
                coef = float(match.group(1))
                met = match.group(2)
            else:
                coef = 1.0
                met = p
            model_met = metabolite_mapping.get(met, met)
            metabolites.append(model_met)
            coefficients.append(-int(coef) if coef == int(coef) else -coef)

        # Format as MATLAB-style cell arrays
        met_str = '{' + ';'.join(metabolites) + '}'
        coef_str = '{' + ';'.join([str(c) for c in coefficients]) + '}'

        # Filter fluxes
        valid_fluxes = [f if pd.notna(f) else '' for f in fluxes]

        output_rows.append({
            'Type': 'flux',
            'RxnName': name,
            'ModelMetabolite': met_str,
            'Coefficient': coef_str,
            'Cond1': valid_fluxes[0],
            'Cond2': valid_fluxes[1],
            'Cond3': valid_fluxes[2],
            'Cond4': valid_fluxes[3],
            'Cond5': valid_fluxes[4],
            'Cond6': valid_fluxes[5],
            'Cond7': valid_fluxes[6],
            'Cond8': valid_fluxes[7],
            'Cond9': valid_fluxes[8],
            'Cond10': valid_fluxes[9],
            'Cond11': valid_fluxes[10],
            'Cond12': valid_fluxes[11],
            'Cond13': valid_fluxes[12],
            'Cond14': valid_fluxes[13],
            'Cond15': valid_fluxes[14],
            'Cond16': valid_fluxes[15],
            'Cond17': valid_fluxes[16],
        })

# ========== 5. Save to TSV ==========
output_df = pd.DataFrame(output_rows)
output_df.to_csv('13CFluxdata.tsv', sep='\t', index=False)

print(f"Generated TSV file with {len(output_df)} rows")
print(f"\nColumns: {list(output_df.columns)}")
print(f"\n=== Sample constraint rows ===")
print(output_df[output_df['Type'] == 'constraint'].head(3).to_string())
print(f"\n=== Sample flux rows ===")
print(output_df[output_df['Type'] == 'flux'].head(5).to_string())
