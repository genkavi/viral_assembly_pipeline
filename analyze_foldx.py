#!/usr/bin/env python3
"""
Analyze FoldX output files (Stability and BuildModel)
Extract energy terms and calculate viral envelope totals
"""

import argparse
import pandas as pd
import glob
from pathlib import Path
import sys
import numpy as np

# FoldX Stability column definitions
FOLDX_STABILITY_COLUMNS = [
    'pdb',
    'total_energy',
    'backbone_hbond',
    'sidechain_hbond',
    'van_der_waals',
    'electrostatics',
    'solvation_polar',
    'solvation_hydrophobic',
    'vdw_clashes',
    'entropy_sidechain',
    'entropy_mainchain',
    'sloop_entropy',
    'mloop_entropy',
    'cis_bond',
    'torsional_clash',
    'backbone_clash',
    'helix_dipole',
    'water_bridge',
    'disulfide',
    'electrostatic_kon',
    'partial_covalent',
    'energy_ionisation',
    'entropy_complex',
    'residue_number'
]

# Assembly type mappings - updated for your directory names
ASSEMBLY_TYPE_MAP = {
    '3-fold': 'threefold',
    '5-fold': 'fivefold',
    'dimer_para': 'side_to_side',
    'dimer_perp': 'side_to_end',
    'dimer': 'dimer'
}

# Icosahedral capsid architecture frequencies
CAPSID_FREQUENCIES = {
    'fivefold': 12,
    'threefold': 20,
    'side_to_side': 60,
    'side_to_end': 60,
    'dimer': 90
}

def parse_fxout_file(fxout_file, file_type='stability', debug=False):
    """
    Parse a single FoldX .fxout file
    
    Args:
        fxout_file: Path to .fxout file
        file_type: 'stability' or 'buildmodel' (determines column names)
        debug: Show debug output
    
    Returns:
        pd.DataFrame with all energy terms
    """
    
    try:
        # Read file
        with open(fxout_file, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            print(f"Warning: File is empty: {fxout_file}")
            return None
        
        if debug:
            print(f"\nDEBUG {Path(fxout_file).name}:")
            print(f"  Total lines: {len(lines)}")
            print(f"  First line: {repr(lines[0][:100])}")
        
        # Find data lines
        data_lines = []
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            if not line or line.startswith('#'):
                continue
            
            # Skip header line if present
            if 'Pdb' in line or 'Total Energy' in line:
                if debug:
                    print(f"  Skipping header: {line[:50]}")
                continue
            
            # Split the line - handle both tab and space delimited
            if '\t' in line:
                parts = line.split('\t')
            else:
                # Use whitespace split
                parts = line.split()
            
            # Fix concatenated numbers like "0-137.807"
            fixed_parts = []
            for part in parts:
                # Check if this part has a digit followed immediately by minus
                # e.g., "0-137.807" should become ["0", "-137.807"]
                import re
                if re.search(r'\d-', part):
                    # Split on the minus that follows a digit
                    subparts = re.split(r'(?<=\d)(?=-)', part)
                    fixed_parts.extend(subparts)
                else:
                    fixed_parts.append(part)
            
            parts = [p for p in fixed_parts if p]  # Remove empty
            
            if debug and i == 0:
                print(f"  Parts count: {len(parts)}")
                print(f"  First few parts: {parts[:5]}")
            
            # Validate this is a data line (has filename and numbers)
            if len(parts) >= 24:  # Expecting 24 columns for stability
                data_lines.append(parts)
            elif len(parts) >= 2:
                # At least filename + some data
                data_lines.append(parts)
        
        if debug:
            print(f"  Data lines found: {len(data_lines)}")
        
        if not data_lines:
            print(f"Warning: No data found in {fxout_file}")
            # Show first non-empty line for debugging
            for line in lines:
                if line.strip() and not line.startswith('#'):
                    print(f"  First data line attempt: {line[:100]}")
                    break
            return None
        
        # Use predefined columns
        columns = FOLDX_STABILITY_COLUMNS
        
        # Make sure we have enough columns
        max_cols = max(len(row) for row in data_lines)
        if len(columns) < max_cols:
            columns = columns + [f'extra_col_{i}' for i in range(len(columns), max_cols)]
        
        # Pad shorter rows
        for row in data_lines:
            while len(row) < max_cols:
                row.append(None)
        
        # Create DataFrame
        df = pd.DataFrame(data_lines, columns=columns[:max_cols])
        
        # Convert numeric columns
        for col in df.columns[1:]:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Add metadata
        df['source_file'] = Path(fxout_file).name
        df['assembly_type_raw'] = Path(fxout_file).parent.name
        
        # Map to standard assembly type
        assembly_raw = df['assembly_type_raw'].iloc[0].lower()
        df['assembly_type'] = ASSEMBLY_TYPE_MAP.get(assembly_raw, assembly_raw)
        
        # Extract frame
        filename = Path(fxout_file).stem.replace('_ST', '').replace('Stability_', '')
        df['frame'] = filename
        
        return df
        
    except Exception as e:
        print(f"Error parsing {fxout_file}: {e}")
        import traceback
        traceback.print_exc()
        return None

def parse_all_fxout_files(input_dir, file_type='stability'):
    """Parse all .fxout files in a directory"""
    
    fxout_files = glob.glob(f"{input_dir}/**/*.fxout", recursive=True)
    
    if not fxout_files:
        print(f"No .fxout files found in {input_dir}")
        return None
    
    print(f"Found {len(fxout_files)} .fxout files")
    
    # Parse all files
    all_data = []
    for idx, fxout_file in enumerate(fxout_files):
        # Debug first 3 failing files
        debug = idx < 3
        df = parse_fxout_file(fxout_file, file_type, debug=debug)
        if df is not None:
            all_data.append(df)
    
    if not all_data:
        print("No valid data extracted")
        return None
    
    combined_df = pd.concat(all_data, ignore_index=True)
    print(f"Extracted {len(combined_df)} data rows")
    
    return combined_df


def calculate_assembly_totals(df, energy_column='total_energy'):
    """
    Calculate total assembly energy by summing across all chains for each frame
    
    Args:
        df: DataFrame with parsed FoldX data
        energy_column: Which energy column to sum (default: 'total_energy')
    
    Returns:
        pd.DataFrame with assembly totals per frame and assembly type
    """
    
    if energy_column not in df.columns:
        print(f"Warning: Column '{energy_column}' not found. Available columns: {df.columns.tolist()}")
        return None
    
    # Group by assembly type, frame, and source file
    # Sum the energy across all chains in each structure
    assembly_totals = df.groupby(['assembly_type', 'frame', 'source_file']).agg({
        energy_column: 'sum',
        'pdb': 'count'  # Count number of chains
    }).reset_index()
    
    assembly_totals.rename(columns={
        energy_column: f'{energy_column}_sum',
        'pdb': 'num_chains'
    }, inplace=True)
    
    return assembly_totals


def calculate_viral_envelope_energy(assembly_df, energy_column='total_energy_sum'):
    """
    Calculate complete viral envelope energy following the paper's formula:
    
    ΔΔG_viral_envelope = 12×vertex_5fold + 20×vertex_3fold + 60×inter_dimer_parallel + 60×inter_dimer_perp + 90×intra_dimer
    
    Where interface-specific energies remove embedded dimer contributions:
    vertex_5fold = fivefold_raw - 5×dimer
    vertex_3fold = threefold_raw - 3×dimer
    inter_dimer_parallel = side_to_side_raw - 2×dimer
    inter_dimer_perp = side_to_end_raw - 2×dimer
    intra_dimer = dimer
    
    Args:
        assembly_df: DataFrame with assembly totals (from calculate_assembly_totals)
        energy_column: Which energy column to use (default: 'total_energy_sum')
    
    Returns:
        pd.DataFrame with viral envelope energies per frame
    """
    
    if energy_column not in assembly_df.columns:
        print(f"Error: Column '{energy_column}' not found in assembly data")
        return None
    
    # DEBUG: Show what assembly types we have
    print(f"\nAssembly types in data: {assembly_df['assembly_type'].unique()}")
    
    # Pivot to get one row per frame with columns for each assembly type
    pivot_df = assembly_df.pivot_table(
        index='frame',
        columns='assembly_type',
        values=energy_column,
        aggfunc='mean'
    ).reset_index()
    
    print(f"Columns after pivot: {pivot_df.columns.tolist()}")
    
    # Check which assembly types are present
    required_types = ['fivefold', 'threefold', 'side_to_side', 'side_to_end', 'dimer']
    missing_types = [t for t in required_types if t not in pivot_df.columns]
    
    if missing_types:
        print(f"\n⚠ Warning: Missing required assembly types: {missing_types}")
        print(f"Available types: {[c for c in pivot_df.columns if c != 'frame']}")
        return None
    
    # Rename to clarify these are raw values (before correction)
    pivot_df.rename(columns={
        'fivefold': 'fivefold_raw',
        'threefold': 'threefold_raw',
        'side_to_side': 'side_to_side_raw',
        'side_to_end': 'side_to_end_raw',
        'dimer': 'intra_dimer'
    }, inplace=True)
    
    # Calculate interface-specific contributions (remove embedded dimer contributions)
    pivot_df['vertex_5fold'] = pivot_df['fivefold_raw'] - 5 * pivot_df['intra_dimer']
    pivot_df['vertex_3fold'] = pivot_df['threefold_raw'] - 3 * pivot_df['intra_dimer']
    pivot_df['inter_dimer_parallel'] = pivot_df['side_to_side_raw'] - 2 * pivot_df['intra_dimer']
    pivot_df['inter_dimer_perp'] = pivot_df['side_to_end_raw'] - 2 * pivot_df['intra_dimer']
    
    # Calculate total viral envelope energy with icosahedral weights
    pivot_df['viral_envelope_energy'] = (
        12 * pivot_df['vertex_5fold'] +           # 12 five-fold vertices
        20 * pivot_df['vertex_3fold'] +           # 20 three-fold vertices
        60 * pivot_df['inter_dimer_parallel'] +   # 60 parallel dimer-dimer interfaces
        60 * pivot_df['inter_dimer_perp'] +       # 60 perpendicular dimer-dimer interfaces
        90 * pivot_df['intra_dimer']              # 90 intra-dimer interfaces
    )
    
    # Select columns in logical order
    result_cols = [
        'frame',
        # Raw assembly energies
        'fivefold_raw', 'threefold_raw', 'side_to_side_raw', 'side_to_end_raw', 'intra_dimer',
        # Interface-specific contributions
        'vertex_5fold', 'vertex_3fold', 'inter_dimer_parallel', 'inter_dimer_perp',
        # Total
        'viral_envelope_energy'
    ]
    
    return pivot_df[result_cols]

def main():
    parser = argparse.ArgumentParser(
        description="Analyze FoldX output files and calculate viral envelope energies",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze Stability results
  python analyze_foldx.py --input stability_results/ --output stability_data.csv
  
  # Calculate viral envelope energies from stability
  python analyze_foldx.py --input stability_results/ \\
      --output stability_data.csv \\
      --assembly-totals assembly_totals.csv \\
      --viral-envelope viral_envelope.csv \\
      --energy-column total_energy
  
  # Analyze BuildModel mutation results (use ΔΔG)
  python analyze_foldx.py --input mutation_results/ \\
      --output mutation_data.csv \\
      --assembly-totals mutation_assembly.csv \\
      --viral-envelope mutation_envelope.csv \\
      --energy-column total_energy \\
      --file-type buildmodel
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input directory containing .fxout files')
    parser.add_argument('-o', '--output', required=True,
                       help='Output CSV file for all individual data')
    parser.add_argument('--assembly-totals', 
                       help='Output CSV file for assembly totals (sum across chains)')
    parser.add_argument('--viral-envelope',
                       help='Output CSV file for complete viral envelope energies')
    parser.add_argument('--energy-column', default='total_energy',
                       help='Energy column to use (default: total_energy for Stability, or specify ddG for mutations)')
    parser.add_argument('--file-type', default='stability', choices=['stability', 'buildmodel'],
                       help='Type of FoldX output files (default: stability)')
    
    args = parser.parse_args()
    
    # Parse all files
    print(f"Parsing FoldX {args.file_type} files from: {args.input}")
    df = parse_all_fxout_files(args.input, args.file_type)
    
    if df is None:
        print("Error: No data extracted")
        sys.exit(1)
    
    # Save all individual data
    df.to_csv(args.output, index=False)
    print(f"\n✓ Saved individual data to: {args.output}")
    print(f"  Total rows: {len(df)}")
    print(f"  Columns: {', '.join(df.columns[:10])}...")
    
    # Calculate assembly totals if requested
    if args.assembly_totals or args.viral_envelope:
        print(f"\nCalculating assembly totals using column: {args.energy_column}")
        assembly_df = calculate_assembly_totals(df, args.energy_column)
        
        if assembly_df is None:
            print("Error: Could not calculate assembly totals")
            sys.exit(1)
        
        if args.assembly_totals:
            assembly_df.to_csv(args.assembly_totals, index=False)
            print(f"\n✓ Saved assembly totals to: {args.assembly_totals}")
            print(f"  Total assemblies: {len(assembly_df)}")
        
        # Calculate viral envelope energies if requested
        if args.viral_envelope:
            print(f"\nCalculating viral envelope energies...")
            energy_col = f'{args.energy_column}_sum'
            envelope_df = calculate_viral_envelope_energy(assembly_df, energy_col)
            
            if envelope_df is None:
                print("Error: Could not calculate viral envelope energies")
                sys.exit(1)
            
            envelope_df.to_csv(args.viral_envelope, index=False)
            print(f"\n✓ Saved viral envelope energies to: {args.viral_envelope}")
            print(f"  Total frames: {len(envelope_df)}")
            
            # Show summary
            print("\nViral Envelope Energy Summary:")
            print(envelope_df[['frame', 'viral_envelope_energy']].describe())
    
    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)


if __name__ == "__main__":
    main()