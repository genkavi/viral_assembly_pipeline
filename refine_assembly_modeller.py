#!/usr/bin/env python3
"""
Refine assembled structures using Modeller optimization
Runs conjugate gradients + molecular dynamics to relax clashes
"""

import argparse
import glob
import sys
import shutil
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import multiprocessing as mp

try:
    from modeller import *
    from modeller.scripts import complete_pdb
    from modeller.optimizers import ConjugateGradients, MolecularDynamics, actions
    MODELLER_AVAILABLE = True
except ImportError:
    MODELLER_AVAILABLE = False
    print("WARNING: Modeller not available. Install with conda or from https://salilab.org/modeller/")


def refine_structure(input_pdb, output_pdb, 
                     cg_steps=50, md_steps=100, md_temp=300,
                     verbose=False):
    """Refine structure using Modeller optimization"""
    
    if not MODELLER_AVAILABLE:
        return False, None, None, "Modeller not available"
    
    # Setup Modeller environment
    env = Environ()
    env.io.atom_files_directory = ['.']
    env.edat.dynamic_sphere = True
    
    # Read topology and parameters
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    
    try:
        # Load structure
        mdl = complete_pdb(env, input_pdb)
        
        # Restore original residue numbering from input file
        original_mdl = Model(env)
        original_mdl.read(file=input_pdb)
        
        # Copy residue numbers back
        for i, (res, orig_res) in enumerate(zip(mdl.residues, original_mdl.residues)):
            res.num = orig_res.num
        
        # Select all atoms
        atmsel = Selection(mdl)
        
        # Generate restraints
        mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
        
        # Initial energy
        initial_energy = atmsel.energy()[0]
        
        # Optimization
        output_level = 'REPORT' if verbose else 'NO_REPORT'
        cg = ConjugateGradients(output=output_level)
        md = MolecularDynamics(output=output_level)
        
        cg.optimize(atmsel, max_iterations=cg_steps)
        md.optimize(atmsel, temperature=md_temp, max_iterations=md_steps, md_return='MINIMAL')
        cg.optimize(atmsel, max_iterations=cg_steps)
        
        # Final energy with breakdown
        edat = atmsel.energy()
        final_energy = edat[0]

        # Write PDB FIRST (before report, so you get output even if report fails)
        mdl.write(file=output_pdb)

        # Write report
        report_file = output_pdb.replace('.pdb', '.out')
        try:
            with open(report_file, 'w') as f:
                f.write(f"Modeller Refinement Report\n")
                f.write(f"=" * 50 + "\n\n")
                f.write(f"Input:  {input_pdb}\n")
                f.write(f"Output: {output_pdb}\n\n")
                f.write(f"Energy (kcal/mol):\n")
                f.write(f"  Initial:         {initial_energy:12.2f}\n")
                f.write(f"  Final:           {final_energy:12.2f}\n")
                f.write(f"  Change (ΔE):     {final_energy-initial_energy:12.2f}\n\n")             
                f.write(f"\nOptimization Protocol:\n")
                f.write(f"  CG steps:        {cg_steps}\n")
                f.write(f"  MD steps:        {md_steps}\n")
                f.write(f"  MD temperature:  {md_temp} K\n")
        except Exception as e:
            # At least we have the PDB
            print(f"Warning: Could not write full report for {output_pdb}: {e}")

        return True, initial_energy, final_energy, None
        
        # Write with original numbering
        mdl.write(file=output_pdb)
        
    except Exception as e:
        return False, None, None, str(e)


def process_single_file(task):
    """Process a single PDB file"""
    
    input_file = task['input_file']
    output_file = task['output_file']
    cg_steps = task['cg_steps']
    md_steps = task['md_steps']
    md_temp = task['md_temp']
    
    try:
        success, e_init, e_final, error = refine_structure(
            input_file, output_file,
            cg_steps=cg_steps,
            md_steps=md_steps,
            md_temp=md_temp,
            verbose=False
        )
        
        if success:
            delta_e = e_final - e_init if (e_init and e_final) else None
            return {
                'status': 'success',
                'input': input_file,
                'output': output_file,
                'energy_initial': e_init,
                'energy_final': e_final,
                'energy_change': delta_e
            }
        else:
            return {
                'status': 'failed',
                'input': input_file,
                'error': error
            }
            
    except Exception as e:
        return {
            'status': 'error',
            'input': input_file,
            'error': str(e)
        }


def get_input_files(input_paths):
    """Get list of PDB files from input paths"""
    
    all_files = []
    
    for path in input_paths:
        path_obj = Path(path)
        
        if path_obj.is_file():
            if path_obj.suffix == '.pdb':
                all_files.append(str(path_obj))
        elif path_obj.is_dir():
            all_files.extend(glob.glob(f"{path}/*.pdb"))
        else:
            # Treat as glob pattern
            all_files.extend(glob.glob(path))
    
    return sorted(set(all_files))


def main():
    parser = argparse.ArgumentParser(
        description="Refine assembled structures using Modeller optimization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Refine 3-fold assemblies
  python refine_assembly_modeller.py \\
      --input assembly/3-fold/*.pdb \\
      --output-dir assembly-refined/3-fold \\
      --workers 8
  
  # Quick refinement (fewer steps)
  python refine_assembly_modeller.py \\
      --input assembly/*/*.pdb \\
      --output-dir assembly-refined \\
      --cg-steps 20 --md-steps 50 \\
      --workers 16
  
  # Thorough refinement (more steps, higher temperature)
  python refine_assembly_modeller.py \\
      --input assembly/*.pdb \\
      --output-dir refined \\
      --cg-steps 100 --md-steps 200 --md-temp 400 \\
      --workers 8

Note: Modeller must be installed (conda install -c salilab modeller)
        """
    )
    
    parser.add_argument('-i', '--input', nargs='+', required=True,
                       help='Input PDB files (files, globs, or directories)')
    parser.add_argument('-o', '--output-dir', required=True,
                       help='Output directory')
    parser.add_argument('--suffix', default='',
                       help='Suffix to add to output filenames (before .pdb)')
    parser.add_argument('--cg-steps', type=int, default=50,
                       help='Conjugate gradient steps per round (default: 50)')
    parser.add_argument('--md-steps', type=int, default=100,
                       help='Molecular dynamics steps (default: 100)')
    parser.add_argument('--md-temp', type=int, default=300,
                       help='MD temperature in Kelvin (default: 300)')
    parser.add_argument('-w', '--workers', type=int,
                       help='Number of parallel workers (default: CPU count - 2)')
    
    args = parser.parse_args()
    
    # Check Modeller
    if not MODELLER_AVAILABLE:
        print("ERROR: Modeller is required but not installed")
        print("Install with: conda install -c salilab modeller")
        print("Or download from: https://salilab.org/modeller/")
        sys.exit(1)
    
    # Get input files
    input_files = get_input_files(args.input)
    
    if not input_files:
        print("ERROR: No PDB files found")
        sys.exit(1)
    
    print(f"Found {len(input_files)} PDB files")
    print(f"Optimization: CG({args.cg_steps}) → MD({args.md_steps}@{args.md_temp}K) → CG({args.cg_steps})")
    
    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Determine number of workers
    if args.workers is None:
        workers = max(1, mp.cpu_count() - 2)
    else:
        workers = args.workers
    
    print(f"Using {workers} parallel workers")
    
    # Create tasks
    tasks = []
    # Detect common input prefix to preserve structure
    input_paths = [Path(f) for f in input_files]
    common_parent = Path(os.path.commonpath([p.parent for p in input_paths]))

    for input_file in input_files:
        input_path = Path(input_file)
        
        # Preserve subdirectory structure
        rel_path = input_path.relative_to(common_parent)
        
        if args.suffix:
            output_name = f"{rel_path.stem}{args.suffix}.pdb"
        else:
            output_name = rel_path.name
        
        output_file = Path(args.output_dir) / rel_path.parent.name / output_name
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        tasks.append({
            'input_file': input_file,
            'output_file': str(output_file),
            'cg_steps': args.cg_steps,
            'md_steps': args.md_steps,
            'md_temp': args.md_temp
        })
    
    # Process in parallel
    results = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(process_single_file, task) for task in tasks]
        
        for future in tqdm(futures, desc="Refining structures"):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"Task failed: {e}")
                results.append({'status': 'error', 'error': str(e)})
    
    # Report results
    successful = sum(1 for r in results if r['status'] == 'success')
    failed = len(results) - successful
    
    print(f"\n✓ Refinement complete: {successful}/{len(results)} successful")
    
    if failed > 0:
        print(f"\n⚠  {failed} files failed:")
        for r in results:
            if r['status'] != 'success':
                print(f"  - {Path(r.get('input', 'unknown')).name}: {r.get('error', 'unknown error')}")
    
    # Energy statistics
    if successful > 0:
        successful_results = [r for r in results if r['status'] == 'success' and r.get('energy_change') is not None]
        
        if successful_results:
            energy_changes = [r['energy_change'] for r in successful_results]
            avg_change = sum(energy_changes) / len(energy_changes)
            
            print(f"\nEnergy change statistics:")
            print(f"  Average ΔE: {avg_change:+.2f}")
            print(f"  Min ΔE: {min(energy_changes):+.2f}")
            print(f"  Max ΔE: {max(energy_changes):+.2f}")
            
            improved = sum(1 for e in energy_changes if e < 0)
            print(f"  Structures improved: {improved}/{len(energy_changes)}")


if __name__ == "__main__":
    main()