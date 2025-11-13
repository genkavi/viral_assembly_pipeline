#!/usr/bin/env python3
"""
Create assemblies by superimposing input structures onto template arms
"""

import argparse
import glob
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore")

from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Superimposer


def get_ca_atoms(chain):
    """Get CA atoms from a chain"""
    ca_atoms = []
    for residue in chain:
        if 'CA' in residue:
            ca_atoms.append(residue['CA'])
    return ca_atoms


def copy_structure(original_structure, new_id):
    """Create a deep copy of structure"""
    new_structure = Structure.Structure(new_id)
    for model in original_structure:
        new_model = Model.Model(model.id)
        for chain in model:
            new_chain = Chain.Chain(chain.id)
            for residue in chain:
                new_chain.add(residue.copy())
            new_model.add(new_chain)
        new_structure.add(new_model)
    return new_structure


def assemble_structure(input_pdb, template_pdb, output_pdb, reference_chain=0):
    """
    Superimpose input structure onto each arm of template structure

    Args:
        input_pdb: Input structure (e.g., dimer)
        template_pdb: Template structure (e.g., 3-fold, 5-fold)
        output_pdb: Output assembly
        reference_chain: Which chain from input to use as reference (default: 0)
    """

    # Load structures
    parser = PDBParser(QUIET=True)
    input_struct = parser.get_structure("input", input_pdb)
    template_struct = parser.get_structure("template", template_pdb)

    # Get chains from both structures
    input_chains = []
    for model in input_struct:
        for chain in model:
            ca_atoms = get_ca_atoms(chain)
            if ca_atoms:
                input_chains.append((chain.id, chain, ca_atoms))

    template_chains = []
    for model in template_struct:
        for chain in model:
            ca_atoms = get_ca_atoms(chain)
            if ca_atoms:
                template_chains.append((chain.id, chain, ca_atoms))

    if not input_chains:
        raise ValueError(f"No chains with CA atoms found in {input_pdb}")
    if not template_chains:
        raise ValueError(f"No chains with CA atoms found in {template_pdb}")

    # Validate reference chain
    if reference_chain >= len(input_chains):
        reference_chain = 0

    # Get reference chain from input
    ref_chain_id, ref_chain, ref_ca_atoms = input_chains[reference_chain]

    # Create assembly
    assembly = Structure.Structure("assembly")
    assembly_model = Model.Model(0)
    assembly.add(assembly_model)

    # Superimpose input onto each template chain
    superimposer = Superimposer()

    for i, (target_chain_id, target_chain, target_ca_atoms) in enumerate(template_chains):
        # Make copy of entire input structure
        input_copy = copy_structure(input_struct, f"copy_{i}")

        # Get reference CA atoms from copy
        copy_ref_ca = None
        for model in input_copy:
            for chain in model:
                if chain.id == ref_chain_id:
                    copy_ref_ca = get_ca_atoms(chain)
                    break

        if copy_ref_ca is None:
            continue

        # Perform superposition
        try:
            # Align copy reference chain to target chain
            min_len = min(len(copy_ref_ca), len(target_ca_atoms))
            superimposer.set_atoms(target_ca_atoms[:min_len], copy_ref_ca[:min_len])

            # Apply transformation to ALL atoms in input copy
            atom_list = []
            for model in input_copy:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atom_list.append(atom)

            superimposer.apply(atom_list)

            # Add transformed chains to assembly
            chain_counter = 0
            for model in input_copy:
                for chain in model:
                    # Give unique chain ID
                    new_chain_id = chr(ord('A') + i * len(input_chains) + chain_counter)
                    chain.id = new_chain_id
                    assembly_model.add(chain)
                    chain_counter += 1

        except Exception as e:
            raise RuntimeError(f"Superposition failed for chain {target_chain_id}: {e}")

    # Save result
    io = PDBIO()
    io.set_structure(assembly)
    io.save(output_pdb)

    total_chains = len(list(assembly.get_chains()))

    return {
        'input_chains': len(input_chains),
        'template_chains': len(template_chains),
        'output_chains': total_chains
    }


def process_single_assembly(task):
    """Process a single assembly"""

    input_file = task['input_file']
    template_file = task['template_file']
    output_file = task['output_file']
    reference_chain = task['reference_chain']

    try:
        info = assemble_structure(
            input_file,
            template_file,
            output_file,
            reference_chain
        )

        return {
            'status': 'success',
            'input': input_file,
            'output': output_file,
            'info': info
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
        description="Create assemblies by superimposing structures onto template arms",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create 3-fold assemblies
  python assemble_structure.py --input frames_opt/*.pdb --template 3-fold.pdb --output-dir 3fold_asm/

  # Create 5-fold assemblies with suffix
  python assemble_structure.py --input frames_opt/ --template 5-fold.pdb --output-dir 5fold_asm/ --suffix _5fold

  # Parallel processing
  python assemble_structure.py --input frames_opt/*.pdb --template 3-fold.pdb --output-dir 3fold/ --workers 8
        """
    )

    parser.add_argument('-i', '--input', nargs='+', required=True,
                       help='Input PDB files (files, globs, or directories)')
    parser.add_argument('-t', '--template', required=True,
                       help='Template PDB file (e.g., 3-fold.pdb, 5-fold.pdb)')
    parser.add_argument('-o', '--output-dir', required=True,
                       help='Output directory')
    parser.add_argument('--suffix', default='',
                       help='Suffix to add to output filenames (before .pdb)')
    parser.add_argument('-r', '--reference-chain', type=int, default=0,
                       help='Reference chain index from input structure (default: 0)')
    parser.add_argument('-w', '--workers', type=int,
                       help='Number of parallel workers (default: CPU count - 2)')

    args = parser.parse_args()

    # Check template exists
    if not Path(args.template).exists():
        print(f"ERROR: Template file not found: {args.template}")
        sys.exit(1)

    # Get input files
    input_files = get_input_files(args.input)

    if not input_files:
        print("ERROR: No PDB files found")
        sys.exit(1)

    print(f"Found {len(input_files)} input PDB files")
    print(f"Template: {args.template}")

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
    for input_file in input_files:
        base_name = Path(input_file).stem
        if args.suffix:
            output_name = f"{base_name}{args.suffix}.pdb"
        else:
            output_name = f"{base_name}.pdb"

        output_file = Path(args.output_dir) / output_name

        tasks.append({
            'input_file': input_file,
            'template_file': args.template,
            'output_file': str(output_file),
            'reference_chain': args.reference_chain
        })

    # Process in parallel
    results = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(process_single_assembly, task) for task in tasks]

        for future in tqdm(futures, desc="Creating assemblies"):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"Task failed: {e}")
                results.append({'status': 'error', 'error': str(e)})

    # Report results
    successful = sum(1 for r in results if r['status'] == 'success')
    failed = len(results) - successful

    print(f"\n✓ Assembly complete: {successful}/{len(results)} successful")

    if failed > 0:
        print(f"\n⚠ {failed} assemblies failed:")
        for r in results:
            if r['status'] != 'success':
                print(f"  - {Path(r.get('input', 'unknown')).name}: {r.get('error', 'unknown error')}")


if __name__ == "__main__":
    main()
