#!/usr/bin/env python3
"""
Wrapper for FoldX commands (mutations, stability analysis, etc.)
"""

import argparse
import glob
import os
import subprocess
import sys
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import multiprocessing as mp

# Import global configuration
try:
    import config
    DEFAULT_FOLDX_PATH = config.FOLDX_PATH
    DEFAULT_ROTABASE_PATH = config.ROTABASE_PATH
    DEFAULT_WORKERS = config.DEFAULT_WORKERS
except ImportError:
    # Fallback if config.py not found
    DEFAULT_FOLDX_PATH = os.path.expanduser("~/foldx_4.1/foldx_20251231")
    DEFAULT_ROTABASE_PATH = os.path.expanduser("~/foldx_4.1/rotabase.txt")
    DEFAULT_WORKERS = None


from Bio.PDB import PDBParser

def get_chains_in_pdb(pdb_file):
    """Extract actual chain IDs present in PDB"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('tmp', pdb_file)
    chains = set()
    for model in structure:
        for chain in model:
            chains.add(chain.id)
    return sorted(chains)


def run_foldx_buildmodel(input_pdb, output_dir, mutations, chains,
                         foldx_path, rotabase_path, combine_mutations=False):
    """Run FoldX BuildModel command
    
    Args:
        combine_mutations: If True, apply all mutations simultaneously in one structure
                          If False, create separate structure for each mutation (default)
    """

    # Detect actual chains present
    available_chains = get_chains_in_pdb(input_pdb)
    
    # Only use chains that exist
    valid_chains = [c for c in chains if c in available_chains]
    missing_chains = [c for c in chains if c not in available_chains]
    
    # Print chain info
    pdb_name = Path(input_pdb).name
    print(f"\n{pdb_name}:")
    print(f"  Detected chains: {', '.join(available_chains)}")
    print(f"  Using chains: {', '.join(valid_chains)}")
    if missing_chains:
        print(f"  ⚠️  Skipping missing chains: {', '.join(missing_chains)}")
    
    if not valid_chains:
        return False, f"None of the specified chains {chains} found in structure"

    # Create working directory
    work_dir = Path(output_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # Copy PDB file to working directory
    work_pdb = work_dir / pdb_name
    shutil.copy2(input_pdb, work_pdb)

    # Create mutation string for valid chains only
    mutation_strings = []
    
    if combine_mutations:
        # Apply ALL mutations simultaneously (combinatorial)
        all_chain_mutations = []
        for mutation in mutations:
            orig_aa = mutation[0]
            pos = mutation[1:-1]
            target_aa = mutation[-1]
            
            # Create mutation for each VALID chain
            for chain in valid_chains:
                all_chain_mutations.append(f"{orig_aa}{chain}{pos}{target_aa}")
        
        mutation_strings.append(",".join(all_chain_mutations) + ";")
    
    else:
        # Apply mutations separately
        for mutation in mutations:
            orig_aa = mutation[0]
            pos = mutation[1:-1]
            target_aa = mutation[-1]

            # Create mutation for each VALID chain
            chain_mutations = []
            for chain in valid_chains:
                chain_mutations.append(f"{orig_aa}{chain}{pos}{target_aa}")

            mutation_strings.append(",".join(chain_mutations) + ";")

    # Write mutation file
    mut_file = work_dir / "individual_list.txt"
    with open(mut_file, 'w') as f:
        for mut_str in mutation_strings:
            f.write(mut_str + "\n")

    # Run FoldX BuildModel
    cmd = [
        foldx_path,
        "--rotabaseLocation", rotabase_path,
        "--command=BuildModel",
        f"--pdb={pdb_name}",
        "--mutant-file=individual_list.txt",
        "--numberOfRuns=1"
    ]

    result = subprocess.run(
        cmd,
        cwd=work_dir,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        return False, result.stderr

    return True, None


def run_foldx_stability(input_pdb, output_dir, foldx_path, rotabase_path):
    """Run FoldX Stability command"""

    # Create working directory
    work_dir = Path(output_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # Copy PDB file to working directory
    pdb_name = Path(input_pdb).name
    work_pdb = work_dir / pdb_name
    shutil.copy2(input_pdb, work_pdb)

    # Run FoldX Stability
    cmd = [
        foldx_path,
        "--rotabaseLocation", rotabase_path,
        "--command=Stability",
        f"--pdb={pdb_name}"
    ]

    result = subprocess.run(
        cmd,
        cwd=work_dir,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        return False, result.stderr

    return True, None


def process_single_foldx_task(task):
    """Process a single FoldX task"""

    input_file = task['input_file']
    output_dir = task['output_dir']
    command = task['command']
    foldx_path = task['foldx_path']
    rotabase_path = task['rotabase_path']

    try:
        if command == 'BuildModel':
            mutations = task['mutations']
            chains = task['chains']
            combine_mutations = task.get('combine_mutations', False)
            success, error = run_foldx_buildmodel(
                input_file, output_dir, mutations, chains,
                foldx_path, rotabase_path, combine_mutations
            )
        elif command == 'Stability':
            success, error = run_foldx_stability(
                input_file, output_dir, foldx_path, rotabase_path
            )
        else:
            return {
                'status': 'error',
                'input': input_file,
                'error': f"Unknown command: {command}"
            }

        if success:
            return {
                'status': 'success',
                'input': input_file,
                'output_dir': output_dir
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

def run_foldx_alascan(input_pdb, output_dir, foldx_path, rotabase_path):
    """Run FoldX AlaScan command"""

    # Create working directory
    work_dir = Path(output_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # Copy PDB file to working directory
    pdb_name = Path(input_pdb).name
    work_pdb = work_dir / pdb_name
    shutil.copy2(input_pdb, work_pdb)

    # Run FoldX AlaScan - no chain specification needed
    cmd = [
        foldx_path,
        "--rotabaseLocation", rotabase_path,
        "--command=AlaScan",
        f"--pdb={pdb_name}"
    ]

    result = subprocess.run(
        cmd,
        cwd=work_dir,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        error_msg = f"Return code: {result.returncode}\n"
        if result.stderr:
            error_msg += f"STDERR: {result.stderr}\n"
        if result.stdout:
            error_msg += f"STDOUT: {result.stdout}"
        return False, error_msg

    return True, None

def process_single_foldx_task(task):
    """Process a single FoldX task"""

    input_file = task['input_file']
    output_dir = task['output_dir']
    command = task['command']
    foldx_path = task['foldx_path']
    rotabase_path = task['rotabase_path']

    try:
        if command == 'BuildModel':
            mutations = task['mutations']
            chains = task['chains']
            combine_mutations = task.get('combine_mutations', False)
            success, error = run_foldx_buildmodel(
                input_file, output_dir, mutations, chains,
                foldx_path, rotabase_path, combine_mutations
            )
        elif command == 'Stability':
            success, error = run_foldx_stability(
                input_file, output_dir, foldx_path, rotabase_path
            )
        elif command == 'AlaScan':
            success, error = run_foldx_alascan(
                input_file, output_dir,
                foldx_path, rotabase_path
            )
        else:
            return {
                'status': 'error',
                'input': input_file,
                'error': f"Unknown command: {command}"
            }

        if success:
            return {
                'status': 'success',
                'input': input_file,
                'output_dir': output_dir
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



def main():
    parser = argparse.ArgumentParser(
        description="Wrapper for FoldX commands",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run mutations on 3-fold assemblies
  python foldx_wrapper.py --command BuildModel \\
      --input 3fold_opt/*.pdb \\
      --mutations E327D E327A E327S \\
      --chains A B C D E F \\
      --output-dir results/mutations/ \\
      --workers 8
      
  python foldx_wrapper.py --command AlaScan \\
      --input assemblies/*.pdb \\
      --chains A B C D E F \\
      --output-dir results/alascan/ \\
      --workers 8

  # Run stability analysis
  python foldx_wrapper.py --command Stability \\
      --input assemblies/*.pdb \\
      --output-dir results/stability/ \\
      --workers 8

Note: FoldX paths are set in config.py by default
        """
    )

    parser.add_argument('-c', '--command', required=True,
                       choices=['BuildModel', 'Stability', 'AlaScan'],
                       help='FoldX command to run')
    parser.add_argument('-i', '--input', nargs='+', required=True,
                       help='Input PDB files (files, globs, or directories)')
    parser.add_argument('-o', '--output-dir', required=True,
                       help='Output directory')
    parser.add_argument('-m', '--mutations', nargs='+',
                       help='Mutations (e.g., E327D E327A) - required for BuildModel')
    parser.add_argument('--chains', nargs='+',
                       help='Chain IDs to mutate (e.g., A B C D E F) - required for BuildModel')
    parser.add_argument('--combine-mutations', action='store_true',
                       help='Apply all mutations simultaneously (combinatorial). Default: test each mutation separately')
    parser.add_argument('--foldx-path', default=DEFAULT_FOLDX_PATH,
                       help=f'Path to FoldX executable (default: from config.py or {DEFAULT_FOLDX_PATH})')
    parser.add_argument('--rotabase-path', default=DEFAULT_ROTABASE_PATH,
                       help=f'Path to FoldX rotabase.txt (default: from config.py or {DEFAULT_ROTABASE_PATH})')
    parser.add_argument('-w', '--workers', type=int, default=DEFAULT_WORKERS,
                       help='Number of parallel workers (default: CPU count - 2)')


    args = parser.parse_args()

    # Validate command-specific arguments
    if args.command == 'BuildModel':
        if not args.mutations:
            print("ERROR: --mutations required for BuildModel")
            sys.exit(1)
        if not args.chains:
            print("ERROR: --chains required for BuildModel")
            sys.exit(1)

    # Check FoldX
    if not Path(args.foldx_path).exists():
        print(f"ERROR: FoldX not found at {args.foldx_path}")
        print(f"Update path in config.py or use --foldx-path")
        sys.exit(1)
    if not Path(args.rotabase_path).exists():
        print(f"ERROR: Rotabase not found at {args.rotabase_path}")
        print(f"Update path in config.py or use --rotabase-path")
        sys.exit(1)

    # Get input files
    input_files = get_input_files(args.input)

    if not input_files:
        print("ERROR: No PDB files found")
        sys.exit(1)

    print(f"Found {len(input_files)} PDB files")
    print(f"FoldX command: {args.command}")
    print(f"FoldX path: {args.foldx_path}")

    if args.command == 'BuildModel':
        print(f"Mutations: {args.mutations}")
        print(f"Chains: {args.chains}")

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

        if args.command == 'BuildModel':
            # Get parent directory (assembly type) from input path
            input_path = Path(input_file)
            parent_dir = input_path.parent.name
            
            if args.combine_mutations:
                # Apply ALL mutations simultaneously in single structure
                mutation_name = "+".join(args.mutations)  # E327D+K155A
                
                if parent_dir and parent_dir != '.':
                    # Preserve assembly type: mutation/assembly_type/frame
                    output_subdir = Path(args.output_dir) / mutation_name / parent_dir / base_name
                else:
                    # No parent dir: mutation/frame
                    output_subdir = Path(args.output_dir) / mutation_name / base_name
                
                tasks.append({
                    'input_file': input_file,
                    'output_dir': str(output_subdir),
                    'command': args.command,
                    'mutations': args.mutations,  # Pass all mutations
                    'chains': args.chains,
                    'combine_mutations': True,
                    'foldx_path': args.foldx_path,
                    'rotabase_path': args.rotabase_path
                })
            else:
                # Create separate output directory for each mutation (default)
                for mutation in args.mutations:
                    if parent_dir and parent_dir != '.':
                        # Preserve assembly type: mutation/assembly_type/frame
                        output_subdir = Path(args.output_dir) / mutation / parent_dir / base_name
                    else:
                        # No parent dir: mutation/frame
                        output_subdir = Path(args.output_dir) / mutation / base_name
                    
                    tasks.append({
                        'input_file': input_file,
                        'output_dir': str(output_subdir),
                        'command': args.command,
                        'mutations': [mutation],
                        'chains': args.chains,
                        'combine_mutations': False,
                        'foldx_path': args.foldx_path,
                        'rotabase_path': args.rotabase_path
                    })
        elif args.command == 'AlaScan':
            # AlaScan: preserve parent directory structure like Stability
            input_path = Path(input_file)
            parent_dir = input_path.parent.name
            
            if parent_dir == '.' or parent_dir == '':
                # Files are in current directory
                output_subdir = Path(args.output_dir) / base_name
            else:
                # Use parent directory name
                output_subdir = Path(args.output_dir) / parent_dir / base_name
                
            tasks.append({
                'input_file': input_file,
                'output_dir': str(output_subdir),
                'command': args.command,
                'chains': args.chains,
                'foldx_path': args.foldx_path,
                'rotabase_path': args.rotabase_path
            })    
        else:
            # Stability: preserve parent directory structure
            input_path = Path(input_file)
            parent_dir = input_path.parent.name
            
            if parent_dir == '.' or parent_dir == '':
                # Files are in current directory
                output_subdir = Path(args.output_dir)
            else:
                # Use parent directory name
                output_subdir = Path(args.output_dir) / parent_dir
                
            tasks.append({
                'input_file': input_file,
                'output_dir': str(output_subdir),
                'command': args.command,
                'foldx_path': args.foldx_path,
                'rotabase_path': args.rotabase_path
            })

    print(f"Total tasks: {len(tasks)}")

    # Process in parallel
    results = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(process_single_foldx_task, task) for task in tasks]

        for future in tqdm(futures, desc=f"Running FoldX {args.command}"):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"Task failed: {e}")
                results.append({'status': 'error', 'error': str(e)})

    # Report results
    successful = sum(1 for r in results if r['status'] == 'success')
    failed = len(results) - successful

    print(f"\nâœ“ FoldX analysis complete: {successful}/{len(results)} successful")

    if failed > 0:
        print(f"\nâš  {failed} tasks failed:")
        for r in results:
            if r['status'] != 'success':
                print(f"  - {Path(r.get('input', 'unknown')).name}: {r.get('error', 'unknown error')}")


if __name__ == "__main__":
    main()
