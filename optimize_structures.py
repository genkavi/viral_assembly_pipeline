#!/usr/bin/env python3
"""
Optimize PDB structures: Fix residue names (CHARMM→PDB) and optionally repair with FoldX
"""

import argparse
import glob
import os
import subprocess
import sys
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


def fix_residue_names(input_file, output_file,
                     his_rename_map=None, fix_ile_cd=True):
    """
    Fix residue names from CHARMM to standard PDB naming
    """

    if his_rename_map is None:
        try:
            from config import HIS_RENAME_MAP
            his_rename_map = HIS_RENAME_MAP
        except ImportError:
            his_rename_map = {
                'HSE': 'HIS', 'HSD': 'HIS', 'HSP': 'HIS',
                'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS',
                'HIS': 'HIS'
            }

    his_changes = 0
    ile_changes = 0

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith(('ATOM', 'HETATM')):
                residue_name = line[17:20].strip()
                atom_name = line[12:16].strip()

                # Fix histidine names (CHARMM → PDB)
                if residue_name in his_rename_map:
                    new_name = his_rename_map[residue_name]
                    new_line = line[:17] + f"{new_name:>3}" + line[20:]
                    outfile.write(new_line)
                    his_changes += 1
                # Fix ILE CD → CD1 (CHARMM → PDB)
                elif fix_ile_cd and residue_name == 'ILE' and atom_name == 'CD':
                    new_line = line[:12] + f"{'CD1':>4}" + line[16:]
                    outfile.write(new_line)
                    ile_changes += 1
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

    return his_changes, ile_changes


def repair_with_foldx(input_file, output_file, foldx_path, rotabase_path):
    """
    Repair structure with FoldX RepairPDB

    Args:
        input_file: Input PDB file
        output_file: Expected output PDB file
        foldx_path: Path to FoldX executable
        rotabase_path: Path to rotabase.txt
    """

    work_dir = Path(input_file).parent
    pdb_name = Path(input_file).name

    # Expand ~ in paths
    foldx_path = os.path.expanduser(foldx_path)
    rotabase_path = os.path.expanduser(rotabase_path)

    cmd = [
        foldx_path,
        "--rotabaseLocation", rotabase_path,
        "--command=RepairPDB",
        f"--pdb={pdb_name}",
        "--numberOfRuns=1"
    ]

    # DEBUG: Uncomment to see what's being run
    # print(f"Working dir: {work_dir}")
    # print(f"Command: {' '.join(cmd)}")

    result = subprocess.run(
        cmd,
        cwd=work_dir,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        return False, result.stderr

    # FoldX creates file with _Repair suffix
    stem = Path(pdb_name).stem
    repaired_file = work_dir / f"{stem}_Repair.pdb"

    if repaired_file.exists():
        return True, None
    else:
        return False, f"Repaired file not created: {repaired_file}"


def process_single_file(task):
    """Process a single PDB file"""

    input_file = task['input_file']
    output_file = task['output_file']
    repair = task['repair']
    foldx_path = task.get('foldx_path')
    rotabase_path = task.get('rotabase_path')

    try:
        # Step 1: Fix residue names
        his_changes, ile_changes = fix_residue_names(input_file, output_file)

        # Step 2: Repair with FoldX (if requested)
        if repair:
            # FoldX expects specific naming
            stem = Path(output_file).stem
            repaired_file = Path(output_file).parent / f"{stem}_Repair.pdb"

            success, error = repair_with_foldx(
                output_file, repaired_file, foldx_path, rotabase_path
            )

            if success:
                # Replace original with repaired
                os.replace(repaired_file, output_file)
                return {
                    'status': 'success',
                    'input': input_file,
                    'output': output_file,
                    'his_changes': his_changes,
                    'ile_changes': ile_changes,
                    'repaired': True
                }
            else:
                return {
                    'status': 'failed',
                    'input': input_file,
                    'error': f"FoldX repair failed: {error}"
                }

        return {
            'status': 'success',
            'input': input_file,
            'output': output_file,
            'his_changes': his_changes,
            'ile_changes': ile_changes,
            'repaired': False
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
        description="Optimize PDB structures: Fix residue names and optionally repair with FoldX",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Fix residue names only
  python optimize_structures.py --input frames/*.pdb --output-dir frames_fixed/

  # Fix residue names and repair with FoldX
  python optimize_structures.py --input frames/ --output-dir frames_opt/ --repair --workers 8

  # With custom suffix
  python optimize_structures.py --input frames/*.pdb --output-dir optimized/ --suffix _opt --repair

Note: FoldX paths are set in config.py by default
        """
    )

    parser.add_argument('-i', '--input', nargs='+', required=True,
                       help='Input PDB files (files, globs, or directories)')
    parser.add_argument('-o', '--output-dir', required=True,
                       help='Output directory')
    parser.add_argument('--suffix', default='',
                       help='Suffix to add to output filenames (before .pdb)')
    parser.add_argument('--repair', action='store_true',
                       help='Repair structures with FoldX after fixing residue names')
    parser.add_argument('--foldx-path', default=DEFAULT_FOLDX_PATH,
                       help=f'Path to FoldX executable (default: from config.py or {DEFAULT_FOLDX_PATH})')
    parser.add_argument('--rotabase-path', default=DEFAULT_ROTABASE_PATH,
                       help=f'Path to FoldX rotabase.txt (default: from config.py or {DEFAULT_ROTABASE_PATH})')
    parser.add_argument('-w', '--workers', type=int, default=DEFAULT_WORKERS,
                       help='Number of parallel workers (default: CPU count - 2)')

    args = parser.parse_args()

    # Get input files - THIS LINE WAS MISSING!
    input_files = get_input_files(args.input)

    if not input_files:
        print("ERROR: No PDB files found")
        sys.exit(1)

    print(f"Found {len(input_files)} PDB files")

    # Check FoldX if repair requested
    if args.repair:
        if not Path(args.foldx_path).exists():
            print(f"ERROR: FoldX not found at {args.foldx_path}")
            print(f"Update path in config.py or use --foldx-path")
            sys.exit(1)
        if not Path(args.rotabase_path).exists():
            print(f"ERROR: Rotabase not found at {args.rotabase_path}")
            print(f"Update path in config.py or use --rotabase-path")
            sys.exit(1)

    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # Determine number of workers
    if args.workers is None:
        workers = max(1, mp.cpu_count() - 2)
    else:
        workers = args.workers

    print(f"Using {workers} parallel workers")
    if args.repair:
        print(f"FoldX repair enabled")
        print(f"FoldX path: {args.foldx_path}")

    # Create tasks - SMART DIRECTORY PRESERVATION
    tasks = []

    # First, check if we have multiple parent directories
    parent_dirs = set(Path(f).parent.name for f in input_files)

    # If only one parent directory across all files, don't preserve it (flatten)
    # If multiple parent directories, preserve them
    preserve_dirs = len(parent_dirs) > 1

    if preserve_dirs:
        print(f"Preserving {len(parent_dirs)} subdirectories: {', '.join(sorted(parent_dirs))}")
    else:
        print(f"Flattening output (single input directory: {parent_dirs.pop() if parent_dirs else 'none'})")

    for input_file in input_files:
        input_path = Path(input_file)
        base_name = input_path.stem

        # Preserve parent directory only if we have multiple different parents
        if preserve_dirs:
            parent_dir = input_path.parent.name
            output_subdir = Path(args.output_dir) / parent_dir
        else:
            # Flatten - output directly to output_dir
            output_subdir = Path(args.output_dir)

        output_subdir.mkdir(parents=True, exist_ok=True)

        if args.suffix:
            output_name = f"{base_name}{args.suffix}.pdb"
        else:
            output_name = f"{base_name}.pdb"

        output_file = output_subdir / output_name

        tasks.append({
            'input_file': input_file,
            'output_file': str(output_file),
            'repair': args.repair,
            'foldx_path': args.foldx_path,
            'rotabase_path': args.rotabase_path
        })


    # Process in parallel
    results = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(process_single_file, task) for task in tasks]

        for future in tqdm(futures, desc="Optimizing structures"):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"Task failed: {e}")
                results.append({'status': 'error', 'error': str(e)})

    # Report results
    successful = sum(1 for r in results if r['status'] == 'success')
    failed = len(results) - successful

    print(f"\n✓ Optimization complete: {successful}/{len(results)} successful")

    if failed > 0:
        print(f"\n⚠ {failed} files failed:")
        for r in results:
            if r['status'] != 'success':
                print(f"  - {r.get('input', 'unknown')}: {r.get('error', 'unknown error')}")

    # Summary statistics
    if successful > 0:
        total_his = sum(r.get('his_changes', 0) for r in results if r['status'] == 'success')
        total_ile = sum(r.get('ile_changes', 0) for r in results if r['status'] == 'success')
        print(f"\nResidue name changes:")
        print(f"  - Histidine (HSE/HSD/HSP → HIS): {total_his}")
        print(f"  - Isoleucine (CD → CD1): {total_ile}")



if __name__ == "__main__":
    main()
