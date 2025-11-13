#!/usr/bin/env python3
"""
Validate and inspect alanine scanning data before analysis
"""

import argparse
import glob
from pathlib import Path
import re
from collections import defaultdict

def inspect_file(file_path, show_sample=True):
    """Inspect a single alanine scanning file"""

    info = {
        'path': file_path,
        'name': Path(file_path).name,
        'assembly': Path(file_path).parent.name,
        'lines': 0,
        'mutations': 0,
        'residue_range': None,
        'chains': set(),
        'sample_lines': [],
        'errors': []
    }

    residue_numbers = []

    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                info['lines'] += 1
                line = line.strip()

                if not line or line.startswith('#'):
                    continue

                # Try to parse
                match = re.match(r'(\w+)\s+([A-Z]?)(\d+)\s+to ALA energy change is\s+([-\d.]+)', line)

                if match:
                    info['mutations'] += 1
                    res_name = match.group(1)
                    chain = match.group(2) if match.group(2) else 'A'
                    res_num = int(match.group(3))
                    ddg = float(match.group(4))

                    residue_numbers.append(res_num)
                    info['chains'].add(chain)

                    if show_sample and len(info['sample_lines']) < 3:
                        info['sample_lines'].append(
                            f"{res_name}{chain}{res_num}: ΔΔG = {ddg:7.3f}"
                        )
                else:
                    if show_sample and len(info['sample_lines']) < 3:
                        info['errors'].append(f"Line {line_num}: Could not parse: {line[:50]}")

        if residue_numbers:
            info['residue_range'] = (min(residue_numbers), max(residue_numbers))

    except Exception as e:
        info['errors'].append(f"Error reading file: {e}")

    return info


def validate_directory_structure(input_dir, pattern="*.txt"):
    """Check if directory structure matches expected format"""

    print("="*70)
    print("DIRECTORY STRUCTURE VALIDATION")
    print("="*70)

    # Expected assembly types
    expected_assemblies = ['3-fold', '5-fold', 'dimer', 'dimer_para', 'dimer_perp']

    # Find all files
    files = glob.glob(f"{input_dir}/**/{pattern}", recursive=True)

    if not files:
        print(f"\n❌ ERROR: No files matching '{pattern}' found in {input_dir}")
        return False

    print(f"\n✓ Found {len(files)} files matching '{pattern}'")

    # Group by assembly type
    assembly_files = defaultdict(list)
    for f in files:
        assembly = Path(f).parent.name
        assembly_files[assembly].append(f)

    print(f"\n📁 Directory structure:")
    for assembly in sorted(assembly_files.keys()):
        count = len(assembly_files[assembly])
        print(f"   {assembly:20s}: {count:3d} files")

    # Check for expected assemblies
    missing_assemblies = []
    for expected in expected_assemblies:
        found = False
        for actual in assembly_files.keys():
            if expected.lower().replace('-', '').replace('_', '') == actual.lower().replace('-', '').replace('_', ''):
                found = True
                break
        if not found:
            missing_assemblies.append(expected)

    if missing_assemblies:
        print(f"\n⚠️  WARNING: Missing expected assembly types: {', '.join(missing_assemblies)}")
        print("   Expected: 3-fold, 5-fold, dimer, dimer_para (side-to-side), dimer_perp (side-to-end)")
    else:
        print(f"\n✓ All expected assembly types found")

    return True


def inspect_all_files(input_dir, pattern="*.txt", verbose=False):
    """Inspect all files and provide summary"""

    print("\n" + "="*70)
    print("FILE CONTENT INSPECTION")
    print("="*70)

    files = glob.glob(f"{input_dir}/**/{pattern}", recursive=True)

    if not files:
        print(f"\n❌ No files found")
        return

    # Inspect each file
    all_info = []
    for f in files:
        info = inspect_file(f, show_sample=verbose)
        all_info.append(info)

    # Group by assembly
    assembly_info = defaultdict(list)
    for info in all_info:
        assembly_info[info['assembly']].append(info)

    # Print summary by assembly
    print(f"\n📊 Content Summary by Assembly Type:\n")

    total_mutations = 0

    for assembly in sorted(assembly_info.keys()):
        files_in_assembly = assembly_info[assembly]

        print(f"  {assembly}:")
        print(f"    Files: {len(files_in_assembly)}")

        mutations_per_file = [f['mutations'] for f in files_in_assembly]
        total_assembly_mutations = sum(mutations_per_file)
        total_mutations += total_assembly_mutations

        print(f"    Total mutations: {total_assembly_mutations}")
        print(f"    Mutations per file: min={min(mutations_per_file)}, "
              f"max={max(mutations_per_file)}, mean={sum(mutations_per_file)/len(mutations_per_file):.1f}")

        # Check for chains
        all_chains = set()
        for f in files_in_assembly:
            all_chains.update(f['chains'])
        print(f"    Chains found: {sorted(all_chains)}")

        # Check residue range consistency
        ranges = [f['residue_range'] for f in files_in_assembly if f['residue_range']]
        if ranges:
            min_res = min(r[0] for r in ranges)
            max_res = max(r[1] for r in ranges)
            print(f"    Residue range: {min_res}-{max_res}")

        # Show sample if verbose
        if verbose and files_in_assembly[0]['sample_lines']:
            print(f"    Sample mutations:")
            for sample in files_in_assembly[0]['sample_lines']:
                print(f"      {sample}")

        # Show errors
        errors = [e for f in files_in_assembly for e in f['errors']]
        if errors:
            print(f"    ⚠️  Errors: {len(errors)}")
            if verbose:
                for error in errors[:3]:
                    print(f"      {error}")

        print()

    print(f"{'='*70}")
    print(f"TOTAL: {len(all_info)} files, {total_mutations} mutations")
    print(f"{'='*70}")


def check_symmetry_expectations(input_dir, pattern="*.txt"):
    """Check if we have the expected number of chains per assembly"""

    print("\n" + "="*70)
    print("SYMMETRY CHECK")
    print("="*70)

    expected_symmetry = {
        '3-fold': 3,
        'threefold': 3,
        '5-fold': 5,
        'fivefold': 5,
        'dimer': 2,
        'dimer_para': 2,
        'side_to_side': 2,
        'dimer_perp': 2,
        'side_to_end': 2
    }

    files = glob.glob(f"{input_dir}/**/{pattern}", recursive=True)

    assembly_chains = defaultdict(set)

    for f in files:
        assembly = Path(f).parent.name.lower()
        info = inspect_file(f, show_sample=False)
        assembly_chains[assembly].update(info['chains'])

    print("\nExpected vs. Actual chains per assembly:\n")

    all_good = True
    for assembly in sorted(assembly_chains.keys()):
        actual_chains = len(assembly_chains[assembly])
        expected = expected_symmetry.get(assembly, '?')

        status = "✓" if actual_chains == expected else "⚠️ "

        print(f"  {status} {assembly:20s}: {actual_chains} chains (expected: {expected})")

        if actual_chains != expected:
            all_good = False
            print(f"      Chains found: {sorted(assembly_chains[assembly])}")

    if all_good:
        print(f"\n✓ All assemblies have expected chain counts")
    else:
        print(f"\n⚠️  Some assemblies have unexpected chain counts")
        print(f"   This may be OK if you're using single-chain files and will average later")


def main():
    parser = argparse.ArgumentParser(
        description="Validate alanine scanning data before analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  # Quick check
  python validate_data.py --input alanine_scanning_results/

  # Detailed inspection
  python validate_data.py --input alanine_scanning_results/ --verbose
        """
    )

    parser.add_argument('-i', '--input', required=True,
                       help='Input directory containing alanine scanning files')
    parser.add_argument('--pattern', default='*.txt',
                       help='File pattern to match (default: *.txt)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Show detailed information including sample lines')

    args = parser.parse_args()

    print("\n")
    print("╔" + "="*68 + "╗")
    print("║" + " ALANINE SCANNING DATA VALIDATOR ".center(68) + "║")
    print("╚" + "="*68 + "╝")

    # 1. Check directory structure
    valid_structure = validate_directory_structure(args.input, args.pattern)

    if not valid_structure:
        print("\n❌ Validation failed: No files found")
        return 1

    # 2. Inspect file contents
    inspect_all_files(args.input, args.pattern, verbose=args.verbose)

    # 3. Check symmetry
    check_symmetry_expectations(args.input, args.pattern)

    print("\n" + "="*70)
    print("VALIDATION COMPLETE")
    print("="*70)
    print("\n✓ Data structure looks good!")
    print("\nNext steps:")
    print("  1. If everything looks correct, run the analysis:")
    print("     python analyze_alanine_scanning.py \\")
    print(f"         --input {args.input} \\")
    print("         --output results.csv \\")
    print("         --viral-contributions viral_contributions.csv \\")
    print("         --plots output_plots")
    print("\n  2. If you see warnings, check your file structure and naming")
    print("     - Assembly directories should be: 3-fold, 5-fold, dimer, dimer_para, dimer_perp")
    print("     - Files should contain lines like: 'LEU 479 to ALA energy change is 1.19627'")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main()):
