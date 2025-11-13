#!/usr/bin/env python3
"""
Extract frames from MD trajectory using gmx trjconv
"""

import subprocess
import argparse
from pathlib import Path
import sys


def extract_frames(trajectory, structure, output_dir, output_prefix,
                   start=None, end=None, dt=None, selection="Protein"):
    """
    Extract frames from MD trajectory

    Args:
        trajectory: Path to trajectory file (.xtc)
        structure: Path to structure file (.tpr or .pdb)
        output_dir: Output directory
        output_prefix: Prefix for output files
        start: Start time (ps)
        end: End time (ps)
        dt: Time step (ps)
        selection: Group to extract (default: Protein)
    """

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Build output path
    output_path = Path(output_dir) / output_prefix

    # Build gmx command
    cmd = [
        'gmx', 'trjconv',
        '-f', str(trajectory),
        '-s', str(structure),
        '-o', f'{output_path}.pdb',
        '-sep'  # Separate output files
    ]

    if start is not None:
        cmd.extend(['-b', str(start)])
    if end is not None:
        cmd.extend(['-e', str(end)])
    if dt is not None:
        cmd.extend(['-dt', str(dt)])

    print(f"Extracting frames from {trajectory}")
    print(f"Command: {' '.join(cmd)}")
    print(f"Selection: {selection}")

    # Run gmx trjconv
    # Echo selection group (e.g., "Protein" = group 1 or 2, typically)
    result = subprocess.run(
        cmd,
        input=f"{selection}\n",
        text=True,
        capture_output=True
    )

    if result.returncode != 0:
        print(f"ERROR: gmx trjconv failed")
        print(f"STDERR: {result.stderr}")
        sys.exit(1)

    print(f"✓ Extraction complete")
    print(f"Output: {output_dir}/{output_prefix}.*.pdb")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Extract frames from MD trajectory using gmx trjconv",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract frames every 100 ns from 500-1000 ns
  python extract_frames.py -t traj.xtc -s topol.tpr -o frames/ -p protein --start 500000 --end 1000000 --dt 100000

  # Extract all frames
  python extract_frames.py -t traj.xtc -s topol.tpr -o frames/ -p protein
        """
    )

    parser.add_argument('-t', '--trajectory', required=True,
                       help='Trajectory file (.xtc)')
    parser.add_argument('-s', '--structure', required=True,
                       help='Structure file (.tpr or .pdb)')
    parser.add_argument('-o', '--output-dir', required=True,
                       help='Output directory')
    parser.add_argument('-p', '--output-prefix', required=True,
                       help='Prefix for output files (e.g., "protein" -> protein.1.pdb, protein.2.pdb)')
    parser.add_argument('--start', type=float,
                       help='Start time (ps)')
    parser.add_argument('--end', type=float,
                       help='End time (ps)')
    parser.add_argument('--dt', type=float,
                       help='Time step between frames (ps)')
    parser.add_argument('--selection', default='Protein',
                       help='Group selection (default: Protein)')

    args = parser.parse_args()

    # Verify input files exist
    if not Path(args.trajectory).exists():
        print(f"ERROR: Trajectory file not found: {args.trajectory}")
        sys.exit(1)
    if not Path(args.structure).exists():
        print(f"ERROR: Structure file not found: {args.structure}")
        sys.exit(1)

    extract_frames(
        args.trajectory,
        args.structure,
        args.output_dir,
        args.output_prefix,
        args.start,
        args.end,
        args.dt,
        args.selection
    )


if __name__ == "__main__":
    main()
