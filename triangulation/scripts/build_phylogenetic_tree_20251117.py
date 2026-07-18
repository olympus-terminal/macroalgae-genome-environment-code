#!/usr/bin/env python3
"""
Build Phylogenetic Tree from rbcL Alignment
Created: 2025-11-17

Uses FastTree to construct a phylogenetic tree from the rbcL alignment.
This tree can be used to visualize triangulation results in phylogenetic context.

Software versions logged.
"""

import subprocess
import os
import sys
from datetime import datetime
from pathlib import Path

print('=' * 80)
print('PHYLOGENETIC TREE CONSTRUCTION')
print('=' * 80)
print(f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
ALIGN_FILE = str(BASE_DIR / 'DATA' / 'rbcL_alignment.fa')
TRIANG_DIR = str(BASE_DIR / 'triangulation')
TREE_DIR = os.path.join(TRIANG_DIR, 'phylogeny')
os.makedirs(TREE_DIR, exist_ok=True)

# Output files
TIMESTAMP = datetime.now().strftime('%Y%m%d_%H%M%S')
TREE_FILE = os.path.join(TREE_DIR, f'rbcL_phylogenetic_tree_{TIMESTAMP}.nwk')
LOG_FILE = os.path.join(TREE_DIR, f'tree_construction_log_{TIMESTAMP}.txt')

print(f'\nInput alignment: {ALIGN_FILE}')
print(f'Output tree: {TREE_FILE}')
print(f'Log file: {LOG_FILE}')

# Check if alignment exists
if not os.path.exists(ALIGN_FILE):
    print(f'\n❌ ERROR: Alignment file not found: {ALIGN_FILE}')
    sys.exit(1)

# Count sequences
print('\nCounting sequences...')
with open(ALIGN_FILE, 'r') as f:
    n_seqs = sum(1 for line in f if line.startswith('>'))
print(f'  Found {n_seqs} sequences')

# Check FastTree version
print('\nChecking FastTree installation...')
try:
    ft_version = subprocess.run(['FastTree'], capture_output=True, text=True, timeout=5)
    ft_version_line = ft_version.stderr.split('\n')[0] if ft_version.stderr else 'Version unknown'
    print(f'  FastTree: {ft_version_line}')
except Exception as e:
    print(f'  ❌ ERROR: FastTree not found: {e}')
    sys.exit(1)

# Build tree with FastTree
print('\n' + '=' * 80)
print('Running FastTree (this may take a few minutes)...')
print('=' * 80)

# FastTree command for protein alignment
cmd = [
    'FastTree',
    '-log', LOG_FILE,
    ALIGN_FILE
]

print(f'\nCommand: {" ".join(cmd)}')
print('\nBuilding tree...')

try:
    with open(TREE_FILE, 'w') as tree_out:
        result = subprocess.run(
            cmd,
            stdout=tree_out,
            stderr=subprocess.PIPE,
            text=True,
            timeout=600  # 10 minute timeout
        )

    if result.returncode == 0:
        print('✓ Tree construction completed successfully')
    else:
        print(f'⚠ Warning: FastTree returned code {result.returncode}')
        print(f'STDERR: {result.stderr}')

except subprocess.TimeoutExpired:
    print('❌ ERROR: Tree construction timed out (>10 minutes)')
    sys.exit(1)
except Exception as e:
    print(f'❌ ERROR: Tree construction failed: {e}')
    sys.exit(1)

# Verify tree file
if os.path.exists(TREE_FILE) and os.path.getsize(TREE_FILE) > 0:
    file_size = os.path.getsize(TREE_FILE)
    print(f'\n✓ Tree file created: {TREE_FILE}')
    print(f'  Size: {file_size} bytes')

    # Read and validate tree
    with open(TREE_FILE, 'r') as f:
        tree_str = f.read()
        n_leaves = tree_str.count(',') + 1  # Rough estimate
        print(f'  Estimated leaves: {n_leaves}')

        # Check if tree looks valid
        if '(' in tree_str and ')' in tree_str and ';' in tree_str:
            print('  ✓ Tree format appears valid (Newick)')
        else:
            print('  ⚠ Warning: Tree format may be invalid')
else:
    print(f'\n❌ ERROR: Tree file not created or empty')
    sys.exit(1)

# Read log file for statistics
print('\n' + '=' * 80)
print('TREE CONSTRUCTION STATISTICS')
print('=' * 80)

if os.path.exists(LOG_FILE):
    with open(LOG_FILE, 'r') as f:
        log_content = f.read()

    # Extract key statistics
    for line in log_content.split('\n'):
        if 'Total time:' in line or 'Gamma20' in line or 'TreeLength' in line:
            print(f'  {line.strip()}')

    print(f'\n✓ Full log saved to: {LOG_FILE}')

print('\n' + '=' * 80)
print('PHYLOGENETIC TREE CONSTRUCTION COMPLETE')
print('=' * 80)
print(f'\nOutput files:')
print(f'  Tree: {TREE_FILE}')
print(f'  Log:  {LOG_FILE}')
print(f'\nNext steps:')
print(f'  1. Visualize tree with sample metadata')
print(f'  2. Map triangulation results onto phylogeny')
print(f'  3. Compare analytical results across the phylogeny')
