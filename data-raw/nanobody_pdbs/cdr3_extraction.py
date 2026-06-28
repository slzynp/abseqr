#!/usr/bin/env python3
"""
Von Mises VAE - CDR3 Extraction Pipeline
- Calculates torsion angles (phi, psi)
- UPDATED: Direct 2D angles in radians
- Ready for Von Mises VAE training
"""

import os
import numpy as np
import pandas as pd
import torch
from pathlib import Path

from Bio import SeqIO, PDB


def parse_fasta_with_pdb_id(fasta_file):
    """
    Parse FASTA file and extract PDB ID from header
    
    Returns:
    --------
    sequences : dict
        {pdb_id: {'chain': 'H', 'sequence': 'QVQL...', 'resolution': '2.59'}}
    """
    
    sequences = {}
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        header = record.id  # e.g., "9rxi_H"
        sequence = str(record.seq)
        
        # Parse header
        parts = header.split()
        pdb_info = parts[0].split('_')
        pdb_id = pdb_info[0].lower()
        chain = pdb_info[1] if len(pdb_info) > 1 else 'A'
        
        resolution = None
        if len(parts) > 1 and 'res:' in parts[1]:
            resolution = parts[1].split(':')[1]
        
        sequences[pdb_id] = {
            'chain': chain,
            'sequence': sequence,
            'resolution': resolution,
        }
    
    return sequences


def calculate_torsion_angles(residue, prev_res=None, next_res=None):
    """
    Calculate phi, psi angles (in radians, not degrees!)
    
    Returns:
    --------
    angles : dict with 'phi', 'psi' (in radians, range [-π, π])
    """
    angles = {}
    
    try:
        # PHI: C(-1) - N - CA - C
        if prev_res and 'C' in prev_res and \
           'N' in residue and 'CA' in residue and 'C' in residue:
            phi = PDB.vectors.calc_dihedral(
                prev_res['C'].get_vector(),
                residue['N'].get_vector(),
                residue['CA'].get_vector(),
                residue['C'].get_vector()
            )
            
            angles['phi'] = float(phi)
        
        # PSI: N - CA - C - N(+1)
        if next_res and 'N' in next_res and \
           'N' in residue and 'CA' in residue and 'C' in residue:
            psi = PDB.vectors.calc_dihedral(
                residue['N'].get_vector(),
                residue['CA'].get_vector(),
                residue['C'].get_vector(),
                next_res['N'].get_vector()
            )
            angles['psi'] = float(psi)
    
    except:
        pass
    
    return angles


def extract_cdr3_from_pdb_chothia(pdb_file, chain):
    """
    Extract CDR3 from PDB file using Chothia numbering (105-117)
    
    Returns:
    --------
    cdr3_residues : list of Bio.PDB.Residue
    """
    
    try:
        parser = PDB.PDBParser(QUIET=True)
        struct = parser.get_structure('pdb', pdb_file)
        chain_obj = struct[0][chain]
    except Exception as e:
        return None
    
    cdr3_residues = []
    
    # Chothia CDR3: residues 105-117
    cdr3_start = 105
    cdr3_end = 117
    
    for res_num in range(cdr3_start, cdr3_end + 1):
        try:
            residue = chain_obj[res_num]
            if PDB.is_aa(residue):
                cdr3_residues.append(residue)
        except:
            pass  # Residue not found
    
    return cdr3_residues if len(cdr3_residues) > 0 else None


def process_nanobodies_2d_angles(fasta_file, pdb_dir, output_csv='cdr3_data.csv'):
    """
    Main pipeline: extract CDR3 and calculate 2D angles [φ, ψ] in radians
    
    UPDATED for Von Mises: Output is 2D angles (φ, ψ) in radians, not 4D sin-cos
    
    Parameters:
    -----------
    fasta_file : str
        Path to FASTA file with sequences
    pdb_dir : str
        Directory containing PDB files (Chothia numbered)
    output_csv : str
        Output CSV file with results
    
    Returns:
    --------
    results : list of dicts
    df : pandas DataFrame
    """
    
    print(f"\n{'='*90}")
    print(f"CDR3 EXTRACTION PIPELINE (Von Mises - 2D angles)")
    print(f"{'='*90}\n")
    
    # Parse FASTA
    print(f"1. Parsing FASTA file: {fasta_file}")
    sequences = parse_fasta_with_pdb_id(fasta_file)
    print(f"   Found {len(sequences)} nanobody sequences\n")
    
    # Find PDB files
    print(f"2. Scanning PDB directory: {pdb_dir}")
    pdb_files = {}
    for pdb_file in Path(pdb_dir).glob('*.pdb'):
        pdb_id = pdb_file.stem.lower().split('_')[0]
        pdb_files[pdb_id] = str(pdb_file)
    print(f"   Found {len(pdb_files)} PDB files\n")
    
    # Process each nanobody
    print(f"3. Extracting CDR3 regions (105-117)...")
    print(f"-{'*'*89}\n")
    
    results = []
    stats = {
        'total': 0,
        'pdb_found': 0,
        'cdr3_extracted': 0,
        'angles_calculated': 0,
        'lengths': [],
    }
    
    for pdb_id, seq_info in sequences.items():
        stats['total'] += 1
        
        # Find PDB file
        if pdb_id not in pdb_files:
            print(f"  ✗ {pdb_id}: PDB file not found")
            continue
        
        stats['pdb_found'] += 1
        pdb_file = pdb_files[pdb_id]
        chain = seq_info['chain']
        
        # Extract CDR3 from PDB
        cdr3_residues = extract_cdr3_from_pdb_chothia(pdb_file, chain)
        
        if not cdr3_residues:
            print(f"  ✗ {pdb_id}: CDR3 not extracted")
            continue
        
        stats['cdr3_extracted'] += 1
        cdr3_length = len(cdr3_residues)
        stats['lengths'].append(cdr3_length)
        
        # Calculate angles for each residue
        phi_list = []
        psi_list = []
        
        for i, res in enumerate(cdr3_residues):
            prev_res = cdr3_residues[i - 1] if i > 0 else None
            next_res = cdr3_residues[i + 1] if i < len(cdr3_residues) - 1 else None
            
            angles = calculate_torsion_angles(res, prev_res, next_res)
            
            # UPDATED: Direct angles in radians
            phi = angles.get('phi')
            psi = angles.get('psi')
            
            if phi is not None:
                phi_list.append(phi)
            else:
                phi_list.append(0.0)  # Default for missing phi
            
            if psi is not None:
                psi_list.append(psi)
            else:
                psi_list.append(0.0)  # Default for missing psi
        
        stats['angles_calculated'] += 1
        
        result = {
            'pdb_id': pdb_id,
            'chain': chain,
            'cdr3_length': cdr3_length,
            'phi': phi_list,  # List of φ angles in radians
            'psi': psi_list,  # List of ψ angles in radians
            'resolution': seq_info['resolution'],
        }
        
        results.append(result)
        print(f"  ✓ {pdb_id}: {cdr3_length} residues, {len(phi_list)} angles")
    
    print(f"\n{'='*90}")
    print(f"STATISTICS:")
    print(f"{'='*90}\n")
    
    print(f"  Total nanobodies: {stats['total']}")
    print(f"  PDB files found: {stats['pdb_found']}")
    print(f"  CDR3 extracted: {stats['cdr3_extracted']}")
    print(f"  Angles calculated: {stats['angles_calculated']}")
    
    if stats['lengths']:
        lengths_arr = np.array(stats['lengths'])
        print(f"\n  CDR3 LENGTH DISTRIBUTION:")
        print(f"    Min: {lengths_arr.min()}")
        print(f"    Max: {lengths_arr.max()}")
        print(f"    Mean: {lengths_arr.mean():.1f}")
        print(f"    Median: {np.median(lengths_arr):.1f}")
        print(f"    95th percentile: {np.percentile(lengths_arr, 95):.0f}")
        print(f"    99th percentile: {np.percentile(lengths_arr, 99):.0f}")
        
        optimal_max_len = int(np.percentile(lengths_arr, 95)) + 2
        print(f"\n  RECOMMENDED max_len: {optimal_max_len}")
    
    print(f"{'='*90}\n")
    
    # Save to CSV
    print(f"4. Saving results to {output_csv}...")
    
    # Create DataFrame
    data_rows = []
    for result in results:
        pdb_id = result['pdb_id']
        
        for i in range(len(result['phi'])):
            data_rows.append({
                'pdb_id': pdb_id,
                'chain': result['chain'],
                'cdr3_length': result['cdr3_length'],
                'position': i + 1,
                'phi': result['phi'][i],  # UPDATED: Direct angle in radians
                'psi': result['psi'][i],  # UPDATED: Direct angle in radians
                'resolution': result['resolution'],
            })
    
    df = pd.DataFrame(data_rows)
    df.to_csv(output_csv, index=False)
    print(f"   ✓ Saved {len(df)} rows to {output_csv}\n")
    
    return results, df, stats['lengths']


def create_vae_dataset_2d(csv_file, max_len=None, train_split=0.8, random_seed=42, 
                          output_train_csv='cdr3_train.csv', output_test_csv='cdr3_test.csv'):
    """
    Create tensor dataset for Von Mises VAE training with train/test split
    
    
    Parameters:
    -----------
    csv_file : str
        Path to CSV file with CDR3 angles
    max_len : int or None
        Maximum sequence length for padding. If None, uses 95th percentile + 2
    train_split : float
        Fraction of PDB structures for training (default 0.8)
    random_seed : int
        Random seed for reproducibility
    output_train_csv : str
        Output path for training CSV
    output_test_csv : str
        Output path for test CSV
    
    Returns:
    --------
    datasets : dict with keys:
        'train_data': torch.Tensor (N_train, max_len, 2)  # UPDATED: 2D not 4D!
        'train_lengths': torch.Tensor (N_train,)
        'train_pdb_ids': list of PDB IDs
        'test_data': torch.Tensor (N_test, max_len, 2)   # UPDATED: 2D not 4D!
        'test_lengths': torch.Tensor (N_test,)
        'test_pdb_ids': list of PDB IDs
        'max_len': int
    """
    
    df = pd.read_csv(csv_file)
    np.random.seed(random_seed)
    
    # Get unique PDB IDs and their lengths
    unique_pdbs = df['pdb_id'].unique()
    pdb_lengths = df.groupby('pdb_id').size().values
    
    if max_len is None:
        max_len = int(np.percentile(pdb_lengths, 95)) + 2
        print(f"\nAuto-selected max_len = {max_len} (95th percentile + 2)")
    
    # Split PDB IDs (not individual samples)
    unique_pdbs = np.array(unique_pdbs)
    np.random.shuffle(unique_pdbs)
    split_idx = int(len(unique_pdbs) * train_split)
    train_pdb_ids = list(unique_pdbs[:split_idx])
    test_pdb_ids = list(unique_pdbs[split_idx:])
    
    print(f"\nTrain/Test split (by PDB structure):")
    print(f"  Total unique PDBs: {len(unique_pdbs)}")
    print(f"  Train PDBs: {len(train_pdb_ids)} ({100*len(train_pdb_ids)/len(unique_pdbs):.1f}%)")
    print(f"  Test PDBs: {len(test_pdb_ids)} ({100*len(test_pdb_ids)/len(unique_pdbs):.1f}%)")
    print(f"  Random seed: {random_seed}")
    print(f"  → Same PDB IDs stay together (no data leakage)")
    
    # Split CSV data
    train_df = df[df['pdb_id'].isin(train_pdb_ids)]
    test_df = df[df['pdb_id'].isin(test_pdb_ids)]
    
    # Save CSV files
    train_df.to_csv(output_train_csv, index=False)
    test_df.to_csv(output_test_csv, index=False)
    
    print(f"\n  Saved train CSV: {output_train_csv} ({len(train_df)} rows)")
    print(f"  Saved test CSV: {output_test_csv} ({len(test_df)} rows)")
    
    # UPDATED: Features are now 2D angles [phi, psi]
    feature_cols = ['phi', 'psi']
    
    def create_split_tensors(pdb_ids, label="Train"):
        """Helper function to create tensors for a set of PDB IDs"""
        all_features = []
        length_list = []
        pdb_list = []
        
        for pdb_id in pdb_ids:
            group = df[df['pdb_id'] == pdb_id].sort_values('position')
            
            # Extract features
            features = group[feature_cols].values  # Shape: (actual_len, 2)
            
            actual_len = len(features)
            length_list.append(actual_len)
            pdb_list.append(pdb_id)
            
            # Pad or truncate
            if actual_len > max_len:
                features = features[:max_len]
            
            # Pad with zeros
            if actual_len < max_len:
                padding = np.zeros((max_len - actual_len, 2))
                features = np.vstack([features, padding])
            
            all_features.append(features)
        
        data_tensor = torch.tensor(np.array(all_features), dtype=torch.float32)
        length_tensor = torch.tensor(length_list, dtype=torch.long)
        
        print(f"\n{label} set:")
        print(f"  Samples: {len(data_tensor)}")
        print(f"  Shape: {data_tensor.shape}")  # Should be (N, max_len, 2)
        print(f"  Features: phi, psi (2D angles in radians)")
        print(f"  Range: [{data_tensor.min():.3f}, {data_tensor.max():.3f}]")
        
        return data_tensor, length_tensor, pdb_list
    
    # Create train and test tensors
    train_data, train_lengths, train_pdbs = create_split_tensors(train_pdb_ids, "Train")
    test_data, test_lengths, test_pdbs = create_split_tensors(test_pdb_ids, "Test")
    
    # Package results
    datasets = {
        'train_data': train_data,
        'train_lengths': train_lengths,
        'train_pdb_ids': train_pdbs,
        'test_data': test_data,
        'test_lengths': test_lengths,
        'test_pdb_ids': test_pdbs,
        'max_len': max_len,
        'features': ['phi', 'psi'],  # UPDATED: 2D not 4D
    }
    
    return datasets


# ============================================================================
# USAGE
# ============================================================================

def main():
    fasta_file = '/Users/slaakyz/Desktop/edu/spring 26\' /project/final_clustered_nanobodies.fasta'
    pdb_dir = '/Users/slaakyz/Desktop/edu/spring 26\' /project/nanobody_pdbs/chothia' 
    output_csv = 'cdr3_data.csv'
    
    # Step 1: Extract CDR3 and calculate 2D angles (φ, ψ) in radians
    print("\n✓ UPDATED for Von Mises: Using 2D angles [φ, ψ] in radians\n")
    _, _, _ = process_nanobodies_2d_angles(fasta_file, pdb_dir, output_csv)
    
    # Step 2: Create VAE dataset with train/test split (2D angles)
    datasets = create_vae_dataset_2d(
        output_csv, 
        train_split=0.8, 
        random_seed=42,
        output_train_csv='cdr3_train.csv',
        output_test_csv='cdr3_test.csv'
    )
    
    # Step 3: Save tensors
    torch.save({
        'train_data': datasets['train_data'],
        'train_lengths': datasets['train_lengths'],
        'train_pdb_ids': datasets['train_pdb_ids'],
        'test_data': datasets['test_data'],
        'test_lengths': datasets['test_lengths'],
        'test_pdb_ids': datasets['test_pdb_ids'],
        'max_len': datasets['max_len'],
        'features': datasets['features'],
    }, 'cdr3_train_test_split.pt')
    
    print(f"\n✓ All done! Train/test split created (VON MISES 2D)")
    print(f"\nFiles generated:")
    print(f"  CSVs:")
    print(f"    - cdr3_train.csv (training data with 2D angles)")
    print(f"    - cdr3_test.csv (test data with 2D angles)")
    print(f"  Tensors:")
    print(f"    - cdr3_train_test_split.pt (PyTorch tensors)")
    print(f"\nUsage in Von Mises VAE script:")
    print(f"  checkpoint = torch.load('cdr3_train_test_split.pt')")
    print(f"  train_data = checkpoint['train_data']        # Shape: (N_train, max_len, 2)")
    print(f"  train_lengths = checkpoint['train_lengths']")
    print(f"  train_pdb_ids = checkpoint['train_pdb_ids']")
    print(f"  test_data = checkpoint['test_data']          # Shape: (N_test, max_len, 2)")
    print(f"  test_lengths = checkpoint['test_lengths']")
    print(f"  test_pdb_ids = checkpoint['test_pdb_ids']")
    print(f"  max_len = checkpoint['max_len']")
    print(f"  features = checkpoint['features']            # ['phi', 'psi']")
    print(f"\nTo load CSVs:")
    print(f"  train_df = pd.read_csv('cdr3_train.csv')     # Columns: phi, psi")
    print(f"  test_df = pd.read_csv('cdr3_test.csv')")


if __name__ == "__main__":
    main()