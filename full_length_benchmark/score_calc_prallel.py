import os.path as osp
import pandas as pd
import os
import re
import numpy as np
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing


def read_plddt(plddt_dir, sid):

    try:
        plddt_path = osp.join(plddt_dir, sid + '.npy')
        plddt = np.load(plddt_path)
        return plddt
    except FileNotFoundError:
        return None


def read_coords(coords_dir, sid):

    with np.load(osp.join(coords_dir, sid + '.npz')) as f:
        coords = f['coords']
        valid_mask = f['valid_mask']
    return coords.copy(), valid_mask.copy()


def step_score(diff):

    score = 0.25 * ((diff < 0.5).astype(np.float32) +
                    (diff < 1.0).astype(np.float32) +
                    (diff < 2.0).astype(np.float32) +
                    (diff < 4.0).astype(np.float32))
    return score


def parse_cigar(cigar_string, ref=1, query=1, include_m=False, order_qt=False):

    ref, query = ref-1, query-1
    matches = []

    for cnt, action in re.findall('([0-9]*)([IDMP])', cigar_string):
        cnt = int(cnt)

        if action == 'D':
            ref += cnt
        elif action == 'I':
            query += cnt
        elif action == 'M':
            if include_m:
                matches += [(ref + i, query + i) for i in range(cnt)]
            ref += cnt
            query += cnt
        elif action == 'P':
            matches += [(ref + i, query + i) for i in range(cnt)]
            ref += cnt
            query += cnt
        else:
            raise ValueError(f'Action {action}')

    if order_qt:
        return np.array(matches)[:, ::-1]
    else:
        return np.array(matches)


def lddt_multidomain(d_ij, d_kl, d_mn, idx1, cutoff):

    d = np.abs(d_kl - d_ij)

    score = 0.25 * ((d < 0.5).astype(np.float32) +
                    (d < 1.0).astype(np.float32) +
                    (d < 2.0).astype(np.float32) +
                    (d < 4.0).astype(np.float32))

    neighbor_mask = ((0 < d_ij) & (d_ij < cutoff)).astype(np.float32) * (1. - np.eye(len(d_ij)))

    score_per_residue = np.sum(neighbor_mask * score, axis=-1)
    norm_per_residue = np.sum((0 < d_mn) & (d_mn < cutoff), axis=-1)[idx1]

    return np.mean(score_per_residue / norm_per_residue), (neighbor_mask * score, norm_per_residue), (score_per_residue / norm_per_residue)


def distances_ij(q_points_aln, t_points_aln, cutoff=15):

    d_ij = np.sqrt(np.sum((q_points_aln[:, None] - q_points_aln[None, :])**2, axis=-1))
    d_kl = np.sqrt(np.sum((t_points_aln[:, None] - t_points_aln[None, :])**2, axis=-1))

    neighbor_mask = (d_ij < cutoff).astype(np.float32) * (1. - np.eye(len(d_ij)))
    return (neighbor_mask * d_ij, neighbor_mask * d_kl)


def align_coords(query, target, alignment):

    coords1, mask1 = query
    coords2, mask2 = target
    cigar, q_start, t_start = alignment

    idx1, idx2 = parse_cigar(cigar, query=q_start, ref=t_start, include_m=True, order_qt=True).T

    mask = mask1[idx1] & mask2[idx2]
    aln_coords1 = coords1[idx1, :3][mask]
    aln_coords2 = coords2[idx2, :3][mask]    
    return (aln_coords1, aln_coords2), (idx1[mask], idx2[mask]), (len(coords1), len(coords2))


def distances_mn(coords, mask):

    d_mn = np.sqrt(np.sum((coords[:, None] - coords[None, :])**2, axis=-1))
    d_mn[~mask, :] = 0
    d_mn[:, ~mask] = 0
    return d_mn


def aligned_distance_mats(aln, coords1, mask1, coords2, mask2):

    (aln_coords1, aln_coords2), (idx1, idx2), (len1, len2) = align_coords(
        (coords1, mask1), (coords2, mask2), (aln['cigar'], aln['qstart'], aln['tstart'])
    )
    mask = mask1[idx1] & mask2[idx2]
    d_mat_ij, d_mat_kl = distances_ij(aln_coords1, aln_coords2)
    d_mat_seq = idx1[None, :] - idx1[:, None]
    d_mat_mn = distances_mn(coords1[:, 0:3], mask1)
    return d_mat_ij, d_mat_kl, d_mat_seq, d_mat_mn, (idx1, idx2, len1, len2)


def lddt_all_neigh(aln, coords1, mask1, coords2, mask2):

    d_mat_ij, d_mat_kl, _, d_mn, (idx1, _, _, _) = aligned_distance_mats(aln, coords1, mask1, coords2, mask2)
    return lddt_multidomain(d_mat_ij, d_mat_kl, d_mn, idx1, 15)



def process_alignment_row(row, coord_dir, plddt_dir, score_config, plddt_threshold=50.0):
    try:
        query = row['query']
        target = row['target']
        cigar = row['cigar']
        qstart = row['qstart']
        tstart = row['tstart']

        coords1, mask1 = read_coords(coord_dir, query)
        coords2, mask2 = read_coords(coord_dir, target)
        
        valid_mask1 = ~np.any(np.isnan(coords1[:, 0:3]), axis=1)
        valid_mask2 = ~np.any(np.isnan(coords2[:, 0:3]), axis=1)
        mask1 = mask1 & valid_mask1
        mask2 = mask2 & valid_mask2

        plddt1 = read_plddt(plddt_dir, query)
        plddt2 = read_plddt(plddt_dir, target)

        (aln_coords1, aln_coords2), (idx1, idx2), (len1, len2) = align_coords(
            (coords1, mask1),
            (coords2, mask2),
            (cigar, qstart, tstart)
        )
        
        example_alignment = {
            'query': query,
            'target': target,
            'cigar': cigar,
            'qstart': qstart,
            'tstart': tstart
        }
        
        qlen = len(coords1[mask1])
        tlen = len(coords2[mask2])
        aln_len = len(aln_coords1)
        
        results = {
            'qlen_calc': qlen,
            'tlen_calc': tlen,
            'aln_len_calc': aln_len,
        }
        
        
        if score_config.get('lddt_all', False) or score_config.get('lddt_all_filtered', False):
            lddt_md, (score_mat, norm_per_residue), lddt_res_all = lddt_all_neigh(
                example_alignment, coords1, mask1, coords2, mask2
            )
            if score_config.get('lddt_all', False):
                results['lddt_all'] = np.round(lddt_md, 6)
            
            if score_config.get('lddt_all_filtered', False):
                if plddt1 is not None:
                    
                    plddt_mask = plddt1 >= plddt_threshold
                    aln_plddt_mask = plddt_mask[idx1]
                    
                    if np.sum(aln_plddt_mask) > 0:
                        lddt_all_filtered = np.mean(lddt_res_all[aln_plddt_mask])
                        results['lddt_all_filtered'] = np.round(lddt_all_filtered, 6)
                        results['qlen_all_filtered'] = np.sum(plddt_mask)
                        results['aln_len_all_filtered'] = np.sum(aln_plddt_mask)
                    else:
                        results['lddt_all_filtered'] = np.nan
                        results['qlen_all_filtered'] = 0
                        results['aln_len_all_filtered'] = 0
                else:
                    results['lddt_all_filtered'] = np.nan
                    results['qlen_all_filtered'] = np.nan
                    results['aln_len_all_filtered'] = np.nan
        
        
        # Explicitly free memory
        del coords1, coords2, plddt1, plddt2, aln_coords1, aln_coords2
        
        return pd.Series(results)
        
    except Exception as e:
        print(f"Error processing row: {e}")
        results = {
            'qlen_calc': np.nan,
            'tlen_calc': np.nan,
            'aln_len_calc': np.nan,
        }
        
        for score_name in score_config:
            if score_config[score_name]:
                results[score_name] = np.nan
        
        return pd.Series(results)


# ============================================================================
# MAIN EXECUTION
# ============================================================================


if __name__ == "__main__":
    print("=" * 60)
    print("Calculating alignment scores (PARALLEL STREAMING MODE)")
    print("=" * 60)

    # ========================================================================
    # CONFIGURATION
    # ========================================================================
    coord_dir = '/alphafold_coords/'
    plddt_dir = '/alphafold_plddt_scores/'    

    score_config = {
        'lddt_all': True,
        'lddt_all_filtered': True,
    }

    input_path ="./dali.tsv"
    output_path ="./dali_calc_full"

    df_aln = pd.read_csv(input_path, sep='\t')
    total = len(df_aln)
    print(f"Processing {total} alignments in parallel...")

    num_workers = min(100, multiprocessing.cpu_count()) # adjust this

    print(f"Using {num_workers} parallel workers")

    header_written = os.path.exists(output_path)
    if not header_written:
        print("Creating header in output file...")
        dummy_result = process_alignment_row(df_aln.iloc[0], coord_dir, plddt_dir, score_config)
        column_order = list(df_aln.columns) + list(dummy_result.index)
        with open(output_path, 'w') as f:
            f.write('\t'.join(column_order) + '\n')
    else:
        print("Appending to existing output file...")
        with open(output_path, 'r') as f:
            first_line = f.readline().strip()
            column_order = first_line.split('\t')

    # ========================================================================
    # PARALLEL PROCESSING WITH ORDER PRESERVATION
    # ========================================================================
    batch_size = 2000
    
    results_dict = {}
    next_index_to_write = 0
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(
                process_alignment_row,
                row,
                coord_dir,
                plddt_dir,
                score_config,
                50.0
            ): i for i, (_, row) in enumerate(df_aln.iterrows())
        }

        for i, future in enumerate(as_completed(futures)):
            idx = futures[future]
            try:
                res = future.result()
            except Exception as e:
                print(f"Error in alignment {idx}: {e}")
                res = pd.Series({})

            results_dict[idx] = res

            while next_index_to_write in results_dict:
                results_buffer = []
                index_buffer = []
                
                while (next_index_to_write in results_dict and 
                       len(results_buffer) < batch_size):
                    results_buffer.append(results_dict[next_index_to_write])
                    index_buffer.append(next_index_to_write)
                    del results_dict[next_index_to_write]
                    next_index_to_write += 1
                
                if results_buffer:
                    results_df = pd.DataFrame(results_buffer)
                    input_df_chunk = df_aln.iloc[index_buffer].reset_index(drop=True)
                    part_df = pd.concat([input_df_chunk, results_df], axis=1)

                    for col in column_order:
                        if col not in part_df.columns:
                            part_df[col] = np.nan
                    part_df = part_df[column_order]
                    part_df.to_csv(output_path, sep='\t', index=False, mode='a', header=False)

            if (i + 1) % 1000 == 0:
                print(f"Processed {i+1}/{total} alignments...")

    if results_dict:
        print(f"Writing {len(results_dict)} remaining results...")
        remaining_indices = sorted(results_dict.keys())
        results_buffer = [results_dict[idx] for idx in remaining_indices]
        
        results_df = pd.DataFrame(results_buffer)
        input_df_chunk = df_aln.iloc[remaining_indices].reset_index(drop=True)
        part_df = pd.concat([input_df_chunk, results_df], axis=1)

        for col in column_order:
            if col not in part_df.columns:
                part_df[col] = np.nan
        part_df = part_df[column_order]

        part_df.to_csv(output_path, sep='\t', index=False, mode='a', header=False)

    print(f"All results written incrementally to: {output_path}")
    print("Done.")