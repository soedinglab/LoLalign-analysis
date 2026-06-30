#!/usr/bin/env python3
"""Superset parquet = v12 + topo15 (5-bin) + qsep6 + qfrac6 (both 6-bin).

Finer long-range binning: SEP_EDGES6 = [5,13,31,51,101]
  -> 1-4, 5-12, 13-30, 31-50, 51-100, >100   (the old 31-100 bin is split)
qsep_b0..b5  : per aligned query residue, hist over ALL 20A neighbours, mean over aligned cols.
qfrac_b0..b5 : per |i-j| bin, fraction of an aligned residue's contacts whose partner is ALSO
               aligned = (aligned-aligned endpoints) / (ALL endpoints), pooled over aligned cols.
               Decoupled from qsep by construction (it is the qaln/qsep ratio at the count level);
               low long-range qfrac flags a trivial-SSE alignment that leaks contacts outside span.
The topology-delta features (tmean/tamax/tstd) stay at the original 5 bins, derived from the
6-bin per-residue hist by merging cols 3+4 (exact) -> v17e tamax byte-identical to before.

Rationale (phase5h): qfrac6 gives a clean monotonic FAM>SFAM>FOLD>CROSS ordering at long range and
~3-4x the within-lddt conditional gain of qaln6, while staying low-correlation with qsep6 (qaln6
re-encoded qsep at long range, corr up to 0.93). qfrac6 replaces qaln6 in the v20 feature set.

Output cols = phase11b + tmean5 + tamax5 + tstd5 + qsep_b0..b5 + qfrac_b0..b5.
v20 = BASE22 + tamax5 + qsep6 + qfrac6 = 39 feats.
"""
from __future__ import annotations
import multiprocessing as mp
import os, re, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd
import pyarrow as pa, pyarrow.parquet as pq

TSV   = "/cbscratch/lolalignBenchmark/scope/allvsall_lol_training.tsv"
SCOP  = "/cbscratch/lolalignBenchmark/scope/data/scop_lookup.fix.tsv"
ADJ   = "/cbscratch/lolalignBenchmark/nn_classifier/outputs/domain_adjacency_20A.npz"
OUT   = "/cbscratch/lolalignBenchmark/nn_classifier/outputs/features_scalar_compact_topo_qsep6_qfrac6_parquet"
CHUNK = 1_000_000
SHARD_ROWS = 5_000_000
N_WORKERS  = int(os.environ.get("SLURM_CPUS_PER_TASK", 16))
BATCH = 5000

C = 5.0
EDGES = np.array([0.0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0001])
NB = len(EDGES) - 1   # 7  (ANF)
Q_ANF = [f"q_anf_b{i}" for i in range(NB)]
T_ANF = [f"t_anf_b{i}" for i in range(NB)]

# finer 6-bin seq-sep scheme for the per-residue histogram
SEP_EDGES6 = np.array([5, 13, 31, 51, 101])   # -> 1-4,5-12,13-30,31-50,51-100,>100
NSEP6 = 6
MERGE_LO, MERGE_HI = 3, 4   # 6-bin cols 3,4 merge back to original "31-100" 5-bin col 3
NSEP5 = 5
TMEAN = [f"tmean_b{i}" for i in range(NSEP5)]
TAMAX = [f"tamax_b{i}" for i in range(NSEP5)]
TSTD  = [f"tstd_b{i}"  for i in range(NSEP5)]
TOPO  = TMEAN + TAMAX + TSTD   # 15
QSEP6  = [f"qsep_b{i}" for i in range(NSEP6)]    # 6-bin all-neighbour
QFRAC6 = [f"qfrac_b{i}" for i in range(NSEP6)]   # 6-bin aligned-internal contact fraction

TSV_COLS = ["query","target","evalue","bits","lddt","alnlen","gapopen",
            "fident","qtmscore","qstart","qend","tstart","tend","cigar"]
KEEP = (["query","target","rel","q_scop","t_scop",
         "evalue","log_evalue","bits","log_bits","lddt","alnlen","log_alnlen",
         "gap_frac","gap_frac_real","fident","qtmscore"] + Q_ANF + T_ANF + TOPO + QSEP6 + QFRAC6)

_CIG = re.compile(r"(\d+)([MID])")
_ADJ_S: dict = {}; _ADJ_D: dict = {}; _ADJ_G: dict = {}
_RES_HIST6: dict = {}   # name -> (n_res x 6) float32 normalised per-residue seq-sep hist
_SEP_BIN6: dict = {}    # name -> int64 per-edge seq-sep bin (0..5)


def relation(qs, ts):
    shared = 0
    for a, b in zip(qs.split("."), ts.split(".")):
        if a == b: shared += 1
        else: break
    if shared >= 4: return "FAM"
    if shared >= 3: return "SFAM"
    if shared >= 2: return "FOLD"
    return "CROSS"


def aligned_indices(cigar, qstart, tstart):
    q_pos, t_pos = qstart-1, tstart-1; q_out, t_out = [], []
    for n_str, op in _CIG.findall(cigar):
        n = int(n_str)
        if op == "M":
            q_out.extend(range(q_pos, q_pos+n)); t_out.extend(range(t_pos, t_pos+n)); q_pos += n; t_pos += n
        elif op == "I": q_pos += n
        elif op == "D": t_pos += n
    return np.asarray(q_out, np.int64), np.asarray(t_out, np.int64)


def build_res_hist6(name):
    """Per-residue normalised seq-sep histogram (n_res x 6)."""
    s = _ADJ_S[name]; d = _ADJ_D[name]
    n = _ADJ_G[name].shape[0]
    cnt = np.zeros((n, NSEP6), np.float64)
    if s.size:
        np.add.at(cnt, (s, _SEP_BIN6[name]), 1.0)
        np.add.at(cnt, (d, _SEP_BIN6[name]), 1.0)
    tot = cnt.sum(1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        h = np.where(tot > 0, cnt / tot, 0.0)
    return h.astype(np.float32)


def _to5(h6):
    """Merge 6-bin hist cols 3+4 -> exact original 5-bin hist."""
    h5 = np.empty((h6.shape[0], NSEP5), h6.dtype)
    h5[:, 0] = h6[:, 0]; h5[:, 1] = h6[:, 1]; h5[:, 2] = h6[:, 2]
    h5[:, 3] = h6[:, MERGE_LO] + h6[:, MERGE_HI]
    h5[:, 4] = h6[:, 5]
    return h5


def topo_feat(qn, tn, qi, ti):
    """15-d (mean5 + absmax5 + std5) topology-delta on the 5-bin scheme."""
    h6q = _RES_HIST6.get(qn); h6t = _RES_HIST6.get(tn)
    if h6q is None or h6t is None:
        return np.zeros(3 * NSEP5, np.float32)
    m = (qi >= 0) & (qi < h6q.shape[0]) & (ti >= 0) & (ti < h6t.shape[0])
    qi, ti = qi[m], ti[m]
    if qi.size == 0:
        return np.zeros(3 * NSEP5, np.float32)
    delta = _to5(h6q[qi]) - _to5(h6t[ti])
    out = np.empty(3 * NSEP5, np.float32)
    out[:NSEP5] = delta.mean(0)
    out[NSEP5:2*NSEP5] = np.abs(delta).max(0)
    out[2*NSEP5:] = delta.std(0)
    return out


def qsep_feat(qn, qi):
    """6-d all-neighbour query |i-j| hist, mean over aligned query columns."""
    h = _RES_HIST6.get(qn)
    if h is None:
        return np.zeros(NSEP6, np.float32)
    qi = qi[(qi >= 0) & (qi < h.shape[0])]
    if qi.size == 0:
        return np.zeros(NSEP6, np.float32)
    return h[qi].mean(0).astype(np.float32)


def qfrac_feat(qn, qi):
    """6-d aligned-internal contact fraction: per |i-j| bin, pooled over aligned query residues,
    (aligned-aligned contact endpoints) / (ALL contact endpoints). The qaln/qsep ratio at the
    count level -> decoupled from the raw qsep shape."""
    s = _ADJ_S.get(qn)
    if s is None:
        return np.zeros(NSEP6, np.float32)
    d = _ADJ_D[qn]; b = _SEP_BIN6[qn]; n = _ADJ_G[qn].shape[0]
    qi = qi[(qi >= 0) & (qi < n)]
    if qi.size == 0:
        return np.zeros(NSEP6, np.float32)
    mask = np.zeros(n, bool); mask[qi] = True
    cs = mask[s].astype(np.float64); cd = mask[d].astype(np.float64)
    keep = (mask[s] & mask[d]).astype(np.float64)
    all_cnt = np.bincount(b, weights=cs + cd, minlength=NSEP6)[:NSEP6]
    aln_cnt = np.bincount(b, weights=keep * 2.0, minlength=NSEP6)[:NSEP6]
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(all_cnt > 0, aln_cnt / all_cnt, 0.0).astype(np.float32)


def anf_hist(name, aln_idx):
    g = _ADJ_G.get(name)
    if g is None: return np.zeros(NB, np.float32)
    N = g.shape[0]
    aln_idx = aln_idx[(aln_idx >= 0) & (aln_idx < N)]
    if aln_idx.size == 0: return np.zeros(NB, np.float32)
    s = _ADJ_S[name]; d = _ADJ_D[name]
    amask = np.zeros(N, bool); amask[aln_idx] = True
    if s.size == 0: return np.zeros(NB, np.float32)
    n_aln = (np.bincount(s, weights=amask[d].astype(np.float64), minlength=N)
             + np.bincount(d, weights=amask[s].astype(np.float64), minlength=N))
    sel = amask & (g > 0)
    if not sel.any(): return np.zeros(NB, np.float32)
    frac = n_aln[sel] / (g[sel].astype(np.float64) + C)
    h = np.histogram(frac, bins=EDGES)[0].astype(np.float32)
    tot = h.sum()
    return h / tot if tot > 0 else h


def gap_positions(cigar):
    g = 0
    for n, op in _CIG.findall(cigar):
        if op != "M": g += int(n)
    return g


def worker(batch):
    n = len(batch)
    gaps = np.zeros(n, np.float32)
    qh = np.zeros((n, NB), np.float32); th = np.zeros((n, NB), np.float32)
    tp = np.zeros((n, 3 * NSEP5), np.float32)
    qsp = np.zeros((n, NSEP6), np.float32)
    qfr = np.zeros((n, NSEP6), np.float32)
    for i, (qn, tn, cg, qs, ts) in enumerate(batch):
        gaps[i] = gap_positions(cg)
        qi, ti = aligned_indices(cg, int(qs), int(ts))
        qh[i] = anf_hist(qn, qi); th[i] = anf_hist(tn, ti)
        tp[i] = topo_feat(qn, tn, qi, ti)
        qsp[i] = qsep_feat(qn, qi)
        qfr[i] = qfrac_feat(qn, qi)
    return gaps, qh, th, tp, qsp, qfr


def main():
    global _ADJ_S, _ADJ_D, _ADJ_G, _RES_HIST6, _SEP_BIN6
    os.makedirs(OUT, exist_ok=True); t0 = time.time()
    scop = pd.read_csv(SCOP, sep="\t", header=None, names=["d","s"])
    d2s = dict(zip(scop["d"], scop["s"]))
    print(f"{len(d2s):,} SCOP domains", flush=True)

    print("loading adjacency npz ...", flush=True)
    npz = np.load(ADJ)
    names = set(k[:-2] for k in npz.files if k.endswith("_g"))
    for nm in names:
        _ADJ_S[nm] = npz[f"{nm}_s"].astype(np.int64)
        _ADJ_D[nm] = npz[f"{nm}_d"].astype(np.int64)
        _ADJ_G[nm] = npz[f"{nm}_g"]
        _SEP_BIN6[nm] = np.digitize(np.abs(_ADJ_S[nm] - _ADJ_D[nm]), SEP_EDGES6).astype(np.int64)
    print(f"  adjacency for {len(_ADJ_G):,} domains (C={C}, anf_bins={NB})", flush=True)

    print("building per-residue 6-bin seq-sep histograms ...", flush=True)
    for nm in names:
        _RES_HIST6[nm] = build_res_hist6(nm)
    print(f"  seq-sep hist for {len(_RES_HIST6):,} domains (sep_bins={NSEP6})", flush=True)

    ctx = mp.get_context("fork")
    pool = ProcessPoolExecutor(max_workers=N_WORKERS, mp_context=ctx)

    reader = pd.read_csv(TSV, sep="\t", header=None, names=TSV_COLS, usecols=TSV_COLS,
                         dtype={"query":"string","target":"string",
                                "evalue":"float32","bits":"float32","lddt":"float32",
                                "alnlen":"int32","gapopen":"int32",
                                "fident":"float32","qtmscore":"float32",
                                "qstart":"int32","tstart":"int32","cigar":"string"},
                         chunksize=CHUNK, na_filter=False)
    shard_idx, buf, buf_rows, n_in, n_kept = 0, [], 0, 0, 0
    rel_counts = {"FAM":0,"SFAM":0,"FOLD":0,"CROSS":0}
    for ci, chunk in enumerate(reader):
        n_in += len(chunk)
        chunk = chunk[chunk["query"] != chunk["target"]]
        chunk["q_scop"] = chunk["query"].map(d2s); chunk["t_scop"] = chunk["target"].map(d2s)
        chunk = chunk.dropna(subset=["q_scop","t_scop"])
        if len(chunk) == 0: continue
        chunk["rel"] = [relation(q,t) for q,t in zip(chunk["q_scop"], chunk["t_scop"])]

        qn = chunk["query"].to_numpy(); tn = chunk["target"].to_numpy()
        cg = chunk["cigar"].to_numpy(); qs = chunk["qstart"].to_numpy(); ts = chunk["tstart"].to_numpy()
        nrows = len(chunk)
        batches = [list(zip(qn[i:i+BATCH], tn[i:i+BATCH], cg[i:i+BATCH].tolist(),
                            qs[i:i+BATCH], ts[i:i+BATCH]))
                   for i in range(0, nrows, BATCH)]
        res = list(pool.map(worker, batches))
        num_gaps = np.concatenate([r[0] for r in res])
        q_hist = np.concatenate([r[1] for r in res]); t_hist = np.concatenate([r[2] for r in res])
        topo   = np.concatenate([r[3] for r in res])
        qsep   = np.concatenate([r[4] for r in res])
        qfrac  = np.concatenate([r[5] for r in res])

        alnlen = chunk["alnlen"].to_numpy(np.float32); gapopen = chunk["gapopen"].to_numpy(np.float32)
        ev = chunk["evalue"].to_numpy(np.float64); bt = chunk["bits"].to_numpy(np.float64)
        derived = pd.DataFrame({
            "log_evalue":   np.log(np.maximum(ev, 1e-300)).astype(np.float32),
            "log_bits":     (np.sign(bt)*np.log1p(np.abs(bt))).astype(np.float32),
            "log_alnlen":   np.log(alnlen + 1.0).astype(np.float32),
            "gap_frac":     (gapopen/np.maximum(alnlen,1.0)).astype(np.float32),
            "gap_frac_real":(num_gaps/np.maximum(alnlen,1.0)).astype(np.float32),
        }, index=chunk.index)
        anf = pd.DataFrame(np.concatenate([q_hist, t_hist], axis=1),
                           columns=Q_ANF + T_ANF, index=chunk.index)
        tdf = pd.DataFrame(topo, columns=TOPO, index=chunk.index)
        qsdf = pd.DataFrame(qsep, columns=QSEP6, index=chunk.index)
        qfdf = pd.DataFrame(qfrac, columns=QFRAC6, index=chunk.index)
        out = pd.concat([chunk[["query","target","rel","q_scop","t_scop",
                                "evalue","bits","lddt","alnlen","fident","qtmscore"]],
                         derived, anf, tdf, qsdf, qfdf], axis=1)[KEEP]
        for r,c in out["rel"].value_counts().items():
            rel_counts[r] = rel_counts.get(r,0)+int(c)
        buf.append(out); buf_rows += len(out); n_kept += len(out)
        if (ci+1) % 10 == 0:
            print(f"  chunk {ci}: in={n_in:,} kept={n_kept:,} rels={rel_counts} "
                  f"{time.time()-t0:.0f}s", flush=True)
        if buf_rows >= SHARD_ROWS:
            sh = pd.concat(buf, ignore_index=True)
            p = os.path.join(OUT, f"shard_{shard_idx:04d}.parquet")
            pq.write_table(pa.Table.from_pandas(sh, preserve_index=False), p, compression="zstd")
            print(f"  >> wrote {p} rows={len(sh):,} ({time.time()-t0:.0f}s)", flush=True)
            shard_idx += 1; buf, buf_rows = [], 0
    if buf:
        sh = pd.concat(buf, ignore_index=True)
        p = os.path.join(OUT, f"shard_{shard_idx:04d}.parquet")
        pq.write_table(pa.Table.from_pandas(sh, preserve_index=False), p, compression="zstd")
        shard_idx += 1
    pool.shutdown()
    print(f"\nDone. {shard_idx} shards, kept {n_kept:,}/{n_in:,} in {time.time()-t0:.0f}s", flush=True)
    print(f"rels: {rel_counts}", flush=True)


if __name__ == "__main__":
    main()
