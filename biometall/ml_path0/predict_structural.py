# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
import logging
import math
from collections import namedtuple
from pathlib import Path

import joblib
import sys
import numpy as np


log = logging.getLogger(__name__)

_NAN = float("nan")

_SHELLS = [
    (0.0, 2.5),
    (2.5, 4.0),
    (4.0, 5.5),
]

_GEO_RADIUS = 3.5
_STRATEGY = "p95"


def _shell_names():
    names = []

    for i in range(3):
        lbl = f"s{i + 1}"

        names += [
            f"n_polar_{lbl}",
            f"n_apolar_{lbl}",
            f"n_chrg_pos_{lbl}",
            f"n_chrg_neg_{lbl}",
            f"n_aromatic_{lbl}",
            f"n_cys_met_{lbl}",
            f"n_his_{lbl}",
            f"n_N_donors_{lbl}",
            f"n_S_donors_{lbl}",
            f"n_O_sc_donors_{lbl}",
            f"n_O_bb_donors_{lbl}",
            f"n_residues_{lbl}",
            f"mean_dist_{lbl}",
            f"std_dist_{lbl}",
        ]

    return names


FEATURE_SCHEMA = (
    _shell_names()
    + [
        "coord_n_atoms",
        "coord_min_dist",
        "coord_mean_dist",
        "coord_std_dist",
        "coord_n_tight",
        "angle_mean",
        "angle_std",
        "angle_min",
        "coord_radius_of_gyration",
        "coord_span",
        "geom_linear",
        "geom_tetrahedral",
    ]
    + [
        "burial_proxy",
        "local_charge_s1",
        "frac_polar_s1",
        "frac_polar_s2",
        "mean_hydrophobicity_s1",
    ]
    + [
        "bm_probe_density_3A",
        "bm_probe_density_5A",
        "bm_is_cluster_center",
    ]
)


Atom = namedtuple(
    "Atom",
    [
        "record",
        "name",
        "resname",
        "chain",
        "resnum",
        "x",
        "y",
        "z",
        "element",
    ],
)


def _parse_pdb_atoms(pdb_path: Path, include_hetatm=False):
    atoms = []

    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM  "):
                record = "ATOM"
            elif include_hetatm and line.startswith("HETATM"):
                record = "HETATM"
            else:
                continue

            try:
                atoms.append(
                    Atom(
                        record=record,
                        name=line[12:16].strip(),
                        resname=line[17:20].strip(),
                        chain=line[21].strip() or "A",
                        resnum=int(line[22:26]),
                        x=float(line[30:38]),
                        y=float(line[38:46]),
                        z=float(line[46:54]),
                        element=(
                            line[76:78].strip()
                            if len(line) >= 78
                            else line[12].strip()
                        ),
                    )
                )
            except (ValueError, IndexError):
                continue

    return atoms


def _get_solution_center(pdb_path: Path):
    atoms = _parse_pdb_atoms(
        pdb_path,
        include_hetatm=True,
    )

    if not atoms:
        return None

    x = sum(a.x for a in atoms) / len(atoms)
    y = sum(a.y for a in atoms) / len(atoms)
    z = sum(a.z for a in atoms) / len(atoms)

    return float(x), float(y), float(z)


def _load_biometall_clusters(probes_dir: Path):
    clusters = []

    for pdb_file in sorted(probes_dir.glob("*.pdb")):
        center = _get_solution_center(pdb_file)

        if center is None:
            continue

        atoms = _parse_pdb_atoms(
            pdb_file,
            include_hetatm=True,
        )

        clusters.append(
            (
                len(atoms),
                center,
            )
        )

    clusters.sort(key=lambda x: -x[0])

    return clusters


_POLAR = frozenset({
    "SER", "THR", "ASN", "GLN"
})

_APOLAR = frozenset({
    "ALA", "VAL", "LEU", "ILE", "PRO", "MET"
})

_CHRG_POS = frozenset({
    "ARG", "LYS", "HIS"
})

_CHRG_NEG = frozenset({
    "ASP", "GLU"
})

_AROMATIC = frozenset({
    "PHE", "TYR", "TRP", "HIS"
})

_CYS_MET = frozenset({
    "CYS", "MET"
})

_HIS = frozenset({
    "HIS"
})

_N_DONOR = frozenset({
    "NE2", "NZ", "ND1"
})

_S_DONOR = frozenset({
    "SG", "SD"
})

_O_SC_DONOR = frozenset({
    "OD1", "OD2", "OE1", "OE2", "OH"
})

_BACKBONE = frozenset({
    "N", "CA", "C", "O", "OXT"
})

_KD = {
    "ILE": 4.5,
    "VAL": 4.2,
    "LEU": 3.8,
    "PHE": 2.8,
    "CYS": 2.5,
    "MET": 1.9,
    "ALA": 1.8,
    "GLY": -0.4,
    "THR": -0.7,
    "SER": -0.8,
    "TRP": -0.9,
    "TYR": -1.3,
    "PRO": -1.6,
    "HIS": -3.2,
    "GLU": -3.5,
    "GLN": -3.5,
    "ASP": -3.5,
    "ASN": -3.5,
    "LYS": -3.9,
    "ARG": -4.5,
}

_CHARGE = {
    "ASP": -1.0,
    "GLU": -1.0,
    "ARG": +1.0,
    "LYS": +1.0,
    "HIS": +0.5,
}

_POLAR_EL = frozenset({
    "N", "O"
})


def _shell_features(atoms, center, shells):
    cx, cy, cz = center
    n = len(shells)

    shell_r2 = [
        (r0 * r0, r1 * r1)
        for r0, r1 in shells
    ]

    polar = [0] * n
    apolar = [0] * n
    chrg_pos = [0] * n
    chrg_neg = [0] * n
    aromatic = [0] * n
    cys_met = [0] * n
    his = [0] * n

    n_don = [0] * n
    s_don = [0] * n
    osc_don = [0] * n
    obb_don = [0] * n

    residues = [set() for _ in range(n)]
    ca_dists = [[] for _ in range(n)]

    for a in atoms:
        dx = a.x - cx
        dy = a.y - cy
        dz = a.z - cz
        d2 = dx * dx + dy * dy + dz * dz

        si = -1

        for idx, (r0_2, r1_2) in enumerate(shell_r2):
            if r0_2 <= d2 < r1_2:
                si = idx
                break

        if si < 0:
            continue

        rn = a.resname
        an = a.name.strip()

        if rn in _POLAR:
            polar[si] += 1

        if rn in _APOLAR:
            apolar[si] += 1

        if rn in _CHRG_POS:
            chrg_pos[si] += 1

        if rn in _CHRG_NEG:
            chrg_neg[si] += 1

        if rn in _AROMATIC:
            aromatic[si] += 1

        if rn in _CYS_MET:
            cys_met[si] += 1

        if rn in _HIS:
            his[si] += 1

        if an in _N_DONOR:
            n_don[si] += 1
        elif an in _S_DONOR:
            s_don[si] += 1
        elif an in _O_SC_DONOR:
            osc_don[si] += 1
        elif an == "O" and an in _BACKBONE:
            obb_don[si] += 1

        if an == "CA":
            residues[si].add((a.chain, a.resnum))
            ca_dists[si].append(math.sqrt(d2))

    feats = {}

    for i in range(n):
        lbl = f"s{i + 1}"

        feats[f"n_polar_{lbl}"] = float(polar[i])
        feats[f"n_apolar_{lbl}"] = float(apolar[i])
        feats[f"n_chrg_pos_{lbl}"] = float(chrg_pos[i])
        feats[f"n_chrg_neg_{lbl}"] = float(chrg_neg[i])
        feats[f"n_aromatic_{lbl}"] = float(aromatic[i])
        feats[f"n_cys_met_{lbl}"] = float(cys_met[i])
        feats[f"n_his_{lbl}"] = float(his[i])
        feats[f"n_N_donors_{lbl}"] = float(n_don[i])
        feats[f"n_S_donors_{lbl}"] = float(s_don[i])
        feats[f"n_O_sc_donors_{lbl}"] = float(osc_don[i])
        feats[f"n_O_bb_donors_{lbl}"] = float(obb_don[i])
        feats[f"n_residues_{lbl}"] = float(len(residues[i]))

        dists = ca_dists[i]

        if dists:
            mu = sum(dists) / len(dists)
            std = math.sqrt(
                sum((d - mu) ** 2 for d in dists) / len(dists)
            )
        else:
            mu = 0.0
            std = 0.0

        feats[f"mean_dist_{lbl}"] = mu
        feats[f"std_dist_{lbl}"] = std

    return feats


def _geometry_features(atoms, center, radius):
    cx, cy, cz = center
    r2 = radius * radius
    donors = []

    for a in atoms:
        an = a.name.strip()

        if not (
            an in _N_DONOR
            or an in _S_DONOR
            or an in _O_SC_DONOR
            or (an == "O" and an in _BACKBONE)
        ):
            continue

        dx = a.x - cx
        dy = a.y - cy
        dz = a.z - cz
        d2 = dx * dx + dy * dy + dz * dz

        if d2 <= r2:
            donors.append(
                (
                    a.x,
                    a.y,
                    a.z,
                    math.sqrt(d2),
                )
            )

    n = len(donors)

    if n == 0:
        return {
            "coord_n_atoms": 0.0,
            "coord_min_dist": 0.0,
            "coord_mean_dist": 0.0,
            "coord_std_dist": 0.0,
            "coord_n_tight": 0.0,
            "angle_mean": _NAN,
            "angle_std": _NAN,
            "angle_min": _NAN,
            "coord_radius_of_gyration": 0.0,
            "coord_span": 0.0,
            "geom_linear": _NAN,
            "geom_tetrahedral": _NAN,
        }

    dists = [d for _, _, _, d in donors]

    min_d = min(dists)
    mean_d = sum(dists) / n

    std_d = (
        math.sqrt(
            sum((d - mean_d) ** 2 for d in dists) / n
        )
        if n > 1
        else 0.0
    )

    n_tight = sum(
        1
        for d in dists
        if 2.0 <= d <= 2.8
    )

    rg = math.sqrt(sum(d * d for d in dists) / n)

    span = 0.0

    if n > 1:
        for i in range(n):
            xi, yi, zi, _ = donors[i]

            for j in range(i + 1, n):
                xj, yj, zj, _ = donors[j]

                dx = xi - xj
                dy = yi - yj
                dz = zi - zj

                span = max(
                    span,
                    math.sqrt(dx * dx + dy * dy + dz * dz),
                )

    if n < 2:
        return {
            "coord_n_atoms": float(n),
            "coord_min_dist": min_d,
            "coord_mean_dist": mean_d,
            "coord_std_dist": std_d,
            "coord_n_tight": float(n_tight),
            "angle_mean": _NAN,
            "angle_std": _NAN,
            "angle_min": _NAN,
            "coord_radius_of_gyration": rg,
            "coord_span": span,
            "geom_linear": _NAN,
            "geom_tetrahedral": _NAN,
        }

    angles = []
    n_lin = 0
    n_tet = 0

    for i in range(n):
        xi, yi, zi, di = donors[i]

        if di < 1e-6:
            continue

        for j in range(i + 1, n):
            xj, yj, zj, dj = donors[j]

            if dj < 1e-6:
                continue

            ux = (xi - cx) / di
            uy = (yi - cy) / di
            uz = (zi - cz) / di

            vx = (xj - cx) / dj
            vy = (yj - cy) / dj
            vz = (zj - cz) / dj

            dot = max(
                -1.0,
                min(
                    1.0,
                    ux * vx + uy * vy + uz * vz,
                ),
            )

            ang = math.degrees(math.acos(dot))
            angles.append(ang)

            if abs(ang - 180.0) <= 20.0:
                n_lin += 1

            if abs(ang - 109.5) <= 15.0:
                n_tet += 1

    if angles:
        am = sum(angles) / len(angles)

        astd = math.sqrt(
            sum((a - am) ** 2 for a in angles) / len(angles)
        )

        amin = min(angles)
        flin = n_lin / len(angles)
        ftet = n_tet / len(angles)
    else:
        am = _NAN
        astd = _NAN
        amin = _NAN
        flin = _NAN
        ftet = _NAN

    return {
        "coord_n_atoms": float(n),
        "coord_min_dist": min_d,
        "coord_mean_dist": mean_d,
        "coord_std_dist": std_d,
        "coord_n_tight": float(n_tight),
        "angle_mean": am,
        "angle_std": astd,
        "angle_min": amin,
        "coord_radius_of_gyration": rg,
        "coord_span": span,
        "geom_linear": flin,
        "geom_tetrahedral": ftet,
    }


def _physicochemistry_features(atoms, center, shells):
    cx, cy, cz = center

    s1_r2 = shells[0][1] ** 2
    s2_r2 = (
        shells[1][1] ** 2
        if len(shells) > 1
        else (shells[0][1] * 2) ** 2
    )

    bur_r2 = 8.0 ** 2

    s1_res = {}
    s2_res = {}

    s1_n = 0
    s1_tot = 0

    s2_n = 0
    s2_tot = 0

    cb = 0

    for a in atoms:
        dx = a.x - cx
        dy = a.y - cy
        dz = a.z - cz
        d2 = dx * dx + dy * dy + dz * dz

        an = a.name.strip()

        if d2 <= bur_r2 and (
            an == "CB"
            or (an == "CA" and a.resname == "GLY")
        ):
            cb += 1

        if d2 < s1_r2:
            s1_res[(a.chain, a.resnum)] = a.resname
            s1_tot += 1

            if a.element.strip().upper() in _POLAR_EL:
                s1_n += 1

        if d2 < s2_r2:
            s2_res[(a.chain, a.resnum)] = a.resname
            s2_tot += 1

            if a.element.strip().upper() in _POLAR_EL:
                s2_n += 1

    charge = sum(
        _CHARGE.get(rn, 0.0)
        for rn in s1_res.values()
    )

    fp1 = s1_n / s1_tot if s1_tot > 0 else _NAN
    fp2 = s2_n / s2_tot if s2_tot > 0 else _NAN

    kd = [
        _KD.get(rn, 0.0)
        for rn in s1_res.values()
    ]

    mkd = sum(kd) / len(kd) if kd else _NAN

    return {
        "burial_proxy": float(cb),
        "local_charge_s1": charge,
        "frac_polar_s1": fp1,
        "frac_polar_s2": fp2,
        "mean_hydrophobicity_s1": mkd,
    }


def _external_features(center, bm_clusters, total_probes):
    if total_probes <= 0 or not bm_clusters:
        return {
            "bm_probe_density_3A": _NAN,
            "bm_probe_density_5A": _NAN,
            "bm_is_cluster_center": _NAN,
        }

    cx, cy, cz = center

    p3 = 0
    p5 = 0
    nearest = float("inf")

    for n_p, (bx, by, bz) in bm_clusters:
        dx = bx - cx
        dy = by - cy
        dz = bz - cz

        d = math.sqrt(dx * dx + dy * dy + dz * dz)

        if d <= 3.0:
            p3 += n_p

        if d <= 5.0:
            p5 += n_p

        nearest = min(nearest, d)

    return {
        "bm_probe_density_3A": p3 / total_probes,
        "bm_probe_density_5A": p5 / total_probes,
        "bm_is_cluster_center": 1.0 if nearest <= 1.5 else 0.0,
    }


def _extract_features(
    atoms,
    center,
    bm_clusters,
    bm_total,
):
    try:
        row = {
            **_shell_features(
                atoms,
                center,
                _SHELLS,
            ),
            **_geometry_features(
                atoms,
                center,
                _GEO_RADIUS,
            ),
            **_physicochemistry_features(
                atoms,
                center,
                _SHELLS,
            ),
            **_external_features(
                center,
                bm_clusters,
                bm_total,
            ),
        }

        return [
            float(row.get(k, _NAN))
            for k in FEATURE_SCHEMA
        ]

    except Exception:
        log.exception(
            "Feature extraction error at %s",
            center,
        )
        return None


def _aggregate(scores, tau):
    if not scores:
        return {
            "is_metalloprotein": False,
            "p_metal": 0.0,
            "n_sites_scored": 0,
        }

    arr = np.array(scores)

    p_metal = float(
        np.percentile(arr, 95)
    )

    return {
        "is_metalloprotein": bool(
            p_metal >= tau
        ),
        "p_metal": round(
            p_metal,
            6,
        ),
        "n_sites_scored": len(scores),
    }


def _deduplicate(
    centers,
    probs,
    thr=2.0,
):
    n = len(centers)
    keep = [True] * n

    for i in range(n):
        if not keep[i]:
            continue

        for j in range(i + 1, n):
            if not keep[j]:
                continue

            dx = centers[i][0] - centers[j][0]
            dy = centers[i][1] - centers[j][1]
            dz = centers[i][2] - centers[j][2]

            d = math.sqrt(
                dx * dx +
                dy * dy +
                dz * dz
            )

            if d < thr:
                if probs[i] >= probs[j]:
                    keep[j] = False
                else:
                    keep[i] = False
                    break

    return [
        i
        for i in range(n)
        if keep[i]
    ]


def predict_protein(
    pdb_path: Path,
    probes_dir: Path,
    model_path: Path,
    tau: float,
):
    biometall_dir = Path(__file__).resolve().parent.parent

    if str(biometall_dir) not in sys.path:
        sys.path.insert(
            0,
            str(biometall_dir),
        )

    model = joblib.load(model_path)

    atoms = _parse_pdb_atoms(pdb_path)

    if not atoms:
        return _empty_result(
            str(pdb_path),
            tau,
        )

    bm_clusters = _load_biometall_clusters(
        probes_dir
    )

    if not bm_clusters:
        return _empty_result(
            str(pdb_path),
            tau,
        )

    bm_total = sum(
        n
        for n, _ in bm_clusters
    )

    site_features = []
    site_centers = []

    for _, xyz in bm_clusters:
        center = (
            float(xyz[0]),
            float(xyz[1]),
            float(xyz[2]),
        )

        feats = _extract_features(
            atoms,
            center,
            bm_clusters,
            bm_total,
        )

        if feats is not None:
            site_features.append(feats)
            site_centers.append(center)

    if not site_features:
        return _empty_result(
            str(pdb_path),
            tau,
        )

    X = np.array(
        site_features,
        dtype="float32",
    )

    proba_all = model.predict_proba(X)
    probs = proba_all[:, 1]

    kept = _deduplicate(
        site_centers,
        probs.tolist(),
    )

    site_centers = [
        site_centers[i]
        for i in kept
    ]

    probs = probs[kept]

    agg = _aggregate(
        probs.tolist(),
        tau,
    )

    top_sites = sorted(
        [
            {
                "center_x": round(cx, 3),
                "center_y": round(cy, 3),
                "center_z": round(cz, 3),
                "p_metal": round(
                    float(probs[i]),
                    4,
                ),
            }
            for i, (cx, cy, cz)
            in enumerate(site_centers)
        ],
        key=lambda s: -s["p_metal"],
    )

    return {
        "pdb_path": str(pdb_path),
        "is_metalloprotein": agg["is_metalloprotein"],
        "p_metal": agg["p_metal"],
        "aggregation_strategy": _STRATEGY,
        "tau_protein": tau,
        "n_sites_scored": agg["n_sites_scored"],
        "top_sites": top_sites,
    }


def _empty_result(
    pdb_path,
    tau,
):
    return {
        "pdb_path": pdb_path,
        "is_metalloprotein": False,
        "p_metal": 0.0,
        "aggregation_strategy": _STRATEGY,
        "tau_protein": tau,
        "n_sites_scored": 0,
        "top_sites": [],
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "MetalPlacer Path 0 structural screener "
            "using precomputed BioMetAll sites."
        )
    )

    parser.add_argument(
        "--pdb",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--probes-dir",
        type=Path,
        required=True,
        help="Directory containing BioMetAll solution PDB files.",
    )

    parser.add_argument(
        "--binary-model",
        type=Path,
        required=True,
        dest="model",
    )

    parser.add_argument(
        "--tau",
        type=float,
        default=0.60,
    )

    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output JSON file.",
    )

    args = parser.parse_args()

    result = predict_protein(
        pdb_path=args.pdb,
        probes_dir=args.probes_dir,
        model_path=args.model,
        tau=args.tau,
    )

    with open(args.output, "w") as fh:
        json.dump(
            result,
            fh,
            indent=2,
        )


if __name__ == "__main__":
    main()