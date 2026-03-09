"""
STIM_formatting_and_sampling.py
--------------------------------------
Generates syndrome data for the rotated surface code with measurement errors
using Stim, decodes with both the Worm decoder (via Julia) and MWPM (via
PyMatching).

The Worm decoder is called as a subprocess running a Julia script
(ML_decoder_new.jl). Communication between Python and Julia is handled through
temporary JSON files containing the detector error model (DEM) and syndrome
data.

Usage:
    Run directly or import individual functions for use in other scripts.

Dependencies:
    stim, stimcircuits, pymatching, numpy, julia (subprocess)
"""

import ast
import csv
import json
import os
import random
import re
import subprocess

import numpy as np
import pymatching
import stim
import stimcircuits


# =============================================================================
# DEM / Graph Utilities
# =============================================================================

def count_detectors_and_observables(dem: stim.DetectorErrorModel) -> tuple[int, int]:
    """
    Count the number of detectors and logical observables in a DEM.

    Parameters
    ----------
    dem : stim.DetectorErrorModel

    Returns
    -------
    (n_detectors, n_observables) : tuple[int, int]
        Total counts, accounting for 0-indexing.
    """
    max_detector = -1
    max_observable = -1

    for instruction in dem:
        if instruction.type != "error":
            continue
        s = str(instruction)
        detector_matches = re.findall(r"D(\d+)", s)
        observable_matches = re.findall(r"L(\d+)", s)
        if detector_matches:
            max_detector = max(max_detector, max(map(int, detector_matches)))
        if observable_matches:
            max_observable = max(max_observable, max(map(int, observable_matches)))

    return max_detector + 1, max_observable + 1


def dem_to_Adj(dem: stim.DetectorErrorModel) -> list:
    """
    Convert a DEM to an adjacency list with log-odds weights.

    Each entry is a tuple (nodes, weight) where weight = p / (1 - p).

    Parameters
    ----------
    dem : stim.DetectorErrorModel

    Returns
    -------
    Adj : list of (list[str], float)
    """
    Adj = []
    for instruction in dem:
        if instruction.type != "error":
            continue
        s = str(instruction)
        head, tail = s.split(")", maxsplit=1)
        p = float(head.split("(")[1])
        nodes = tail.strip().split()
        Adj.append((nodes, p / (1 - p)))
    return Adj


def Adj_to_a_surface_code(Adj_: list, n_nodes: int, n_logicals: int) -> tuple:
    """
    Build adjacency lists for the surface code graph from a DEM adjacency list.

    Single-detector edges are connected to a virtual boundary node at index
    n_nodes + 1.

    Parameters
    ----------
    Adj_      : list of (nodes, weight) from dem_to_Adj
    n_nodes   : int, total number of detectors
    n_logicals: int, number of logical observables (unused, kept for API consistency)

    Returns
    -------
    (Adj, Adj_weights) : tuple of lists-of-lists
    """
    Adj = [[] for _ in range(n_nodes + 1)]
    Adj_weights = [[] for _ in range(n_nodes + 1)]

    for nodes, p in Adj_:
        indices = [int(n[1:]) for n in nodes if n.startswith("D")]

        if len(indices) == 2:
            # Bulk edge between two detectors
            Adj[indices[0]].append(indices[1] + 1)
            Adj[indices[1]].append(indices[0] + 1)
            Adj_weights[indices[0]].append(p)
            Adj_weights[indices[1]].append(p)

        elif len(indices) == 1:
            # Boundary edge: connect to virtual boundary node
            Adj[indices[0]].append(n_nodes + 1)
            Adj[n_nodes].append(indices[0] + 1)
            Adj_weights[indices[0]].append(p)
            Adj_weights[n_nodes].append(p)

    return Adj, Adj_weights



def dem_logical_edges_surface_code(dem: stim.DetectorErrorModel, n_nodes: int) -> list:
    """
    Extract logical edges from a DEM for the surface code.

    A logical edge is any error that flips observable L0. Single-detector
    logical edges are connected to the virtual boundary node at n_nodes + 1.

    Parameters
    ----------
    dem     : stim.DetectorErrorModel
    n_nodes : int, total number of detectors

    Returns
    -------
    logical_edges : list of [int, int]
    """
    logical_edges = []
    for instruction in dem:
        if instruction.type != "error":
            continue
        nodes = str(instruction).split(")")[1].strip().split()
        if "L0" not in nodes:
            continue
        detectors = [int(n[1:]) for n in nodes if n.startswith("D")]
        if len(detectors) == 2:
            logical_edges.append([detectors[0] + 1, detectors[1] + 1])
        elif len(detectors) == 1:
            logical_edges.append([detectors[0] + 1, n_nodes + 1])
    return logical_edges





# =============================================================================
# Syndrome Formatting
# =============================================================================

def syndromes_and_observables_formatting(
    detection_events: np.ndarray,
    observable_flips: np.ndarray,
    num_shots: int,
    N_nodes: int,
) -> tuple:
    """
    Convert raw Stim detection events into the sparse syndrome format expected
    by the Julia Worm decoder.

    Each syndrome is a list of detector indices (1-indexed) that fired. If an
    odd number of detectors fired, a virtual boundary node (N_nodes + 1) is
    appended to make the syndrome even-weight.

    Parameters
    ----------
    detection_events : np.ndarray, shape (num_shots, n_detectors)
    observable_flips : np.ndarray, shape (num_shots, n_observables)
    num_shots        : int
    N_nodes          : int, total number of detectors

    Returns
    -------
    (syndromes_list, obs_flip_binary) : tuple
        syndromes_list  : list of lists of int (1-indexed detector indices)
        obs_flip_binary : list of lists of int (+1 = no flip, -1 = flip)
    """
    syndromes_list = []
    for i in range(num_shots):
        # Prepend a False so detector indices become 1-indexed after flatnonzero
        de_i = [False] + detection_events[i].tolist()
        ind_i = np.flatnonzero(de_i).tolist()
        if len(ind_i) % 2 == 1:
            ind_i.append(N_nodes + 1)
        syndromes_list.append(ind_i)

    # Convert bool observable flips to ±1 (matching Julia convention)
    obs_flip_binary = np.where(observable_flips, -1, 1).tolist()
    return syndromes_list, obs_flip_binary


# =============================================================================
# MWPM Decoder
# =============================================================================

def count_logical_errors_MWPM(
    circuit: stim.Circuit,
    detection_events: np.ndarray,
    detector_error_model: stim.DetectorErrorModel,
) -> list:
    """
    Decode a batch of syndromes using PyMatching (MWPM) and return per-shot
    correction outcomes as ±1.

    Parameters
    ----------
    circuit              : stim.Circuit (unused, kept for API consistency)
    detection_events     : np.ndarray, shape (num_shots, n_detectors)
    detector_error_model : stim.DetectorErrorModel

    Returns
    -------
    predictions : list of [int]
        Each entry is [1] (no logical error predicted) or [-1] (logical error
        predicted).
    """
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)
    predictions = matcher.decode_batch(detection_events)
    return [[1] if pred[0] == 0 else [-1] for pred in predictions]


# =============================================================================
# Simulation Parameters
# =============================================================================

L_dist   = 23        # Code distance (also used as number of measurement rounds)
t_auto   = 10000     # Worm autocorrelation time (MCMC steps between samples)
t_therm  = 10000     # Worm thermalisation time (MCMC steps before sampling)
N_samples = 5000     # Number of Worm MCMC samples per syndrome
num_shots = 1        # Number of syndrome shots per p value

p_list = [0.025,0.03]    # Physical error rates to simulate

# Path to the Julia Worm decoder script and project environment.
# Update these to match your local setup before running.
JULIA_SCRIPT   = "path/to/worm_decoder.jl"
JULIA_PROJECT  = "path/to/.julia/environments/v1.11"

# =============================================================================
# Main Simulation Loop
# =============================================================================

for p in p_list:

    # ------------------------------------------------------------------
    # 1. Build Stim circuit (rotated surface code, phenomenological noise)
    # ------------------------------------------------------------------
    circuit_SC = stimcircuits.generate_circuit(
        "surface_code:rotated_memory_z",
        rounds=L_dist,
        distance=L_dist,
        after_clifford_depolarization=0.00,
        after_reset_flip_probability=0.00,
        before_measure_flip_probability=p,
        before_round_data_depolarization=p * 1.5,
        exclude_other_basis_detectors=True,
    )

    # ------------------------------------------------------------------
    # 2. Extract the detector error model and build the decoding graph
    # ------------------------------------------------------------------
    dem = circuit_SC.detector_error_model(decompose_errors=True)
    n_detectors, n_observables = count_detectors_and_observables(dem.flattened())

    Adj_  = dem_to_Adj(dem.flattened())
    adj, adj_w = Adj_to_a_surface_code(Adj_, n_detectors, n_observables)
    logicals_ = dem_logical_edges_surface_code(dem.flattened(), n_detectors)

    # ------------------------------------------------------------------
    # 3. Sample syndromes
    # ------------------------------------------------------------------
    sampler_SC = circuit_SC.compile_detector_sampler()
    detection_events, obs_flips = sampler_SC.sample(
        num_shots, separate_observables=True
    )
    synd_l, obs_ = syndromes_and_observables_formatting(
        detection_events, obs_flips, num_shots, n_detectors
    )

    # ------------------------------------------------------------------
    # 4. Write DEM file (shared across shots for the same p / L)
    #    and syndrome file (unique per batch of shots)
    # ------------------------------------------------------------------
    random_num   = random.randint(1, 10**14)
    filename_DEM  = f"DEM_data_rotated_surface_code_ME__L_{L_dist}_p_{p}.json"
    filename_synd = f"synd_data_rotated_surface_code_ME_L_{L_dist}_p_{p}_{random_num}.json"

    # DEM file is static for a given (L, p) — reuse if it already exists
    if os.path.isfile(filename_DEM):
        print(f"[INFO] Reusing existing DEM file: {filename_DEM}")
    else:
        dem_data = {
            "logicals":        logicals_,
            "adj":             adj,
            "adj_w":           adj_w,
            "t_autocorrelation": t_auto,
            "t_thermalize":    t_therm,
            "N_samples":       N_samples,
        }
        with open(filename_DEM, "w") as f:
            json.dump(dem_data, f, indent=2)
        print(f"[INFO] Created DEM file: {filename_DEM}")

    # Syndrome file is unique to each batch — fall back to error.csv on collision
    if os.path.isfile(filename_synd):
        filename_synd = "error.csv"

    with open(filename_synd, "w") as f:
        json.dump({"syndromes": synd_l}, f, indent=2)

    # ------------------------------------------------------------------
    # 5. Run the Julia Worm decoder as a subprocess
    # ------------------------------------------------------------------
    result = subprocess.run(
        [
            "julia",
            f"--project={JULIA_PROJECT}",
            JULIA_SCRIPT,
            os.path.abspath(filename_DEM),
            os.path.abspath(filename_synd),
        ],
        capture_output=True,
        text=True,
    )

    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
    print("Return code:", result.returncode)

    # ------------------------------------------------------------------
    # 6. Decode with MWPM (PyMatching) for comparison
    # ------------------------------------------------------------------
    results_MWPM = count_logical_errors_MWPM(circuit_SC, detection_events, dem)

    # ------------------------------------------------------------------
    # 7. Parse Julia output
    #    The Julia script prints two Julia Any[] arrays on separate lines:
    #      Line 1: per-shot correction decisions  (±1)
    #      Line 2: per-shot MLD probability estimates
    # ------------------------------------------------------------------
    lines = result.stdout.strip().split("\n")
    lists = []
    for line in lines:
        if "Any" in line:
            lst = ast.literal_eval(line.replace("Any", ""))
            lists.append(lst)

    corrections       = lists[0]   # Worm decoder correction decisions
    corrections_probs = lists[1]   # Worm decoder MLD probability estimates

    # ------------------------------------------------------------------
    # 8. Evaluate decoding success and collect results
    # ------------------------------------------------------------------
    data_list = []
    for j in range(num_shots):
        MLD_status  = 1 if corrections[j]      == obs_[j][0] else 0
        MWPM_status = 1 if results_MWPM[j][0]  == obs_[j][0] else 0
        data_list.append([MWPM_status, MLD_status, corrections_probs[j]])
