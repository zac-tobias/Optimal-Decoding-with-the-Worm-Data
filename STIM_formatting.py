import stim
import numpy as np
import pymatching
import re
import ast
import stimcircuits
import julia
import subprocess
import json
import os
import csv
import random



def Adj_to_a_surface_code(Adj_, n_nodes, n_logicals):
    Adj = [[] for _ in range(n_nodes + 1)]
    Adj_weights = [[] for _ in range(n_nodes + 1)]
    B_nodes = []
    for a in Adj_:
        p = a[1]
        indices = [int(n[1:]) for n in a[0] if n.startswith('D')]
        if len(indices) == 2:
            Adj[indices[0]].append(indices[1] + 1)
            Adj[indices[1]].append(indices[0] + 1)

            Adj_weights[indices[1]].append(p)
            Adj_weights[indices[0]].append(p)

        if len(indices) == 1:
            Adj[indices[0]].append(n_nodes + 1)
            Adj[n_nodes].append(indices[0] + 1)

            Adj_weights[indices[0]].append(p)
            Adj_weights[n_nodes].append(p)

    return (Adj, Adj_weights)


def dem_logical_edges_surface_code(dem: stim.DetectorErrorModel, n_nodes):
    logical_edges = []
    for instruction in dem:
        if instruction.type != "error":
            continue

        terms = str(instruction).split(")")
        nodes = terms[1].strip().split()

        if "L0" in nodes:
            detectors = [int(n[1:]) for n in nodes if n.startswith("D")]
            if len(detectors) == 2:
                logical_edges.append([detectors[0] + 1, detectors[1] + 1])
            if len(detectors) == 1:
                logical_edges.append([detectors[0] + 1, n_nodes + 1])

    return logical_edges


def syndromes_and_observables_formatting(detection_events, observable_flips, num_shots, N_nodes):
    syndromes_list = []
    for i in range(num_shots):
        de_i = detection_events[i].tolist()
        de_i.insert(0, False)
        indices_i = np.flatnonzero(de_i)
        ind_i = indices_i.tolist()
        if len(ind_i) % 2 == 1:
            ind_i.append(N_nodes + 1)
        syndromes_list.append(ind_i)

    obs_flip_binary_ = np.where(observable_flips, -1, 1)
    obs_flip_binary = obs_flip_binary_.tolist()
    return (syndromes_list, obs_flip_binary)


def count_detectors_and_observables(dem: stim.DetectorErrorModel):
    max_detector = -1
    max_observable = -1

    for instruction in dem:
        if instruction.type != 'error':
            continue

        s = str(instruction)

        # Match detectors (D<number>) and logicals (L<number>)
        detector_matches = re.findall(r'D(\d+)', s)

        observable_matches = re.findall(r'L(\d+)', s)

        if detector_matches:
            max_detector = max(max_detector, max(map(int, detector_matches)))
        if observable_matches:
            max_observable = max(max_observable, max(map(int, observable_matches)))

    return max_detector + 1, max_observable + 1  # account for 0-indexing


def dem_to_Adj(dem: stim.DetectorErrorModel):
    Adj = []

    for instruction in dem:
        if instruction.type != 'error':
            continue

        s = str(instruction)  # e.g. "error(0.00666667) D13 D14"

        # Extract probability and node list
        head, tail = s.split(")", maxsplit=1)
        p_str = head.split("(")[1]
        p = float(p_str)

        nodes = tail.strip().split()
        Adj.append((nodes, p / (1 - p)))
    return Adj


def Adj_to_a(Adj_, n_nodes, n_logicals):
    Adj = [[] for _ in range(n_nodes)]
    Adj_weights = [[] for _ in range(n_nodes)]
    B_nodes = []
    for a in Adj_:
        # print(a[0])
        p = a[1]
        # print(a)
        indices = [int(n[1:]) for n in a[0] if n.startswith('D')]

        if len(indices) == 2:
            Adj[indices[0]].append(indices[1] + 1)
            Adj[indices[1]].append(indices[0] + 1)

            Adj_weights[indices[1]].append(p)
            Adj_weights[indices[0]].append(p)

    return (Adj, Adj_weights)


def dem_logical_edges(dem: stim.DetectorErrorModel):
    logical_edges = []
    for instruction in dem:
        if instruction.type != "error":
            continue

        terms = str(instruction).split(")")
        nodes = terms[1].strip().split()

        if "L0" in nodes:
            detectors = [int(n[1:]) for n in nodes if n.startswith("D")]
            if len(detectors) == 2:
                logical_edges.append([detectors[0] + 1, detectors[1] + 1])

    return logical_edges


def count_logical_errors_MWPM(circuit: stim.Circuit, detection_events, detector_error_model):
    # Configure a decoder using the circuit.
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

    # Run the decoder.
    predictions = matcher.decode_batch(detection_events)
    predictions_ = []
    for pred in predictions:
        if pred[0] == 0:
            predictions_.append([1])
        if pred[0] == 1:
            predictions_.append([-1])

    return predictions_



def prepare_shared_file(shared_filename, logicals_, adj, adj_w, t_auto, t_therm, N_samples):
    """Create a shared JSON file if it doesn't exist yet."""
    if os.path.exists(shared_filename):
        print(f"[INFO] Shared file already exists: {shared_filename}")
        return shared_filename

    data = {
        "logicals": logicals_,
        "adj": adj,
        "adj_w": adj_w,
        "t_autocorrelation": t_auto,
        "t_thermalize": t_therm,
        "N_samples": N_samples
    }

    with open(shared_filename, "w") as f:
        json.dump(data, f, indent=2)

    print(f"[INFO] Created shared file: {shared_filename}")
    return shared_filename



L_dist = 23
t_auto = 800
t_therm = 5000
N_samples = 5000
num_shots = 5

#p_pheno_list =  [0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425, 0.045]
#p_pheno_list =  [0.025, 0.0275, 0.03, 0.0325, 0.035]
p_pheno_list =  [0.0225]


for p_pheno in p_pheno_list:
    circuit_SC = stimcircuits.generate_circuit(
        "surface_code:rotated_memory_z",
        rounds=L_dist,
        distance=L_dist,
        after_clifford_depolarization=0.00,
        after_reset_flip_probability=0.00,
        before_measure_flip_probability=p_pheno,
        before_round_data_depolarization=p_pheno * 1.5,
        exclude_other_basis_detectors=True)

    dem = circuit_SC.detector_error_model(decompose_errors=True)
    n_detectors, n_observables = count_detectors_and_observables(dem.flattened())
    Adj_ = dem_to_Adj(dem.flattened())
    adj, adj_w = Adj_to_a_surface_code(Adj_, n_detectors, n_observables)
    logicals_ = dem_logical_edges_surface_code(dem.flattened(), n_detectors)
    sampler_SC = circuit_SC.compile_detector_sampler()

    detection_events, obs_flips = sampler_SC.sample(num_shots, separate_observables=True)
    synd_l, obs_ = syndromes_and_observables_formatting(detection_events, obs_flips, num_shots, n_detectors)

    # Generate one random number for both files
    random_num = random.randint(1, 10 ** 14)
    filename_DEM = f"DEM_data_rotated_surface_code_ME__L_{L_dist}_p_{p_pheno}.json"
    filename_synd = f"synd_data_rotated_surface_code_ME_L_{L_dist}_p_{p_pheno}_{random_num}.json"

    # --- DEM file (shared/static data) ---
    # ✅ If it already exists, that’s good — reuse it.
    if os.path.isfile(filename_DEM):
        print(f"[INFO] DEM file already exists, reusing: {filename_DEM}")
    else:
        data_full = {
            "logicals": logicals_,
            "adj": adj,
            "adj_w": adj_w,
            "t_autocorrelation": t_auto,
            "t_thermalize": t_therm,
            "N_samples": N_samples
        }
        with open(filename_DEM, "w") as f:
            json.dump(data_full, f, indent=2)
        print(f"[INFO] Created new DEM file: {filename_DEM}")

    # --- Syndrome file (always new) ---
    if os.path.isfile(filename_synd):
        filename_synd = f"error.csv"

    data_synd = {"syndromes": synd_l}
    with open(filename_synd, "w") as f:
        json.dump(data_synd, f, indent=2)


    # --- Run Julia ---
    result = subprocess.run(
        [
            "julia",
            "--project=/home/users/ztobias/.julia/environments/v1.11",
            "/home/users/ztobias/ML_decoder_new.jl",
            os.path.abspath(filename_DEM),
            os.path.abspath(filename_synd)
        ],
        capture_output=True,
        text=True,
    )

    #os.remove(filename_synd)

    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
    print("Return code:", result.returncode)

    results_MWPM = count_logical_errors_MWPM(circuit_SC, detection_events, dem)

    results = result.stdout

    lines = results.strip().split("\n")
    lists = []
    for line in lines:
        if "Any" in line:
            # Remove "Any"
            cleaned = line.replace("Any", "")
            # Evaluate the string to a Python list
            lst = ast.literal_eval(cleaned)
            lists.append(lst)


    corrections = lists[0]
    corrections_probs = lists[1]

    success_rate_MLD = 0

    data_list = []
    for j in range(num_shots):
        if corrections[j] == obs_[j][0]:
            MLD_status = 1
        else:
            MLD_status = 0

        if results_MWPM[j][0] == obs_[j][0]:
            MWPM_status = 1

        else:
            MWPM_status = 0

        data_list.append([MWPM_status, MLD_status, corrections_probs[j]])



    folder = "/home/users/ztobias/Decoding_Project_Data/Data_for_paper/RSC_ME_prescision_data"

    # Create the folder if it doesn't exist
    os.makedirs(folder, exist_ok=True)

    # Build the filename
    filename_only = f"data_prec_RSC_ME_L_{L_dist}_p_{str(p_pheno).replace('.', '_')}.csv"

    # Full path to the file
    filename = os.path.join(folder, filename_only)

    # Check if the file exists
    file_exists = os.path.isfile(filename)

    # Open the file in append mode
    with open(filename, mode="a", newline="") as f:
        writer = csv.writer(f)

        # Write header if file didn't exist
        if not file_exists:
            writer.writerow(["MWPM_status", "MLD_status", "Z_ratio"])

        # Write the data rows
        writer.writerows(data_list)

print("-----------------DONE-----------------------------")



