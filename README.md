# Ancillary data and code for "Optimal Decoding with the Worm"

This repository contains simulation data and decoding code accompanying the paper:

> **Optimal Decoding with the Worm**  
> Zac Tobias et al. (2026)  
> [arXiv:2603.05428](https://arxiv.org/abs/2603.05428)

---

## Contents

- [`data/`](#data) — Logical error rate data used to produce the paper's figures
- [`Worm_decoder.jl`](#julia-code) — Worm MCMC decoder for matchable qLDPC codes
- [`worm_decoder_hyperbolic_sc.jl`](#julia-code) — Worm MCMC decoder for hyperbolic surface codes
- [`surface_code_correlated_decoding.jl`](#julia-code) — Correlated worm decoder for the RSC under depolarizing noise
- [`STIM_formatting_and_sampling.py`](#python-code) — Stim circuit generation and syndrome sampling

---

## Data

```
data/
├── hyperbolic_decoding_graphs/
├── Rotated_surface_code_measure...csv
├── correl_decoding_threshold_data...csv
├── correlated_decoding_comparisi...csv
├── correlated_decoding_parsing_da...csv
├── hyp_MWPM_comp_data.csv
└── hyperbolic_sc_threshold_data.csv
```

### CSV columns

| Column | Description |
|--------|-------------|
| `L` | Code distance (surface code) |
| `n` | Number of data qubits |
| `k` | Number of encoded logical qubits |
| `p` | Physical error rate |
| `N_samples` | Number of Monte Carlo samples |
| `N_errors` | Number of logical errors observed |
| `logical_error_rate` | Logical error rate estimate |
| `logical_error_rate_2sigma` | 2σ uncertainty on logical error rate |

The `hyperbolic_decoding_graphs/` subdirectory contains adjacency lists and logical
operator edge sets for each hyperbolic surface code instance, labelled by `[[n,k,d]]`
parameters.

---

## Julia code

All scripts use PyCall to access NetworkX for MWPM (via the Blossom algorithm).

### Requirements

Julia ≥ 1.9 with the following packages: `Random`, `LinearAlgebra`, `CSV`, `DataFrames`,
`Graphs`, `PyCall`. Python ≥ 3.9 with `networkx` installed.

### Files

| File | Description |
|------|-------------|
| `Worm_decoder.jl` | Worm MCMC decoder for matchable qLDPC codes. Takes a DEM JSON and syndrome JSON as command-line arguments; outputs per-syndrome corrections and MLD probability estimates to stdout. |
| `worm_decoder_hyperbolic_sc.jl` | Worm MCMC decoder for hyperbolic surface codes. Decoding graph and logical operators are specified directly as adjacency lists (see `data/hyperbolic_decoding_graphs/`). |
| `surface_code_correlated_decoding.jl` | Correlated worm decoder for the rotated surface code (RSC) under depolarizing noise, with iterative X↔Z parsing. Benchmarks against uncorrelated MWPM, correlated MWPM, and uncorrelated worm. |

### Running a simulation

Each script is self-contained. Simulation parameters (`L_list`, `p_list`, `N_real`,
`N_samp`, `t_auto`, etc.) are set as named constants near the top of each file.
Set `output_folder` to your desired output path before running.

**Worm decoder (surface code / generic qLDPC):**
```bash
julia Worm_decoder.jl path/to/dem.json path/to/syndromes.json
```

**Worm decoder (hyperbolic surface codes):**

Paste the appropriate `adj` and `loops_list` from `data/hyperbolic_decoding_graphs/`
into the designated section of `worm_decoder_hyperbolic_sc.jl`, then run:
```bash
julia worm_decoder_hyperbolic_sc.jl
```

**Correlated worm decoder (RSC, depolarizing noise):**
```bash
julia surface_code_correlated_decoding.jl
```

Output is one CSV file per `(L, p)`, written to `output_folder`.

---

## Python code

### Requirements

Python ≥ 3.9 with `stim` and `numpy`. Install via:
```bash
pip install stim numpy
```

### Files

| File | Description |
|------|-------------|
| `STIM_formatting_and_sampling.py` | Generates rotated surface code Stim circuits, samples syndromes, and calls `Worm_decoder.jl` via subprocess. Outputs logical error rates for a sweep over code distances and error rates. |

---

## Citation

```bibtex
@article{tobias2026worm,
  title   = {Optimal Decoding with the Worm},
  author  = {Tobias, Zac and others},
  journal = {arXiv preprint arXiv:2603.05428},
  year    = {2026},
  url     = {https://arxiv.org/abs/2603.05428}
}
```
