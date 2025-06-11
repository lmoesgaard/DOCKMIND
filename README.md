# DOCKMIND
DOCK Model from Inferred Normal Distributions. This repository provides tools for fitting a hit rate model to virtual screening data using genetic algorithm optimization. The model:
- Uses a bivariate normal model to account for correlation between docking score and pKi.
- Supports artifact modeling using an additional Gaussian distribution.

---

## Contents

- `model.py`: Implements the HitRateModel class for simulating hit rates.
- `optimize.py`: Uses a genetic algorithm to fit model parameters to hit rate data.

---

## Installation

This code requires Python 3.7+ and the following packages:

    pip install numpy pandas scipy deap

---

## Fitting

### 1. Prepare your input data

The experimental data for fitting the model should be supplied as a .csv file. An example is availble in `data/hitrates_conf_intervals.csv`. Ensure your input CSV file uses comma separation and has the following columns:

- `target`: identifier for each target (e.g., a protein)
- `pprop`: log-scaled proportion of top-ranked molecules (e.g., 3 means top 1 in 1000)
- `definition`: experimental threshold (e.g., pKi cutoff)
- `hit_rate_mean`: mean observed hit rate at that pProp

---

### 2. Define parameter search ranges

Create or edit the `param_ranges.json` JSON file like this:

    {
        "artifact_mean": [-5, -3],
        "artifact_std": [0.5, 1],
        "artifact_freq": [1e-7, 1e-2],
        "rho": [-0.9, -0.4],
        "exp_mean": [-4, 4],
        "exp_std": [0.5, 3]
    }
---

### 3. Run the optimization

    python optimize.py \
      --input_csv data/hitrates_conf_intervals.csv \
      --output output/fitted_params.json \
      --param_ranges config/param_ranges.json \
      --generations 100 \
      --population_size 5000 \
      --cxpb 0.4 \
      --mutpb 0.4 \
      --n_cpus -1

#### Arguments

| Argument           | Description                                                     | Default                     |
|--------------------|------------------------------------------------------------------|-----------------------------|
| --input_csv        | Path to the input CSV file with hit rate data                   | (Required)                  |
| --output           | Path to the output JSON file with fitted parameters             | `output/fitted_params.json`   |
| --param_ranges     | JSON file with parameter search ranges                          | `config/param_ranges.json`           |
| --generations      | Number of generations for genetic algorithm                     | `100`                           |
| --population_size  | Size of population for genetic algorithm                        | `5000`                        |
| --cxpb             | Crossover probability                                           | `0.4`                         |
| --mutpb            | Mutation probability                                            | `0.4`                         |
| --n_cpus           | Number of CPU cores to use for multiprocessing                  | `8`                           |

### Output

The final output JSON file contains best-fit parameters for each target, for example:

    {
        "target1": {
            "rho": -0.73,
            "exp_mean": 5.9,
            "exp_std": 1.1,
            "artifact_freq": 0.001,
            "artifact_mean": -4.5,
            "artifact_std": 0.8,
            "MSE": 0.0031,
            "target": "target1"
        }
    }

---

## Model Overview

`model.py` implements the HitRateModel class, which simulates virtual screening results based on the fitted parameters.

### Key methods

- `get_pprop_hr(pprop, definition):`  
  Estimates hit rate for a given top-ranked proportion (pprop) and hit definition threshold.

- `sample(pprop):`  
  Samples a pKi value for a molecule based on its rank position.


---

### Example: Using HitRateModel directly
The notebook, `demonstration.ipynb`, is supplemented to show how the HitRateModel functions can be utilized. Additionally, a simple example is available here:

    from model import HitRateModel

    params = {
        "rho": -0.7,
        "exp_mean": 6.0,
        "exp_std": 1.0,
        "artifact_freq": 0.001,
        "artifact_mean": -4.0,
        "artifact_std": 0.7
    }

    model = HitRateModel(params)
    pprop = 3  # Top 1 in 1000
    definition = 6.5  # Experimental pKi threshold

    hit_rate = model.get_pprop_hr(pprop, definition)
    print("Predicted hit rate:", hit_rate)

    sample = model.sample(pprop)
    print("pKi of sampled molecule:", sample)

---
