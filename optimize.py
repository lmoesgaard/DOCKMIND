import argparse
import json
import pandas as pd
import numpy as np
from deap import base, creator, tools, algorithms
from multiprocessing import Pool
from model import HitRateModel

def load_param_ranges(path):
    """Load parameter ranges from a JSON file."""
    with open(path, "r") as f:
        param_ranges = json.load(f)
    return param_ranges

def _individual_to_params(individual, param_ranges):
    """Translate an individual (list of parameters) into a dictionary."""
    return {key: individual[i] for i, key in enumerate(param_ranges.keys())}

def evaluate(individual, data, param_ranges):
    """Evaluate the fitness of an individual by calculating the mean squared error of predicted hit rates."""
    params = _individual_to_params(individual, param_ranges)
    model = HitRateModel(params)
    data["predicted_hit_rate"] = model.get_pprop_hr(pprop=data.pprop, definition=data.definition)
    return ((data.hit_rate_mean - data.predicted_hit_rate).abs() ** 2).mean(),

def check_bounds(mins, maxs):
    """Decorator to ensure that parameters remain within specified bounds after crossover and mutation."""
    def decorator(func):
        def wrapper(*args, **kwargs):
            offspring = func(*args, **kwargs)
            for child in offspring:
                for i in range(len(child)):
                    child[i] = max(min(child[i], maxs[i]), mins[i])
            return offspring
        return wrapper
    return decorator

def run(input_csv, output, param_ranges, generations, population_size, cxpb, mutpb, n_cpus):
    """Run the genetic algorithm to optimize model parameters."""
    # Create DEAP structures
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)
    
    toolbox = base.Toolbox()
    mins = [param_ranges[key][0] for key in param_ranges]
    maxs = [param_ranges[key][1] for key in param_ranges]

    def init_individual():
        return [np.random.uniform(mins[i], maxs[i]) for i in range(len(param_ranges))]

    toolbox.register("individual", tools.initIterate, creator.Individual, init_individual)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("mate", tools.cxBlend, alpha=0.5)
    toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.2)
    toolbox.register("select", tools.selTournament, tournsize=3)
    toolbox.decorate("mate", check_bounds(mins, maxs))
    toolbox.decorate("mutate", check_bounds(mins, maxs))

    # Define dictionary to store best parameters for each target
    best_params = {}
    # Read input CSV and group by target
    for target, data in pd.read_csv(input_csv).groupby("target"):
        toolbox.register("evaluate", evaluate, data=data, param_ranges=param_ranges)
        population = toolbox.population(n=population_size)

        # Use multiprocessing for parallel fitness evaluation
        with Pool(processes=n_cpus) as pool:
            toolbox.register("map", pool.map)
            algorithms.eaSimple(
                population, toolbox, cxpb=cxpb, mutpb=mutpb,
                ngen=generations, stats=None, halloffame=None, verbose=True
            )

        # Save the best individual for this target
        best = tools.selBest(population, k=1)[0]
        best_params[target] = _individual_to_params(best, param_ranges)
        best_params[target]["MSE"] = evaluate(best, data=data, param_ranges=param_ranges)[0]
        best_params[target]["target"] = target

    # Save the best parameters to the output JSON file
    with open(output, "w") as f:
        json.dump(best_params, f, indent=4)
    print(f"Results saved to {output}")

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Run genetic algorithm optimization.")
    parser.add_argument("--input_csv", required=True, help="Input CSV file path.")
    parser.add_argument("--output", default="output/fitted_params.json", help="Output JSON file path.")
    parser.add_argument("--param_ranges", default="config/param_ranges.json", help="Path to JSON file with parameter ranges.")
    parser.add_argument("--n_cpus", type=int, default=8, help="Number of CPUs to use for multiprocessing.")
    parser.add_argument("--generations", type=int, default=100)
    parser.add_argument("--population_size", type=int, default=5000)
    parser.add_argument("--cxpb", type=float, default=0.4)
    parser.add_argument("--mutpb", type=float, default=0.4)
    return parser.parse_args()

if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()
    # Load parameter ranges from JSON file
    param_ranges = load_param_ranges(args.param_ranges)
    # Run the genetic algorithm optimization
    run(
        args.input_csv,
        args.output,
        param_ranges,
        args.generations,
        args.population_size,
        args.cxpb,
        args.mutpb,
        args.n_cpus
    )
