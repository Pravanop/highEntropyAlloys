import datetime
from itertools import combinations
from initial_config import InitialConfig
import numpy as np
import random
from tqdm import tqdm
import pickle
from collections import Counter
from nearest_neighbour import create_neighbor_list
import wandb


class MonteCarlo:

    def __init__(self, initial_conf, config_dict):

        self.initial_config = initial_conf
        self.lookup = initial_config.lookup
        self.config_dict = config_dict
        self.energy_trajectory = []
        self.structure_trajectory = []
        self.steps = []
        self.kb = 8.612e-5
        self.neighbour_list = create_neighbor_list(self.initial_config.final_bcclattice, flag = 2)

    def hamiltonian(self, point, arr):
        site = arr[point[0]]
        return sum([self.lookup[str(sorted([site, arr[i]]))] for i in point[1]])

    def hamiltonian2(self, point, arr, T):
        composition = [arr[i] for i in point[1]] + [arr[point[0]]]
        counter = dict(Counter(composition))
        total = sum(counter.values())
        counter = {key: round(value / total, 2) for key, value in counter.items()}
        combs = list(combinations(counter.keys(), 2))
        return sum([self.lookup[str(sorted([i[0], i[1]]))] * counter[i[0]] * counter[i[1]] for i in combs])
        # entropy = -self.kb * sum([np.log(i) * i for i in counter.values()])
        # return round(enthalpy - T * entropy, 3)

    def energy_finder(self, arr):
        return sum([self.hamiltonian(point, arr) for point in self.neighbour_list]) / 2

    def energy_finder1(self, arr, T):
        return sum([self.hamiltonian2(point, arr, T) for point in self.neighbour_list]) / 8

    def pair_swapper(self, arr):

        rand = random.choices(self.initial_config.non_zero, k=2)
        arr[rand[0][0], rand[0][1], rand[0][2]], arr[rand[1][0], rand[1][1], rand[1][2]] = (
            arr[rand[1][0], rand[1][1], rand[1][2]], arr[rand[0][0], rand[0][1], rand[0][2]])
        return arr

    def n_pair_swapper(self, n, arr):
        for i in range(n):
            arr = self.pair_swapper(arr)
        return arr

    def boltzmann_probability(self, delta_e, temperature):
        return np.exp(-delta_e / (self.kb * temperature))

    def logger(self):
        dump_dict = {
            'steps_trajectory': np.array(self.steps),
            'energy_trajectory': np.array(self.energy_trajectory),
            'structure_trajectory': np.array(self.structure_trajectory),
            'Config': self.config_dict,
        }

        now = datetime.datetime.now()
        now.isoformat()
        system = '-'.join(list(self.initial_config.ele_list))
        with open(f'./dump/{system}_{self.initial_config.atoms}_{self.config_dict["T"]}_{now}.pickle', 'wb') as handle:
            pickle.dump(dump_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def mc_single_temp(self, n_trails, temp, lattice, log):
        x_i = lattice.copy()
        if log["ham"] == "enthalpy":
            e_i = self.energy_finder1(x_i, temp)
        elif log["ham"] == "bonds":
            e_i = self.energy_finder(x_i)
        # self.steps = []
        count = 0
        swaps = 6
        wandb.log({"No of Swaps": swaps})
        for i in tqdm(range(n_trails), desc=f"Running at {temp} K with {swaps} swap"):

            x_iplus1 = self.n_pair_swapper(n=swaps, arr=x_i.copy())
            try:
                assert np.array_equal(x_iplus1, x_i) is not True
            except AssertionError:
                x_iplus1 = self.n_pair_swapper(n=1, arr=x_i.copy())

            if log["ham"] == "enthalpy":
                e_iplus1 = self.energy_finder1(x_iplus1, temp)
            elif log["ham"] == "bonds":
                e_iplus1 = self.energy_finder(x_iplus1)
            wandb.log({"Energy_at_Each_point": e_iplus1/initial_config.atoms})
            if e_iplus1 < e_i or self.boltzmann_probability(e_iplus1 - e_i, temp) >= random.random():
                x_i = x_iplus1
                e_i = e_iplus1
                wandb.log({'Structure': x_i, 'Energy': e_iplus1/initial_config.atoms})
                if count % 2 == 0:
                    self.steps.append(i)
                    self.energy_trajectory.append(e_iplus1)
                    self.structure_trajectory.append(x_i)

                count += 1

        return x_i

    @property
    def protocol_single_temp(self) -> np.array:
        system = '-'.join(list(self.initial_config.ele_list))
        log = {"ham": "bonds"}
        wandb.init(
            project="Monte-Carlo",
            group=log["ham"],
            name = str(system) + str(initial_config.atoms),
            tags = [str(self.config_dict["T"]), str(initial_config.ele_dict)]
        )

        print("Warmup Run:")
        x_warm = mc.mc_single_temp(n_trails=self.config_dict["n_warm"],
                                   temp=self.config_dict["warm_T"],
                                   lattice=initial_config.final_bcclattice,
                                   log=log)
        self.steps.append(100)
        self.energy_trajectory.append(100)
        self.structure_trajectory.append(np.zeros_like(x_warm))

        print("Main Run:")
        x_final = mc.mc_single_temp(n_trails=self.config_dict["n_trails"],
                                    temp=self.config_dict["T"],
                                    lattice=x_warm,
                                    log=log)
        self.logger()

        return x_final


lookup = {'Cr-W': -11.129, 'Cr-Ta': -10.562, 'Cr-V': -9.288, 'Ta-W': -12.454, 'V-Ta': -10.332, 'V-W': -11.015, 'Cr-Ti': -8.57, 'Cr-Hf': -9.538, 'Ta-Hf': -10.763, 'Ta-Ti': 0, 'V-Hf': -9.253, 'V-Ti': -8.291, 'W-Hf': -11.313, 'W-Ti': -10.382, 'Ti-Hf': 0, 'Cr-Cr': -9.52, 'W-W': -12.96, 'Ta-Ta': -11.82, 'Hf-Hf': -9.91, 'Ti-Ti': -7.79, 'V-V': -8.97}
rep_unit = 8
k = {
    'V': 0.5,
    'W': 0.5,
}

ele_dict = k
initial_config = InitialConfig(rep_unit, ele_dict, lookup_dict=lookup)
print(f"Creating Initial Configuration of {'-'.join(list(ele_dict.keys()))}")
print(f"{initial_config.atoms} Atoms in total")

mc = MonteCarlo(initial_config, {'n_warm': 50000, 'warm_T': 2500, 'T': 300, 'n_trails': 1000000})
_ = mc.protocol_single_temp
