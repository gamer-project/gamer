import os
import shutil

import h5py
import numpy as np


src = "CloudyData_UVB=HM2012.h5"
dst = "CloudyData_NoCooling_mu06.h5"

mmw_dataset = "CoolingRates/Primordial/MMW"
thermal_datasets = [
    "CoolingRates/Metals/Cooling",
    "CoolingRates/Metals/Heating",
    "CoolingRates/Primordial/Cooling",
    "CoolingRates/Primordial/Heating",
]


def print_range(label, data):
    print(f"{label}:")
    print(f"  min = {np.nanmin(data):.15E}")
    print(f"  max = {np.nanmax(data):.15E}")


def set_dataset_value(h5file, name, value):
    if name not in h5file:
        raise RuntimeError(f"Dataset not found: {name}")

    dset = h5file[name]
    old_data = dset[...]

    print(f"\nModify: {name}")
    print_range("Before", old_data)

    dset[...] = value

    new_data = dset[...]
    print_range("After", new_data)

if not os.path.isfile(src):
    raise FileNotFoundError(f"Input file not found: {src}")

shutil.copy2(src, dst)

with h5py.File(dst, "r+") as f:
    print("\n========== MMW / mu ==========")
    set_dataset_value(f, mmw_dataset, 0.6)

    print("\n======= Cooling/Heating =======")
    for dataset_name in thermal_datasets:
        set_dataset_value(f, dataset_name, 0.0)

print(f"\nCreated: {dst}")