from typing import Dict
import numpy as np
import numpy.typing as npt
import netCDF4 as nc
import networkx as nx
import funmixer
import numpy as np
import numpy.typing as npt
import netCDF4 as nc
import networkx as nx
import funmixer
from typing import Dict
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# Build sample network and labels
# --------------------------------------------------

flowdir_file = "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Original/welsh_d8.nc"
sample_csv = "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Modelling/phosphate_averages.csv"

sample_network, labels = funmixer.get_sample_graph(
    flowdirs_filename=flowdir_file,
    sample_data_filename=sample_csv,
)

print(list(sample_network.nodes(data=True))[:1])

def get_unique_upstream_areas(
    sample_network: nx.DiGraph,
    labels: npt.NDArray[np.int_],
) -> Dict[str, npt.NDArray[np.bool_]]:
    return {
        node: labels == node_data["data"].label
        for node, node_data in sample_network.nodes(data=True)
    }

def average_runoff_per_subcatchment(
    sample_network: nx.DiGraph,
    labels: npt.NDArray[np.int_],
    runoff_nc_path: str,
    runoff_var: str = "runoff",
    nodata: float | None = None,
) -> Dict[str, float]:

    with nc.Dataset(runoff_nc_path) as ds:
        runoff = ds.variables[runoff_var][:]

    # Collapse time dimension if present
    if runoff.ndim == 3:
        runoff = np.nanmean(runoff, axis=0)

    # flip runoff on y-axis to match labels orientation
    runoff = np.flipud(runoff)
        
    plt.imshow(runoff)
    plt.colorbar(label='Runoff')
    plt.show()

    plt.imshow(labels)
    plt.colorbar(label='Subcatchment Labels')
    plt.show()

    print("Runoff shape:", runoff.shape)
    print("Labels shape:", labels.shape)

    upstream_masks = get_unique_upstream_areas(sample_network, labels)

    mean_runoff: Dict[str, float] = {}

    for sample_id, mask in upstream_masks.items():
        values = np.array(runoff[mask], copy=True)


        if nodata is not None:
            values = values[values != nodata]

        mean_runoff[sample_id] = (
            float(np.nanmean(values)) if values.size > 0 else np.nan
        )

    #plt.imshow(mean_runoff)
    return mean_runoff

""" Calculate mean runoff per subcatchment for Spring 2024 data """

mean_runoff = average_runoff_per_subcatchment(
    sample_network=sample_network,
    labels=labels,
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2024/runoff_spring2024_mean_wales.nc",
    runoff_var="rowe",  
    nodata=-9999,   
)        

print(mean_runoff)

runoff_df = pd.DataFrame(
    {
        "Site_IDs": list(mean_runoff.keys()),
        "Mean_Runoff": list(mean_runoff.values()),
    }
)

out_csv = (
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2024/runoff_spring2024_mean_by_subcatchment.csv"
)

runoff_df.to_csv(out_csv, index=False)

print(f"✅ Saved Spring 2024 mean runoff per subcatchment to:\n{out_csv}")

""" Calculate mean runoff per subcatchment for Spring 2025 data """

mean_runoff = average_runoff_per_subcatchment(
    sample_network=sample_network,
    labels=labels,
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2025/runoff_spring2025_mean_wales.nc",
    runoff_var="rowe",  
    nodata=-9999,   
)        

print(mean_runoff)

runoff_df = pd.DataFrame(
    {
        "Site_IDs": list(mean_runoff.keys()),
        "Mean_Runoff": list(mean_runoff.values()),
    }
)

out_csv = (
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Spring_2025/runoff_spring2025_mean_by_subcatchment.csv"
)

runoff_df.to_csv(out_csv, index=False)

print(f"✅ Saved Spring 2025mean runoff per subcatchment to:\n{out_csv}")

""" 
Calculate mean runoff per subcatchment for Summer 2024 data 
"""
    
mean_runoff = average_runoff_per_subcatchment(
    sample_network=sample_network,
    labels=labels,
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2024/runoff_summer2024_mean_wales.nc",
    runoff_var="rowe",
    nodata=-9999,   
)      

print(mean_runoff)

runoff_df = pd.DataFrame(
    {
        "Site_IDs": list(mean_runoff.keys()),
        "Mean_Runoff": list(mean_runoff.values()),
    }
)

out_csv = (
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2024/runoff_summer2024_mean_by_subcatchment.csv"
)

runoff_df.to_csv(out_csv, index=False)

print(f"✅ Saved Summer 2024 mean runoff per subcatchment to:\n{out_csv}")

""" 
Calculate mean runoff per subcatchment for Summer 2025 data 
"""
    
mean_runoff = average_runoff_per_subcatchment(
    sample_network=sample_network,
    labels=labels,
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2025/runoff_summer2025_mean_wales.nc",
    runoff_var="rowe",
    nodata=-9999,   
)      

print(mean_runoff)

runoff_df = pd.DataFrame(
    {
        "Site_IDs": list(mean_runoff.keys()),
        "Mean_Runoff": list(mean_runoff.values()),
    }
)

out_csv = (
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Summer_2025/runoff_summer2025_mean_by_subcatchment.csv"
)

runoff_df.to_csv(out_csv, index=False)

print(f"✅ Saved Summer 2025 mean runoff per subcatchment to:\n{out_csv}")

""" 
Calculate mean runoff per subcatchment for Autumn 2023 data 
"""
    
mean_runoff = average_runoff_per_subcatchment(
    sample_network=sample_network,
    labels=labels,
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2023/runoff_autumn2023_mean_wales.nc",
    runoff_var="rowe",
    nodata=-9999,   
)      

print(mean_runoff)

runoff_df = pd.DataFrame(
    {
        "Site_IDs": list(mean_runoff.keys()),
        "Mean_Runoff": list(mean_runoff.values()),
    }
)

out_csv = (
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2023/runoff_autumn2023_mean_by_subcatchment.csv"
)

runoff_df.to_csv(out_csv, index=False)

print(f"✅ Save Autumn 2023 mean runoff per subcatchment to:\n{out_csv}")

""" 
Calculate mean runoff per subcatchment for Autumn 2024 data 
"""
    
mean_runoff = average_runoff_per_subcatchment(
    sample_network=sample_network,
    labels=labels,
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2024/runoff_autumn2024_mean_wales.nc",
    runoff_var="rowe",
    nodata=-9999,   
)      

print(mean_runoff)

runoff_df = pd.DataFrame(
    {
        "Site_IDs": list(mean_runoff.keys()),
        "Mean_Runoff": list(mean_runoff.values()),
    }
)

out_csv = (
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Autumn_2024/runoff_autumn2024_mean_by_subcatchment.csv"
)

runoff_df.to_csv(out_csv, index=False)

print(f"✅ Saved Autumn 2024 mean runoff per subcatchment to:\n{out_csv}")

""" 
Calculate mean runoff per subcatchment for Winter 2023-2024 data 
"""
    
mean_runoff = average_runoff_per_subcatchment(
    sample_network=sample_network,
    labels=labels,
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2023_2024/runoff_winter_2023_2024_mean_wales.nc",
    runoff_var="rowe",
    nodata=-9999,   
)      

print(mean_runoff)

runoff_df = pd.DataFrame(
    {
        "Site_IDs": list(mean_runoff.keys()),
        "Mean_Runoff": list(mean_runoff.values()),
    }
)

out_csv = (
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2023_2024/runoff_winter_2023_2024_mean_by_subcatchment.csv"
)

runoff_df.to_csv(out_csv, index=False)

print(f"✅ Saved Winter 2023-2024 mean runoff per subcatchment to:\n{out_csv}")

""" 
Calculate mean runoff per subcatchment for Winter 2024 - 2025 data 
"""
    
mean_runoff = average_runoff_per_subcatchment(
    sample_network=sample_network,
    labels=labels,
    runoff_nc_path="/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2024_2025/runoff_winter_2024_2025_mean_wales.nc",
    runoff_var="rowe",
    nodata=-9999,   
)      

print(mean_runoff)

runoff_df = pd.DataFrame(
    {
        "Site_IDs": list(mean_runoff.keys()),
        "Mean_Runoff": list(mean_runoff.values()),
    }
)

out_csv = (
    "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Runoff/Winter_2024_2025/runoff_winter_2024_2025_mean_by_subcatchment.csv"
)

runoff_df.to_csv(out_csv, index=False)

print(f"✅ Saved Winter 2024-2025 mean runoff per subcatchment to:\n{out_csv}")