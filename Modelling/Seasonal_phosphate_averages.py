import pandas as pd
import os

# -----------------------------
# 1. LOAD DATA
# -----------------------------
file_path = "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Data_PP/Output/Final_With_Dates_Snapped.csv"
df = pd.read_csv(file_path)
df["Date"] = pd.to_datetime(df["Date"])  # ISO format

# -----------------------------
# 2. REMOVE INVALID PHOSPHATE VALUES
# -----------------------------
df["Phosphate (Hanna)"] = df["Phosphate (Hanna)"].replace(0, pd.NA)
df = df.dropna(subset=["Phosphate (Hanna)"])

# -----------------------------
# 3. ASSIGN SEASONS
# -----------------------------
def assign_season(date):
    month = date.month
    year = date.year
    if month in [9, 10, 11]:
        return f"Autumn {year}"
    elif month == 12:
        return f"Winter {year}/{year+1}"
    elif month in [1, 2]:
        return f"Winter {year-1}/{year}"
    elif month in [3, 4, 5]:
        return f"Spring {year}"
    elif month in [6, 7, 8]:
        return f"Summer {year}"

df["Season"] = df["Date"].apply(assign_season)

# -----------------------------
# 4. CALCULATE PER-SITE MEAN & STD
#    Merging by snapped coordinates
# -----------------------------
summary = (
    df.groupby(["Season", "snapped_x", "snapped_y"])
      .agg(
          phosphate_mean=("Phosphate (Hanna)", "mean"),
          phosphate_std=("Phosphate (Hanna)", "std"),
          n_samples=("Phosphate (Hanna)", "count"),
          Site_IDs=("Site ID", lambda x: ", ".join(sorted(x.unique())))
      )
      .reset_index()
)

# Remove invalid site-seasons
summary = summary[
    (summary["phosphate_mean"].notna()) &
    (summary["phosphate_std"].notna()) &
    (summary["phosphate_mean"] > 0) &
    (summary["phosphate_std"] > 0)
]

# -----------------------------
# 5. SEASONAL METRICS ACROSS ALL SITES
# -----------------------------
season_stats = []

for season, season_df in df.groupby("Season"):
    # Remove 0 / NaN values
    season_df = season_df[season_df["Phosphate (Hanna)"].notna() & (season_df["Phosphate (Hanna)"] > 0)]
    
    if season_df.empty:
        continue
    
    # 1. Max, Min, Mean phosphate
    max_phos = season_df["Phosphate (Hanna)"].max()
    min_phos = season_df["Phosphate (Hanna)"].min()
    mean_phos = season_df["Phosphate (Hanna)"].mean()
    
    # 2. Site ID with most samples
    site_counts = season_df.groupby("Site ID")["Phosphate (Hanna)"].count()
    max_site_id = site_counts.idxmax()
    max_count = site_counts.max()
    
    # 3. Site ID with least samples
    min_site_id = site_counts.idxmin()
    min_count = site_counts.min()
    
    season_stats.append({
        "Season": season,
        "Max_Phosphate": max_phos,
        "Min_Phosphate": min_phos,
        "Mean_Phosphate": mean_phos,
        "SiteID_Most_Samples": max_site_id,
        "Count_Most_Samples": max_count,
        "SiteID_Least_Samples": min_site_id,
        "Count_Least_Samples": min_count
    })

season_stats_df = pd.DataFrame(season_stats)

# -----------------------------
# 6. SAVE FILES (WITHOUT SEASON COLUMN)
# -----------------------------
output_dir = "/Users/mahikabhandari/Desktop/Earth Science Year 4/MSci Project/Analysis/Modelling/Phosphate_averages"
os.makedirs(output_dir, exist_ok=True)

target_seasons = [
    "Autumn 2023",
    "Autumn 2024",
    "Spring 2024",
    "Spring 2025",
    "Summer 2024",
    "Summer 2025",
    "Winter 2023/2024",
    "Winter 2024/2025"
]

for season in target_seasons:
    season_df = summary[summary["Season"] == season]
    
    # Select only required columns (drop Season + n_samples)
    season_df = season_df[[
        "Site_IDs",
        "snapped_x",
        "snapped_y",
        "phosphate_mean",
        "phosphate_std"
    ]]
    
    filename = f"{season.replace('/','-')}.csv"
    season_df.to_csv(os.path.join(output_dir, filename), index=False)

print("Seasonal summary saved successfully.")