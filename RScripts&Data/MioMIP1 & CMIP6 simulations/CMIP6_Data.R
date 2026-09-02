
library(reticulate)
library(data.table)
library(ggplot2)

use_virtualenv("r-cmip6", required = TRUE)


intake <- import("intake")
col <- intake$open_esm_datastore("https://storage.googleapis.com/cmip6/pangeo-cmip6.json")

py_run_string("r_csv_path = 'catalog_dump.csv'")
py_run_string("r.col.df.to_csv(r_csv_path, index=False)")

catalog_df <- fread("catalog_dump.csv")
setDT(catalog_df)


candidates_vec <- c("CESM2", "CESM2-WACCM", "MPI-ESM1-2-LR", "MPI-ESM1-2-HR",              
  "AWI-CM-1-1-MR", "HadGEM3-GC31-LL", "UKESM1-0-LL", "IPSL-CM6A-LR", "NorESM2-LM", "NorESM2-MM", "MIROC6", "INM-CM4-8", "GFDL-ESM4")

vars_needed <- c("pr", "evspsbl")
exps_needed <- c("historical", "ssp245", "ssp585")


filtered <- catalog_df[
  source_id %in% candidates_vec &
    experiment_id %in% exps_needed &
    variable_id %in% vars_needed &
    table_id == "Amon"]

wide <- dcast(filtered, source_id + member_id ~ experiment_id + variable_id,
              value.var = "zstore", fun.aggregate = length)

wide[, complete := historical_pr == 1 & historical_evspsbl == 1 &
       ssp245_pr == 1     & ssp245_evspsbl == 1 &
       ssp585_pr == 1     & ssp585_evspsbl == 1]

complete_members <- wide[complete == TRUE]
complete_members[, member_num := as.integer(sub("^r([0-9]+).*", "\\1", member_id))]

final_members <- complete_members[
  order(source_id, member_num)
][, .SD[1], by = source_id, .SDcols = "member_id"]

final_members <- rbind(final_members,
                       data.table(source_id = "NorESM2-LM", member_id = "r1i1p1f1"))


lookup <- merge(catalog_df, final_members, by = c("source_id", "member_id"))

final_paths <- lookup[
  experiment_id %in% exps_needed &
    variable_id %in% vars_needed &
    table_id == "Amon",
  .(source_id, member_id, experiment_id, variable_id, table_id, zstore)
]

noresm_daily <- catalog_df[
  source_id == "NorESM2-LM" & experiment_id == "ssp585" &
    variable_id == "pr" & table_id == "day" & member_id == "r1i1p1f1",
  .(source_id, member_id, experiment_id, variable_id, table_id, zstore)
]

final_paths <- rbind(final_paths, noresm_daily)
final_paths[, needs_monthly_aggregation := table_id == "day"]

fwrite(final_paths, "cmip6_final_paths.csv")


py_script <- '
import xarray as xr
import pandas as pd
import numpy as np
import gcsfs
import regionmask
import warnings
warnings.filterwarnings("ignore")

fs = gcsfs.GCSFileSystem(token="anon")

LAT_MIN, LAT_MAX = -35, -15

paths = pd.read_csv("cmip6_final_paths.csv")

def open_zarr(path):
    mapper = fs.get_mapper(path)
    return xr.open_zarr(mapper, consolidated=True)

land_110 = regionmask.defined_regions.natural_earth_v5_0_0.land_110

results = []
failed_combos = []

for idx, row in paths.iterrows():
    source_id = row["source_id"]
    experiment_id = row["experiment_id"]
    variable_id = row["variable_id"]
    zstore = row["zstore"]
    needs_agg = bool(row["needs_monthly_aggregation"])

    print(f"[{idx+1}/{len(paths)}] {source_id} / {experiment_id} / {variable_id}")

    try:
        ds = open_zarr(zstore)
        da = ds[variable_id]

        years_all = da["time"].dt.year
        if experiment_id == "historical":
            da = da.sel(time=(years_all >= 1950) & (years_all <= 2014))
        else:
            da = da.sel(time=(years_all >= 2015) & (years_all <= 2100))

        if needs_agg:
            da = da.resample(time="MS").mean()

        lat_name = "lat" if "lat" in da.dims else "latitude"
        lon_name = "lon" if "lon" in da.dims else "longitude"

        lat_vals = da[lat_name].values
        if lat_vals[0] < lat_vals[-1]:
            da_sub = da.sel({lat_name: slice(LAT_MIN, LAT_MAX)})
        else:
            da_sub = da.sel({lat_name: slice(LAT_MAX, LAT_MIN)})

        mask = land_110.mask(da_sub[lon_name], da_sub[lat_name])
        land_frac = float((~np.isnan(mask.values)).mean())

        if idx < 3:
            print(f"  lon range: {float(da_sub[lon_name].min())} to {float(da_sub[lon_name].max())}")
            print(f"  land fraction in domain: {land_frac:.3f}")

        if land_frac == 0.0 or land_frac == 1.0:
            print(f"  -> WARNING: suspicious land fraction ({land_frac}) - check lon convention!")

        da_masked = da_sub.where(~np.isnan(mask))

        weights = np.cos(np.deg2rad(da_masked[lat_name]))
        weighted_mean = da_masked.weighted(weights).mean(dim=(lat_name, lon_name), skipna=True)

        years = weighted_mean["time"].dt.year.values
        months = weighted_mean["time"].dt.month.values
        vals_mm_day = weighted_mean.values * 86400.0

        df_ts = pd.DataFrame({
            "source_id": source_id,
            "experiment_id": experiment_id,
            "variable_id": variable_id,
            "year": years,
            "month": months,
            "value_mm_day": vals_mm_day
        })
        results.append(df_ts)

    except Exception as e:
        print(f"  -> FAILED: {e}")
        failed_combos.append((source_id, experiment_id, variable_id, str(e)))
        continue

final_df = pd.concat(results, ignore_index=True)
final_df.to_csv("cmip6_pr_evspsbl_timeseries.csv", index=False)
print("Saved:", final_df.shape)

if failed_combos:
    fail_df = pd.DataFrame(failed_combos, columns=["source_id","experiment_id","variable_id","error"])
    fail_df.to_csv("failed_combos.csv", index=False)
    print("Failures logged to failed_combos.csv:", len(failed_combos))
else:
    print("No failures.")
'

writeLines(py_script, "extract_cmip6_timeseries.py")
py_run_file("extract_cmip6_timeseries.py")


ts <- fread("cmip6_pr_evspsbl_timeseries.csv")

ts_wide <- dcast(ts, source_id + experiment_id + year + month ~ variable_id,
                 value.var = "value_mm_day")

ts_wide[, pe := pr - evspsbl]

ts_wide[, monsoon_year := fifelse(month == 12, year + 1, year)]

djf <- ts_wide[month %in% c(12, 1, 2),
               .(djf_pe = mean(pe, na.rm = TRUE)),
               by = .(source_id, experiment_id, monsoon_year)]

jja <- ts_wide[month %in% c(6, 7, 8),
               .(jja_pe = mean(pe, na.rm = TRUE)),
               by = .(source_id, experiment_id, year)]

setnames(jja, "year", "monsoon_year")

monsoon <- merge(djf, jja, by = c("source_id", "experiment_id", "monsoon_year"))
monsoon[, monsoon_index := djf_pe - jja_pe]

baseline <- monsoon[
  experiment_id == "historical" & monsoon_year >= 1995 & monsoon_year <= 2014,
  .(baseline_mean = mean(monsoon_index, na.rm = TRUE)),
  by = source_id
]

monsoon <- merge(monsoon, baseline, by = "source_id")
monsoon[, anomaly := monsoon_index - baseline_mean]

hist_part   <- monsoon[experiment_id == "historical", .(source_id, monsoon_year, anomaly)]
ssp245_part <- monsoon[experiment_id == "ssp245",     .(source_id, monsoon_year, anomaly)]
ssp585_part <- monsoon[experiment_id == "ssp585",     .(source_id, monsoon_year, anomaly)]

ssp245_full <- rbind(hist_part, ssp245_part); ssp245_full[, scenario := "SSP2-4.5"]
ssp585_full <- rbind(hist_part, ssp585_part); ssp585_full[, scenario := "SSP5-8.5"]

combined <- rbind(ssp245_full, ssp585_full)

fwrite(combined, "monsoon_index_per_model.csv")

hist_only   <- combined[monsoon_year <= 2014 & scenario == "SSP2-4.5"]  
hist_only[, scenario := "Historical"]

ssp245_only <- combined[monsoon_year >= 2015 & scenario == "SSP2-4.5"]
ssp585_only <- combined[monsoon_year >= 2015 & scenario == "SSP5-8.5"]

combined_clean <- rbind(hist_only, ssp245_only, ssp585_only)

ens <- combined_clean[, .(
  ens_mean = mean(anomaly, na.rm = TRUE),
  ens_p05  = quantile(anomaly, 0.05, na.rm = TRUE),
  ens_p95  = quantile(anomaly, 0.95, na.rm = TRUE),
  n_models = sum(!is.na(anomaly))
), by = .(scenario, monsoon_year)]

fwrite(ens, "monsoon_index_ensemble_spread.csv")

ens[, scenario := factor(scenario, levels = c("Historical", "SSP2-4.5", "SSP5-8.5"))]

scenario_colors <- c("Historical" = "black", "SSP2-4.5" = "#1f78b4", "SSP5-8.5" = "#e31a1c")

p <- ggplot(ens, aes(x = monsoon_year)) +
  geom_ribbon(aes(ymin = ens_p05, ymax = ens_p95, fill = scenario), alpha = 0.25, color = NA) +
  geom_line(aes(y = ens_mean, color = scenario), linewidth = 1) +
  geom_vline(xintercept = 2015, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  scale_color_manual(values = scenario_colors, name = NULL) +
  scale_fill_manual(values = scenario_colors, name = NULL) +
  labs(
    x = "Year",
    y = "P-E Monsoon Index Anomaly (mm/day)\nrelative to 1995-2014",
    title = "Southern Hemisphere Subtropical Land Monsoon Index (15-35°S)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(4, "pt"),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

print(p)

setwd("RScripts&Data/MioMIP1 & CMIP6 simulations/")
ggsave("Extended Figure 8.png", plot = p, width = 10, height = 6.5, dpi = 600)
