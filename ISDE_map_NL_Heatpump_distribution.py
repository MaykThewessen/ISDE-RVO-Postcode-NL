"""
ISDE Heat Pump Distribution Map for the Netherlands

Generates choropleth maps of heat pump installations per 4-digit postcode area,
using ISDE subsidy data (RVO.nl) and CBS postcode boundary geodata (PDOK).

Pipeline:
  1. (Optional) Convert raw ISDE Excel → filtered CSV
  2. Load ISDE CSV + postcode boundaries + population data
  3. Aggregate and merge datasets
  4. Produce matplotlib (PNG/PDF) and interactive folium (HTML) maps

Author: Mayk Thewessen
"""

import argparse
import os
from io import BytesIO
from urllib.parse import urlencode

import branca.colormap as cm
import folium
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import requests
from folium import plugins

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ISDE_EXCEL = "ISDE downloadbestand augustus 2025.xlsx"
ISDE_CSV = "ISDE_Warmtepomp_all_rows.csv"
POSTCODE_CACHE = "nl_postcode4_boundaries.geojson"
POPULATION_CACHE = "nl_postcode4_population.csv"

PDOK_WFS_URL = "https://service.pdok.nl/cbs/postcode4/2024/wfs/v1_0"
PDOK_BATCH_SIZE = 1000

FOOTNOTE = (
    "Bron: RVO.nl van ISDE downloadbestand augustus 2025.xlsx versie 9-Sept-2025\n"
    "Figuur gemaakt door Mayk Thewessen, alle rechten behouden, "
    "naamsvermelding bij hergebruik"
)

# Color-scale caps
QUANTILE_CAP = 0.95
KM2_CAP = 50  # heat pumps per km²

FOLIUM_TOOLTIP_STYLE = """
    background-color: white;
    border: 2px solid black;
    border-radius: 3px;
    box-shadow: 3px;
"""


# ---------------------------------------------------------------------------
# 1. Data preparation (Excel → CSV)
# ---------------------------------------------------------------------------
def excel_to_csv(excel_path: str = ISDE_EXCEL, csv_path: str = ISDE_CSV) -> str:
    """Convert the raw ISDE Excel file to a CSV filtered for heat pumps only."""
    print(f"Reading {excel_path} ...")
    df = pd.read_excel(excel_path, engine="openpyxl")
    df_hp = df[df["TECHNIEK"] == "Warmtepomp"]
    df_hp.to_csv(csv_path, index=False)
    print(f"  Wrote {len(df_hp)} heat-pump rows → {csv_path}")
    return csv_path


# ---------------------------------------------------------------------------
# 2. Load datasets
# ---------------------------------------------------------------------------
def load_isde(csv_path: str = ISDE_CSV) -> pd.DataFrame:
    """Load and aggregate ISDE heat-pump data by 4-digit postcode."""
    df = pd.read_csv(csv_path, delimiter=",", encoding="utf-8")
    df["Postcode"] = df["POSTCODE_AGG"].astype(str).str[:4].str.zfill(4)
    print(f"Loaded {len(df)} ISDE rows from {csv_path}")
    return df


def load_postcode_boundaries(cache_path: str = POSTCODE_CACHE) -> gpd.GeoDataFrame:
    """Load PC4 boundaries from cache or download from PDOK WFS."""
    if os.path.exists(cache_path):
        print(f"Loading cached postcode data from {cache_path} ...")
        gdf = gpd.read_file(cache_path)
        print(f"  {len(gdf)} postcode areas loaded")
        return gdf

    print("Downloading postcode areas from PDOK (this may take a moment) ...")
    batches: list[gpd.GeoDataFrame] = []
    start = 0
    while True:
        params = {
            "service": "WFS",
            "version": "2.0.0",
            "request": "GetFeature",
            "typeName": "postcode4",
            "outputFormat": "json",
            "count": str(PDOK_BATCH_SIZE),
            "startIndex": str(start),
        }
        batch = gpd.read_file(f"{PDOK_WFS_URL}?{urlencode(params)}")
        if batch.empty:
            break
        batches.append(batch)
        print(f"  Downloaded {sum(len(b) for b in batches)} postcodes so far")
        if len(batch) < PDOK_BATCH_SIZE:
            break
        start += PDOK_BATCH_SIZE

    gdf = gpd.GeoDataFrame(pd.concat(batches, ignore_index=True))
    gdf.to_file(cache_path, driver="GeoJSON")
    print(f"  Cached {len(gdf)} postcode areas → {cache_path}")
    return gdf


def load_population(
    geo: gpd.GeoDataFrame, cache_path: str = POPULATION_CACHE
) -> pd.DataFrame:
    """Load population per postcode from cache or extract from geodata."""
    if os.path.exists(cache_path):
        df = pd.read_csv(cache_path, dtype={"Postcode": str})
        df["Postcode"] = df["Postcode"].str.zfill(4)
        print(f"Loaded population data for {len(df)} postcodes from {cache_path}")
        return df

    pop_col = next(
        (c for c in ("aantalInwoners", "aantal_inwoners", "population") if c in geo.columns),
        None,
    )
    if pop_col is None:
        print("  No population column found in geodata")
        return pd.DataFrame(columns=["Postcode", "Aantal_inwoners"])

    df = geo[geo[pop_col] > 0][["Postcode", pop_col]].copy()
    df = df.rename(columns={pop_col: "Aantal_inwoners"})
    df["Postcode"] = df["Postcode"].astype(str).str.zfill(4)
    df = df.groupby("Postcode")["Aantal_inwoners"].sum().reset_index()
    df.to_csv(cache_path, index=False)
    print(f"  Extracted & cached population for {len(df)} postcodes → {cache_path}")
    return df


# ---------------------------------------------------------------------------
# 3. Merge & compute derived columns
# ---------------------------------------------------------------------------
def prepare_data(
    isde: pd.DataFrame,
    geo: gpd.GeoDataFrame,
    population: pd.DataFrame,
) -> tuple[gpd.GeoDataFrame, str]:
    """Aggregate ISDE data, merge with geometry + population, compute metrics."""
    year_min = isde["SUBSIDIEJAAR"].min()
    year_max = isde["SUBSIDIEJAAR"].max()
    year_range = (
        f"{int(year_min)}-{int(year_max)}" if year_min != year_max else f"{int(year_max)}"
    )

    agg = (
        isde.groupby("Postcode")
        .agg({"AANTAL_APP": "sum", "PLAATS": "first", "GEMEENTENAAM": "first"})
        .reset_index()
        .rename(columns={"AANTAL_APP": "Aantal warmtepompen"})
    )

    # Normalise postcode column in geodata
    pc_col = next(
        (c for c in ("pc4", "PC4", "postcode", "postcode4") if c in geo.columns), None
    )
    if pc_col is None:
        raise ValueError(f"No postcode column found. Available: {geo.columns.tolist()}")
    geo["Postcode"] = geo[pc_col].astype(str).str.zfill(4)

    merged = geo.merge(agg, on="Postcode", how="left")
    merged["Aantal warmtepompen"] = merged["Aantal warmtepompen"].fillna(0)
    merged["Oppervlakte_km2"] = merged.geometry.area / 1_000_000

    if not population.empty:
        merged = merged.merge(population, on="Postcode", how="left")
        merged["Aantal_inwoners"] = merged["Aantal_inwoners"].fillna(0)
        mask_pop = merged["Aantal_inwoners"] > 0
        merged["Warmtepompen_per_1000_inwoners"] = 0.0
        merged.loc[mask_pop, "Warmtepompen_per_1000_inwoners"] = (
            merged.loc[mask_pop, "Aantal warmtepompen"]
            / merged.loc[mask_pop, "Aantal_inwoners"]
            * 1000
        )
        mask_area = merged["Oppervlakte_km2"] > 0
        merged["Warmtepompen_per_km2"] = 0.0
        merged.loc[mask_area, "Warmtepompen_per_km2"] = (
            merged.loc[mask_area, "Aantal warmtepompen"]
            / merged.loc[mask_area, "Oppervlakte_km2"]
        ).round(1)

    _print_summary(merged, year_range)
    return merged, year_range


def _print_summary(merged: gpd.GeoDataFrame, year_range: str) -> None:
    n_with_hp = (merged["Aantal warmtepompen"] > 0).sum()
    total_hp = int(merged["Aantal warmtepompen"].sum())
    print(f"\nData summary ({year_range}):")
    print(f"  Postcodes total:         {len(merged)}")
    print(f"  Postcodes with heat pumps: {n_with_hp}")
    print(f"  Total heat pumps:        {total_hp:,}")
    if "Aantal_inwoners" in merged.columns:
        print(f"  Total population:        {merged['Aantal_inwoners'].sum():,.0f}")


# ---------------------------------------------------------------------------
# 4a. Matplotlib maps
# ---------------------------------------------------------------------------
def _save_matplotlib(fig: plt.Figure, stem: str, year_range: str) -> None:
    for ext in ("png", "pdf"):
        path = f"{stem}_{year_range}.{ext}"
        fig.savefig(path, bbox_inches="tight", facecolor="white", dpi=150)
        print(f"  Saved {path}")
    plt.close(fig)


def plot_absolute(merged: gpd.GeoDataFrame, year_range: str) -> None:
    vmax = merged["Aantal warmtepompen"].quantile(QUANTILE_CAP)
    fig, ax = plt.subplots(figsize=(12, 15))
    merged.plot(
        column="Aantal warmtepompen",
        cmap="RdYlGn_r",
        linewidth=0.1,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "Aantal warmtepompen per postcode", "shrink": 0.7},
        ax=ax,
        vmin=0,
        vmax=vmax,
    )
    ax.set_title(
        f"ISDE Warmtepompen per Postcodegebied (Nederland, {year_range})",
        fontsize=16,
        pad=20,
    )
    ax.axis("off")
    fig.text(0.5, 0.01, FOOTNOTE, ha="center", fontsize=8, style="italic", wrap=True)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    _save_matplotlib(fig, "ISDE_Warmtepompen_Matplotlib", year_range)


def plot_per_capita(merged: gpd.GeoDataFrame, year_range: str) -> None:
    if "Warmtepompen_per_1000_inwoners" not in merged.columns:
        return
    subset = merged[merged["Aantal_inwoners"] > 0].copy()
    vmax = subset["Warmtepompen_per_1000_inwoners"].quantile(QUANTILE_CAP)

    fig, ax = plt.subplots(figsize=(12, 15))
    subset.plot(
        column="Warmtepompen_per_1000_inwoners",
        cmap="RdYlGn_r",
        linewidth=0.1,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "Warmtepompen per 1000 inwoners", "shrink": 0.7},
        ax=ax,
        vmin=0,
        vmax=vmax,
    )
    ax.set_title(
        f"ISDE Warmtepompen per 1000 Inwoners (Nederland, {year_range})",
        fontsize=16,
        pad=20,
    )
    ax.axis("off")
    fig.text(0.5, 0.01, FOOTNOTE, ha="center", fontsize=8, style="italic", wrap=True)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    _save_matplotlib(fig, "ISDE_Warmtepompen_PerCapita_Matplotlib", year_range)

    _print_stats("Heat Pumps per 1000 Inhabitants", subset["Warmtepompen_per_1000_inwoners"])


def plot_per_km2(merged: gpd.GeoDataFrame, year_range: str) -> None:
    if "Warmtepompen_per_km2" not in merged.columns:
        return
    subset = merged[merged["Oppervlakte_km2"] > 0].copy()

    fig, ax = plt.subplots(figsize=(12, 15))
    subset.plot(
        column="Warmtepompen_per_km2",
        cmap="RdYlGn_r",
        linewidth=0.1,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "Warmtepompen per km²", "shrink": 0.7},
        ax=ax,
        vmin=0,
        vmax=KM2_CAP,
    )
    ax.set_title(
        f"ISDE Warmtepompen per km² (Nederland, {year_range})",
        fontsize=16,
        pad=20,
    )
    ax.axis("off")
    fig.text(0.5, 0.01, FOOTNOTE, ha="center", fontsize=8, style="italic", wrap=True)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    _save_matplotlib(fig, "ISDE_Warmtepompen_PerKm2_Matplotlib", year_range)

    _print_stats("Heat Pumps per km²", subset["Warmtepompen_per_km2"])


def _print_stats(title: str, series: pd.Series) -> None:
    print(f"\n{'=' * 50}")
    print(f"  {title}")
    print(f"{'=' * 50}")
    print(f"  Mean:   {series.mean():.2f}")
    print(f"  Median: {series.median():.2f}")
    print(f"  Min:    {series.min():.2f}")
    print(f"  Max:    {series.max():.2f}")
    print(f"  Std:    {series.std():.2f}")
    print(f"{'=' * 50}")


# ---------------------------------------------------------------------------
# 4b. Interactive folium maps
# ---------------------------------------------------------------------------
def _nl_center(gdf: gpd.GeoDataFrame) -> tuple[float, float]:
    bounds = gdf.total_bounds
    return (bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2


def _make_folium_map(
    gdf_wgs84: gpd.GeoDataFrame,
    value_col: str,
    tooltip_fields: list[str],
    tooltip_aliases: list[str],
    colormap: cm.LinearColormap,
    center: tuple[float, float],
) -> folium.Map:
    """Create a folium Map with a GeoJson choropleth layer."""
    m = folium.Map(
        location=list(center),
        zoom_start=7,
        tiles="OpenStreetMap",
        control_scale=True,
    )
    folium.GeoJson(
        gdf_wgs84,
        style_function=lambda feat, vc=value_col, cmap=colormap: {
            "fillColor": cmap(feat["properties"][vc])
            if feat["properties"][vc] > 0
            else "lightgray",
            "color": "black",
            "weight": 0.5,
            "fillOpacity": 0.7,
        },
        tooltip=folium.GeoJsonTooltip(
            fields=tooltip_fields,
            aliases=tooltip_aliases,
            localize=True,
            sticky=False,
            labels=True,
            style=FOLIUM_TOOLTIP_STYLE,
            max_width=300,
        ),
    ).add_to(m)
    colormap.add_to(m)
    folium.LayerControl().add_to(m)
    plugins.Fullscreen(
        position="topright",
        title="Volledig scherm",
        title_cancel="Sluit volledig scherm",
        force_separate_button=True,
    ).add_to(m)
    return m


def _to_wgs84(gdf: gpd.GeoDataFrame, simplify_tol: float = 0.001) -> gpd.GeoDataFrame:
    gdf_wgs = gdf.to_crs(epsg=4326)
    gdf_wgs["geometry"] = gdf_wgs["geometry"].simplify(
        tolerance=simplify_tol, preserve_topology=True
    )
    return gdf_wgs


def interactive_per_capita(merged: gpd.GeoDataFrame, year_range: str) -> None:
    if "Warmtepompen_per_1000_inwoners" not in merged.columns:
        return
    subset = merged[merged["Aantal_inwoners"] > 0].copy()
    gdf = _to_wgs84(subset)
    gdf["Warmtepompen_per_1000_inwoners"] = gdf["Warmtepompen_per_1000_inwoners"].round(1)
    gdf["Aantal_inwoners"] = gdf["Aantal_inwoners"].round(0).astype(int)

    center = _nl_center(gdf)
    vmax = gdf["Warmtepompen_per_1000_inwoners"].quantile(QUANTILE_CAP)
    cmap = cm.LinearColormap(
        colors=["#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027"],
        vmin=0,
        vmax=vmax,
        caption="Warmtepompen per 1000 inwoners",
    )

    m = _make_folium_map(
        gdf,
        "Warmtepompen_per_1000_inwoners",
        ["Postcode", "Warmtepompen_per_1000_inwoners", "Aantal warmtepompen", "Aantal_inwoners", "PLAATS", "GEMEENTENAAM"],
        ["Postcode:", "Per 1000 inwoners:", "Totaal warmtepompen:", "Inwoners:", "Plaats:", "Gemeente:"],
        cmap,
        center,
    )
    path = f"ISDE_Warmtepompen_PerCapita_Interactive_{year_range}.html"
    m.save(path)
    print(f"  Saved {path}")


def interactive_choropleth(merged: gpd.GeoDataFrame, year_range: str) -> None:
    if "Warmtepompen_per_1000_inwoners" not in merged.columns:
        return
    subset = merged[merged["Aantal_inwoners"] > 0].copy()
    gdf = _to_wgs84(subset)
    gdf["Warmtepompen_per_1000_inwoners"] = gdf["Warmtepompen_per_1000_inwoners"].round(1)
    gdf["Aantal_inwoners"] = gdf["Aantal_inwoners"].round(0).astype(int)

    center = _nl_center(gdf)
    actual_max = gdf["Warmtepompen_per_1000_inwoners"].max()

    m = folium.Map(location=list(center), zoom_start=7, tiles="CartoDB positron", control_scale=True)
    folium.Choropleth(
        geo_data=gdf,
        name="Warmtepompen per 1000 inwoners",
        data=gdf[["Postcode", "Warmtepompen_per_1000_inwoners"]],
        columns=["Postcode", "Warmtepompen_per_1000_inwoners"],
        key_on="feature.properties.Postcode",
        fill_color="RdYlGn_r",
        fill_opacity=0.7,
        line_opacity=0.5,
        legend_name="Warmtepompen per 1000 inwoners",
        smooth_factor=0,
        threshold_scale=[0, 5, 15, 30, 50, actual_max],
        highlight=True,
        nan_fill_color="lightgray",
    ).add_to(m)

    # Transparent tooltip overlay
    folium.GeoJson(
        gdf,
        style_function=lambda _: {"fillColor": "transparent", "color": "transparent", "weight": 0, "fillOpacity": 0},
        tooltip=folium.GeoJsonTooltip(
            fields=["Postcode", "Warmtepompen_per_1000_inwoners", "Aantal warmtepompen", "Aantal_inwoners", "PLAATS", "GEMEENTENAAM"],
            aliases=["Postcode:", "Per 1000 inwoners:", "Totaal warmtepompen:", "Inwoners:", "Plaats:", "Gemeente:"],
            localize=True,
            sticky=False,
            labels=True,
            style=FOLIUM_TOOLTIP_STYLE,
            max_width=300,
        ),
    ).add_to(m)

    plugins.MiniMap(toggle_display=True).add_to(m)
    plugins.Fullscreen(position="topright", title="Volledig scherm", title_cancel="Sluit volledig scherm", force_separate_button=True).add_to(m)
    plugins.MeasureControl(position="topleft", primary_length_unit="kilometers").add_to(m)
    folium.LayerControl().add_to(m)

    path = f"ISDE_Warmtepompen_Choropleth_{year_range}.html"
    m.save(path)
    print(f"  Saved {path}")


def interactive_per_km2(merged: gpd.GeoDataFrame, year_range: str) -> None:
    if "Warmtepompen_per_km2" not in merged.columns:
        return
    subset = merged[merged["Oppervlakte_km2"] > 0].copy()
    gdf = _to_wgs84(subset)
    gdf["Oppervlakte_km2"] = gdf["Oppervlakte_km2"].round(1)

    center = _nl_center(gdf)
    cmap = cm.LinearColormap(
        colors=["#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027"],
        vmin=0,
        vmax=KM2_CAP,
        caption="Warmtepompen per km²",
    )

    m = _make_folium_map(
        gdf,
        "Warmtepompen_per_km2",
        ["Postcode", "Warmtepompen_per_km2", "Aantal warmtepompen", "Oppervlakte_km2", "PLAATS", "GEMEENTENAAM"],
        ["Postcode:", "Per km²:", "Totaal warmtepompen:", "Oppervlakte (km²):", "Plaats:", "Gemeente:"],
        cmap,
        center,
    )
    path = f"ISDE_Warmtepompen_PerKm2_Interactive_{year_range}.html"
    m.save(path)
    print(f"  Saved {path}")


# ---------------------------------------------------------------------------
# 5. Main entry-point
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="ISDE Heat Pump Distribution Maps")
    parser.add_argument(
        "--from-excel",
        action="store_true",
        help="Convert raw ISDE Excel to CSV first (slow, only needed once)",
    )
    args = parser.parse_args()

    # Step 1 — optional Excel conversion
    if args.from_excel:
        excel_to_csv()

    # Step 2 — load data
    isde = load_isde()
    geo = load_postcode_boundaries()
    population = load_population(geo)

    # Step 3 — merge
    merged, year_range = prepare_data(isde, geo, population)

    # Step 4 — matplotlib maps
    print("\nCreating matplotlib plots ...")
    plt.close("all")
    plot_absolute(merged, year_range)
    plot_per_capita(merged, year_range)
    plot_per_km2(merged, year_range)

    # Step 5 — interactive maps
    print("\nCreating interactive HTML maps ...")
    interactive_per_capita(merged, year_range)
    interactive_choropleth(merged, year_range)
    interactive_per_km2(merged, year_range)

    print("\nAll maps created successfully!")


if __name__ == "__main__":
    main()
