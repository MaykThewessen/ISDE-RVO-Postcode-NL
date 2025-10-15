# ISDE Heat Pump Distribution Map for the Netherlands
# Author: Mayk
# Description:
#   Automatically generates a map of the Netherlands with heat pump distribution per postcode
#   using ISDE subsidy data (RVO.nl) and open postcode boundary data from CBS.

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import requests
from io import BytesIO
import os
import folium
from folium import plugins
import branca.colormap as cm
clear = os.system('cls' if os.name == 'nt' else 'clear')

# ---------------------------
# 1. Read ISDE subsidy data
# ---------------------------
# Ensure file path or URL to ISDE Excel file
isde_path = "ISDE downloadbestand augustus 2025.xlsx"
isde = pd.read_csv("ISDE_Warmtepomp_all_rows.csv", delimiter=',', encoding='utf-8')
#isde = pd.DataFrame(isde)

print(isde)
print(isde.describe())


# Example column assumption (adjust as needed):
# Columns: ['Postcode', 'Aantal warmtepompen', 'Provincie', 'Gemeente']
# Ensure 4-digit postcodes
isde['Postcode'] = isde['POSTCODE_AGG'].astype(str).str[:4]

# Get year range for title
year_min = isde['SUBSIDIEJAAR'].min()
year_max = isde['SUBSIDIEJAAR'].max()
year_range = f"{int(year_min)}-{int(year_max)}" if year_min != year_max else f"{int(year_max)}"

# Aggregate by postcode to count heat pumps (ALL YEARS COMBINED)
# Use AANTAL_APP column to count number of heat pump installations
isde_agg = isde.groupby('Postcode').agg({
    'AANTAL_APP': 'sum',  # Sum the number of heat pumps across ALL years
    'PLAATS': 'first',    # Keep first place name
    'GEMEENTENAAM': 'first'  # Keep first municipality name
}).reset_index()

isde_agg.rename(columns={'AANTAL_APP': 'Aantal warmtepompen'}, inplace=True)
print("\nAggregated ISDE data (ALL YEARS COMBINED):")
print(isde_agg.head(20))
print(f"\nTotal postcodes with heat pumps: {len(isde_agg)}")
print(f"Total heat pumps ({year_range}): {isde_agg['Aantal warmtepompen'].sum()}")

# ---------------------------
# 2. Load NL postcode boundaries automatically
# ---------------------------
# Download from PDOK WFS service (official Dutch government geodata)
print("Loading postcode boundaries from PDOK...")

# Local cache file
postcode_cache_file = "nl_postcode4_boundaries.geojson"

# Check if cached file exists
if os.path.exists(postcode_cache_file):
    print(f"Loading cached postcode data from {postcode_cache_file}...")
    nl_postcodes = gpd.read_file(postcode_cache_file)
    print(f"✓ Successfully loaded {len(nl_postcodes)} postcode areas from cache")
else:
    # PDOK WFS service for PC4 postcode areas (2024)
    # Note: WFS has pagination limits, so we'll download in batches
    print("Downloading all postcode areas from PDOK (this may take a moment)...")
    
    try:
        all_postcodes = []
        batch_size = 1000
        start_index = 0
        
        while True:
            wfs_url = "https://service.pdok.nl/cbs/postcode4/2024/wfs/v1_0"
            params = {
                'service': 'WFS',
                'version': '2.0.0',
                'request': 'GetFeature',
                'typeName': 'postcode4',
                'outputFormat': 'json',
                'count': str(batch_size),
                'startIndex': str(start_index)
            }
            
            from urllib.parse import urlencode
            full_url = f"{wfs_url}?{urlencode(params)}"
            batch = gpd.read_file(full_url)
            
            if len(batch) == 0:
                break
                
            all_postcodes.append(batch)
            print(f"  Downloaded {len(batch)} postcodes (total so far: {sum(len(df) for df in all_postcodes)})")
            
            if len(batch) < batch_size:
                break
                
            start_index += batch_size
        
        nl_postcodes = gpd.GeoDataFrame(pd.concat(all_postcodes, ignore_index=True))
        print(f"✓ Successfully loaded {len(nl_postcodes)} postcode areas from PDOK")
        
        # Save to cache file for future use
        print(f"Saving postcode data to {postcode_cache_file} for future use...")
        nl_postcodes.to_file(postcode_cache_file, driver='GeoJSON')
        print(f"✓ Cached postcode data saved successfully")
    except Exception as e:
        print(f"Error loading from PDOK: {e}")
        print("\nTrying alternative sources...")
        
        # Try alternative URLs
        alternatives = [
            "https://geodata.nationaalgeoregister.nl/cbspostcode4/wfs?request=GetFeature&service=WFS&version=2.0.0&typeName=cbs_pc4_2021&outputFormat=json",
            "https://raw.githubusercontent.com/joostdenijs/nl-postcode-boundaries/main/pc4.geojson"
        ]
        
        success = False
        for alt_url in alternatives:
            try:
                nl_postcodes = gpd.read_file(alt_url)
                print(f"Successfully loaded postcode data from alternative source")
                success = True
                break
            except:
                continue
        
        if not success:
            # Try local file
            try:
                nl_postcodes = gpd.read_file("pc4.geojson")
                print("Loaded postcode data from local file")
            except:
                print("\n❌ Could not load postcode data from any source.")
                print("Please download manually from:")
                print("https://www.pdok.nl/datasets/cbs-postcode-4-vlakken")
                raise

# Ensure same key name for merge
# Handle different possible column names from different sources
if 'pc4' in nl_postcodes.columns:
    nl_postcodes['Postcode'] = nl_postcodes['pc4'].astype(str)
elif 'PC4' in nl_postcodes.columns:
    nl_postcodes['Postcode'] = nl_postcodes['PC4'].astype(str)
elif 'postcode' in nl_postcodes.columns:
    nl_postcodes['Postcode'] = nl_postcodes['postcode'].astype(str)
elif 'postcode4' in nl_postcodes.columns:
    nl_postcodes['Postcode'] = nl_postcodes['postcode4'].astype(str)
else:
    print("Available columns:", nl_postcodes.columns.tolist())
    raise ValueError("Could not find postcode column in geodata")

# Ensure both Postcode columns are strings for proper merging
nl_postcodes['Postcode'] = nl_postcodes['Postcode'].astype(str).str.zfill(4)
isde_agg['Postcode'] = isde_agg['Postcode'].astype(str).str.zfill(4)

# ---------------------------
# 3. Load population data per postcode
# ---------------------------
print("\nLoading population data per postcode...")

# Local cache file for population data
population_cache_file = "nl_postcode4_population.csv"

if os.path.exists(population_cache_file):
    print(f"Loading cached population data from {population_cache_file}...")
    population_df = pd.read_csv(population_cache_file, dtype={'Postcode': str})
    population_df['Postcode'] = population_df['Postcode'].astype(str).str.zfill(4)
    print(f"✓ Successfully loaded population data for {len(population_df)} postcodes from cache")
else:
    # Try to extract population data from the existing geojson file
    print("Extracting population data from local postcode boundaries file...")
    
    try:
        if os.path.exists(postcode_cache_file):
            # Extract population data from the geojson file
            all_population = []
            
            for idx, row in nl_postcodes.iterrows():
                # Check for population data in different possible column names
                population = None
                if 'aantalInwoners' in nl_postcodes.columns:
                    population = row.get('aantalInwoners', 0)
                elif 'aantal_inwoners' in nl_postcodes.columns:
                    population = row.get('aantal_inwoners', 0)
                elif 'population' in nl_postcodes.columns:
                    population = row.get('population', 0)
                
                if population is not None and population > 0:
                    all_population.append({
                        'Postcode': row['Postcode'],
                        'Aantal_inwoners': population
                    })
            
            if len(all_population) > 0:
                population_df = pd.DataFrame(all_population)
                population_df['Postcode'] = population_df['Postcode'].astype(str).str.zfill(4)
                population_df = population_df.groupby('Postcode')['Aantal_inwoners'].sum().reset_index()
                
                # Save to cache file
                print(f"✓ Successfully extracted population data for {len(population_df)} postcodes")
                print(f"Saving population data to {population_cache_file} for future use...")
                population_df.to_csv(population_cache_file, index=False)
                print(f"✓ Cached population data saved successfully")
            else:
                print("  ⚠️ No population data found in geojson file")
                population_df = pd.DataFrame({'Postcode': [], 'Aantal_inwoners': []})
        else:
            print("  ⚠️ Postcode boundaries file not found")
            population_df = pd.DataFrame({'Postcode': [], 'Aantal_inwoners': []})
            
    except Exception as e:
        print(f"Error extracting population data: {e}")
        population_df = pd.DataFrame({'Postcode': [], 'Aantal_inwoners': []})

# ---------------------------
# 4. Merge ISDE data with postcode geometry and population
# ---------------------------
merged = nl_postcodes.merge(isde_agg, on="Postcode", how="left")

# Calculate area in km² for each postcode
merged['Oppervlakte_km2'] = merged.geometry.area / 1_000_000  # Convert from m² to km²

# Merge with population data
if len(population_df) > 0:
    merged = merged.merge(population_df, on="Postcode", how="left")
    
    # Fill NaN values with 0 for postcodes without heat pumps
    merged['Aantal warmtepompen'] = merged['Aantal warmtepompen'].fillna(0)
    merged['Aantal_inwoners'] = merged['Aantal_inwoners'].fillna(0)
    
    # Calculate heat pumps per 1000 inhabitants
    merged['Warmtepompen_per_1000_inwoners'] = 0.0
    mask = merged['Aantal_inwoners'] > 0
    merged.loc[mask, 'Warmtepompen_per_1000_inwoners'] = (
        merged.loc[mask, 'Aantal warmtepompen'] / merged.loc[mask, 'Aantal_inwoners'] * 1000
    )
    
    # Calculate heat pumps per km² and round to 1 decimal
    merged['Warmtepompen_per_km2'] = 0.0
    mask_area = merged['Oppervlakte_km2'] > 0
    merged.loc[mask_area, 'Warmtepompen_per_km2'] = (
        merged.loc[mask_area, 'Aantal warmtepompen'] / merged.loc[mask_area, 'Oppervlakte_km2']
    ).round(1)
    
    print(f"\nMerged data: {len(merged)} postcodes")
    print(f"Postcodes with heat pumps: {(merged['Aantal warmtepompen'] > 0).sum()}")
    print(f"Postcodes with population data: {(merged['Aantal_inwoners'] > 0).sum()}")
    print(f"Total population: {merged['Aantal_inwoners'].sum():,.0f}")
    print(f"Average heat pumps per 1000 inhabitants: {merged['Warmtepompen_per_1000_inwoners'].mean():.2f}")
else:
    merged['Aantal warmtepompen'] = merged['Aantal warmtepompen'].fillna(0)
    print(f"\nMerged data: {len(merged)} postcodes")
    print(f"Postcodes with heat pumps: {(merged['Aantal warmtepompen'] > 0).sum()}")

# ---------------------------
# 5. Create matplotlib plots
# ---------------------------
print("\nCreating matplotlib plots...")

# Close all existing figures first
plt.close('all')

# Footnote text for all plots
footnote_text = "Bron: RVO.nl van ISDE downloadbestand augustus 2025.xlsx versie 9-Sept-2025\nFiguur gemaakt door Mayk Thewessen, alle rechten behouden, naamsvermelding bij hergebruik"

# Calculate better color scale based on actual data
color_max = merged['Aantal warmtepompen'].quantile(0.95)

# Plot 1: Absolute numbers - Matplotlib
fig, ax = plt.subplots(figsize=(12, 15))
merged.plot(
    column='Aantal warmtepompen',
    cmap='RdYlGn_r',  # Reversed: groen=weinig, rood=veel
    linewidth=0.1,
    edgecolor='black',
    legend=True,
    legend_kwds={'label': "Aantal warmtepompen per postcode", 'shrink': 0.7},
    ax=ax,
    vmin=0,
    vmax=color_max  # Use same scale as interactive map
)

ax.set_title(f"ISDE Warmtepompen per Postcodegebied (Nederland, {year_range})", fontsize=16, pad=20)
ax.axis("off")

# Add footnote
fig.text(0.5, 0.01, footnote_text, ha='center', fontsize=8, style='italic', wrap=True)

plt.tight_layout(rect=[0, 0.03, 1, 1])  # Make room for footnote

# Save matplotlib plot
matplotlib_pdf_1 = f"ISDE_Warmtepompen_Matplotlib_{year_range}.pdf"
plt.savefig(matplotlib_pdf_1, bbox_inches='tight', facecolor='white')
print(f"✓ Saved matplotlib plot to: {matplotlib_pdf_1}")
plt.close()

# Plot 2: Heat pumps per capita (if population data is available)
if 'Warmtepompen_per_1000_inwoners' in merged.columns:
    fig, ax = plt.subplots(figsize=(12, 15))
    
    # Filter out postcodes with no population data for better visualization
    merged_with_pop = merged[merged['Aantal_inwoners'] > 0].copy()
    
    # Calculate percentile for color scale
    vmax = merged_with_pop['Warmtepompen_per_1000_inwoners'].quantile(0.95)
    
    merged_with_pop.plot(
        column='Warmtepompen_per_1000_inwoners',
        cmap='RdYlGn_r',  # Reversed: groen=weinig, rood=veel
        linewidth=0.1,
        edgecolor='black',
        legend=True,
        legend_kwds={'label': "Warmtepompen per 1000 inwoners", 'shrink': 0.7},
        ax=ax,
        vmin=0,
        vmax=vmax
    )
    
    ax.set_title(f"ISDE Warmtepompen per 1000 Inwoners (Nederland, {year_range})", fontsize=16, pad=20)
    ax.axis("off")
    
    # Add footnote
    fig.text(0.5, 0.01, footnote_text, ha='center', fontsize=8, style='italic', wrap=True)
    
    plt.tight_layout(rect=[0, 0.03, 1, 1])  # Make room for footnote
    
    # Save matplotlib plot
    matplotlib_pdf_2 = f"ISDE_Warmtepompen_PerCapita_Matplotlib_{year_range}.pdf"
    plt.savefig(matplotlib_pdf_2, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved matplotlib plot to: {matplotlib_pdf_2}")
    plt.close()
    
    # Print statistics
    print("\n" + "="*60)
    print("STATISTICS - Heat Pumps per 1000 Inhabitants")
    print("="*60)
    stats_data = merged_with_pop['Warmtepompen_per_1000_inwoners']
    print(f"Mean:   {stats_data.mean():.2f}")
    print(f"Median: {stats_data.median():.2f}")
    print(f"Min:    {stats_data.min():.2f}")
    print(f"Max:    {stats_data.max():.2f}")
    print(f"Std:    {stats_data.std():.2f}")
    print("="*60)
    
    # Plot 3: Heat pumps per km²
    if 'Warmtepompen_per_km2' in merged.columns:
        fig, ax = plt.subplots(figsize=(12, 15))
        
        # Filter out postcodes with no area data for better visualization
        merged_with_area = merged[merged['Oppervlakte_km2'] > 0].copy()
        
        # Set maximum value for color scale
        vmax_km2 = 50  # Cap at 50 warmtepompen per km²
        
        merged_with_area.plot(
            column='Warmtepompen_per_km2',
            cmap='RdYlGn_r',  # Reversed: groen=laag, rood=hoog
            linewidth=0.1,
            edgecolor='black',
            legend=True,
            legend_kwds={'label': "Warmtepompen per km²", 'shrink': 0.7},
            ax=ax,
            vmin=0,
            vmax=vmax_km2
        )
        
        ax.set_title(f"ISDE Warmtepompen per km² (Nederland, {year_range})", fontsize=16, pad=20)
        ax.axis("off")
        
        # Add footnote
        fig.text(0.5, 0.01, footnote_text, ha='center', fontsize=8, style='italic', wrap=True)
        
        plt.tight_layout(rect=[0, 0.03, 1, 1])  # Make room for footnote
        
        # Save matplotlib plot
        matplotlib_pdf_3 = f"ISDE_Warmtepompen_PerKm2_Matplotlib_{year_range}.pdf"
        plt.savefig(matplotlib_pdf_3, bbox_inches='tight', facecolor='white')
        print(f"✓ Saved matplotlib plot to: {matplotlib_pdf_3}")
        plt.close()
        
        # Print statistics
        print("\n" + "="*60)
        print("STATISTICS - Heat Pumps per km²")
        print("="*60)
        stats_data_km2 = merged_with_area['Warmtepompen_per_km2']
        print(f"Mean:   {stats_data_km2.mean():.2f}")
        print(f"Median: {stats_data_km2.median():.2f}")
        print(f"Min:    {stats_data_km2.min():.2f}")
        print(f"Max:    {stats_data_km2.max():.2f}")
        print(f"Std:    {stats_data_km2.std():.2f}")
        print("="*60)
    
    # Create interactive HTML map for per capita visualization
    print("\nCreating interactive HTML map for per capita data...")
    
    # Convert to WGS84 (lat/lon) for folium
    merged_with_pop_wgs84 = merged_with_pop.to_crs(epsg=4326)
    
    # Simplify geometries to reduce file size (tolerance in degrees, ~200m)
    # This significantly reduces the number of coordinate points
    merged_with_pop_wgs84['geometry'] = merged_with_pop_wgs84['geometry'].simplify(tolerance=0.001, preserve_topology=True)
    #merged_with_pop_wgs84['geometry'] = merged_with_pop_wgs84['geometry']
    
    # Round numeric values to reduce file size
    merged_with_pop_wgs84['Warmtepompen_per_1000_inwoners'] = merged_with_pop_wgs84['Warmtepompen_per_1000_inwoners'].round(1)
    merged_with_pop_wgs84['Aantal_inwoners'] = merged_with_pop_wgs84['Aantal_inwoners'].round(0).astype(int)
    
    # Calculate center of the Netherlands
    bounds = merged_with_pop_wgs84.total_bounds
    center_lat = (bounds[1] + bounds[3]) / 2
    center_lon = (bounds[0] + bounds[2]) / 2
    
    # Create folium map
    m = folium.Map(
        location=[center_lat, center_lon],
        zoom_start=7,
        tiles='OpenStreetMap',
        control_scale=True
    )
    
    # Create colormap (groen=weinig, rood=veel)
    vmax_interactive = merged_with_pop_wgs84['Warmtepompen_per_1000_inwoners'].quantile(0.95)
    colormap = cm.LinearColormap(
        colors=['#1a9850', '#91cf60', '#d9ef8b', '#fee08b', '#fc8d59', '#d73027'],
        vmin=0,
        vmax=vmax_interactive,
        caption='Warmtepompen per 1000 inwoners'
    )
    
    # Add choropleth layer
    folium.GeoJson(
        merged_with_pop_wgs84,
        style_function=lambda feature: {
            'fillColor': colormap(feature['properties']['Warmtepompen_per_1000_inwoners']) 
                if feature['properties']['Warmtepompen_per_1000_inwoners'] > 0 else 'lightgray',
            'color': 'black',
            'weight': 0.5,
            'fillOpacity': 0.7
        },
        tooltip=folium.GeoJsonTooltip(
            fields=['Postcode', 'Warmtepompen_per_1000_inwoners', 'Aantal warmtepompen', 'Aantal_inwoners', 'PLAATS', 'GEMEENTENAAM'],
            aliases=['Postcode:', 'Per 1000 inwoners:', 'Totaal warmtepompen:', 'Inwoners:', 'Plaats:', 'Gemeente:'],
            localize=True,
            sticky=False,
            labels=True,
            style="""
                background-color: white;
                border: 2px solid black;
                border-radius: 3px;
                box-shadow: 3px;
            """,
            max_width=300,
        )
    ).add_to(m)
    
    # Add colormap to map
    colormap.add_to(m)
    
    # Add layer control
    folium.LayerControl().add_to(m)
    
    # Add fullscreen button
    plugins.Fullscreen(
        position='topright',
        title='Volledig scherm',
        title_cancel='Sluit volledig scherm',
        force_separate_button=True
    ).add_to(m)
    
    # Save interactive map
    interactive_html = f"ISDE_Warmtepompen_PerCapita_Interactive_{year_range}.html"
    m.save(interactive_html)
    print(f"✓ Saved interactive HTML map to: {interactive_html}")
    
    # Create third map using Choropleth plugin (different visualization style)
    print("\nCreating third interactive HTML map using Choropleth plugin...")
    
    # Create new folium map
    m2 = folium.Map(
        location=[center_lat, center_lon],
        zoom_start=7,
        tiles='CartoDB positron',
        control_scale=True
    )
    
    # Prepare data for choropleth
    # Choropleth requires GeoJSON and a dataframe with matching keys
    choropleth_data = merged_with_pop_wgs84[['Postcode', 'Warmtepompen_per_1000_inwoners']].copy()
    
    # Create choropleth layer using folium.Choropleth
    # Get actual max value to ensure all data is covered
    actual_max = merged_with_pop_wgs84['Warmtepompen_per_1000_inwoners'].max()
    
    folium.Choropleth(
        geo_data=merged_with_pop_wgs84,
        name='Warmtepompen per 1000 inwoners',
        data=choropleth_data,
        columns=['Postcode', 'Warmtepompen_per_1000_inwoners'],
        key_on='feature.properties.Postcode',
        fill_color='RdYlGn_r',  # Reversed: groen=weinig, rood=veel
        fill_opacity=0.7,
        line_opacity=0.5,
        legend_name='Warmtepompen per 1000 inwoners',
        smooth_factor=0,
        threshold_scale=[0, 5, 15, 30, 50, actual_max],
        highlight=True,
        nan_fill_color='lightgray'
    ).add_to(m2)
    
    # Add tooltips using a separate GeoJson layer (Choropleth doesn't support tooltips directly)
    folium.GeoJson(
        merged_with_pop_wgs84,
        style_function=lambda feature: {
            'fillColor': 'transparent',
            'color': 'transparent',
            'weight': 0,
            'fillOpacity': 0
        },
        tooltip=folium.GeoJsonTooltip(
            fields=['Postcode', 'Warmtepompen_per_1000_inwoners', 'Aantal warmtepompen', 'Aantal_inwoners', 'PLAATS', 'GEMEENTENAAM'],
            aliases=['Postcode:', 'Per 1000 inwoners:', 'Totaal warmtepompen:', 'Inwoners:', 'Plaats:', 'Gemeente:'],
            localize=True,
            sticky=False,
            labels=True,
            style="""
                background-color: white;
                border: 2px solid black;
                border-radius: 3px;
                box-shadow: 3px;
            """,
            max_width=300,
        )
    ).add_to(m2)
    
    # Add minimap plugin
    minimap = plugins.MiniMap(toggle_display=True)
    m2.add_child(minimap)
    
    # Add fullscreen button
    plugins.Fullscreen(
        position='topright',
        title='Volledig scherm',
        title_cancel='Sluit volledig scherm',
        force_separate_button=True
    ).add_to(m2)
    
    # Add measure control
    plugins.MeasureControl(position='topleft', primary_length_unit='kilometers').add_to(m2)
    
    # Add layer control
    folium.LayerControl().add_to(m2)
    
    # Save the choropleth map
    interactive_html_choropleth = f"ISDE_Warmtepompen_Choropleth_{year_range}.html"
    m2.save(interactive_html_choropleth)
    print(f"✓ Saved choropleth HTML map to: {interactive_html_choropleth}")
    
    # Create interactive HTML map for heat pumps per km²
    if 'Warmtepompen_per_km2' in merged.columns:
        print("\nCreating interactive HTML map for heat pumps per km²...")
        
        # Filter areas with valid data
        merged_km2_wgs84 = merged[merged['Oppervlakte_km2'] > 0].to_crs(epsg=4326)
        
        # Simplify geometries to reduce file size (tolerance in degrees, ~200m)
        merged_km2_wgs84['geometry'] = merged_km2_wgs84['geometry'].simplify(tolerance=0.001, preserve_topology=True)
        #merged_km2_wgs84['geometry'] = merged_km2_wgs84['geometry']
        
        # Round numeric values to reduce file size
        merged_km2_wgs84['Oppervlakte_km2'] = merged_km2_wgs84['Oppervlakte_km2'].round(1)
        # Warmtepompen_per_km2 is already rounded to 1 decimal
        
        # Create folium map
        m3 = folium.Map(
            location=[center_lat, center_lon],
            zoom_start=7,
            tiles='OpenStreetMap',
            control_scale=True
        )
        
        # Create colormap for km² data (reversed: groen=laag, rood=hoog)
        vmax_km2_interactive = 50  # Cap at 50 warmtepompen per km²
        colormap_km2 = cm.LinearColormap(
            colors=['#1a9850', '#91cf60', '#d9ef8b', '#fee08b', '#fc8d59', '#d73027'],
            vmin=0,
            vmax=vmax_km2_interactive,
            caption='Warmtepompen per km²'
        )
        
        # Add choropleth layer
        folium.GeoJson(
            merged_km2_wgs84,
            style_function=lambda feature: {
                'fillColor': colormap_km2(feature['properties']['Warmtepompen_per_km2']) 
                    if feature['properties']['Warmtepompen_per_km2'] > 0 else 'lightgray',
                'color': 'black',
                'weight': 0.5,
                'fillOpacity': 0.7
            },
            tooltip=folium.GeoJsonTooltip(
                fields=['Postcode', 'Warmtepompen_per_km2', 'Aantal warmtepompen', 'Oppervlakte_km2', 'PLAATS', 'GEMEENTENAAM'],
                aliases=['Postcode:', 'Per km²:', 'Totaal warmtepompen:', 'Oppervlakte (km²):', 'Plaats:', 'Gemeente:'],
                localize=True,
                sticky=False,
                labels=True,
                style="""
                    background-color: white;
                    border: 2px solid black;
                    border-radius: 3px;
                    box-shadow: 3px;
                """,
                max_width=300,
            )
        ).add_to(m3)
        
        # Add colormap to map
        colormap_km2.add_to(m3)
        
        # Add layer control
        folium.LayerControl().add_to(m3)
        
        # Add fullscreen button
        plugins.Fullscreen(
            position='topright',
            title='Volledig scherm',
            title_cancel='Sluit volledig scherm',
            force_separate_button=True
        ).add_to(m3)
        
        # Save interactive map
        interactive_html_km2 = f"ISDE_Warmtepompen_PerKm2_Interactive_{year_range}.html"
        m3.save(interactive_html_km2)
        print(f"✓ Saved interactive HTML map to: {interactive_html_km2}")

print("\n✓ All maps created successfully!")
print(f"  - {matplotlib_pdf_1}")
if 'Warmtepompen_per_1000_inwoners' in merged.columns:
    print(f"  - {matplotlib_pdf_2}")
if 'Warmtepompen_per_km2' in merged.columns:
    print(f"  - {matplotlib_pdf_3}")
    print(f"  - {interactive_html_km2}")
if 'Warmtepompen_per_1000_inwoners' in merged.columns:
    print(f"  - {interactive_html}")
    print(f"  - {interactive_html_choropleth}")

