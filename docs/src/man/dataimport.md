# Data importing

We have a number of ways to import data, besides using any of the additional packages in julia to read files.

```@docs
GeophysicalModelGenerator.screenshot_to_GeoData
GeophysicalModelGenerator.screenshot_to_CartData
GeophysicalModelGenerator.screenshot_to_UTMData
GeophysicalModelGenerator.import_topo
GeophysicalModelGenerator.import_GeoTIFF
GeophysicalModelGenerator.tomo_2_GeoData
```

## InSAR ground motion (EGMS)

The [European Ground Motion Service](https://egms.land.copernicus.eu) provides InSAR-derived ground motion for Europe. The files have to be downloaded by hand from the EGMS Explorer (which requires a free EU Login account), after which they can be read, converted to velocities over a period of choice, and put on a topography or a regular grid. See the [EGMS tutorial](tutorial_EGMS.md) for the download instructions and a worked example for Etna.

```@docs
GeophysicalModelGenerator.import_EGMS
GeophysicalModelGenerator.EGMSData
GeophysicalModelGenerator.EGMS_velocity
GeophysicalModelGenerator.EGMS_on_grid
GeophysicalModelGenerator.EGMS_points
GeophysicalModelGenerator.laea_to_lonlat
GeophysicalModelGenerator.lonlat_to_laea
```
