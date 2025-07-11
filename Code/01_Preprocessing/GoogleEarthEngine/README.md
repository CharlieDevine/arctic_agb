## Google Earth Engine (GEE) multispectral vegetation indices (VI) preprocessing scripts

Links to GEE scripts for generating annual files of different MVIs from MODIS (VARI, NDVI, NIRv, LAI).

Each script performs the following series of steps:
1) Import the image collection containing either the pre-generated MVI or the reflectance bands needed for its computation.
2) Temporally subset MVI image collection to the growing season (June - October) of the year of interest.
3) Spatially subset by clipping to bounding box of extended ABoVE domain using shapefile asset.
4) Compute mean of the spatial/temporal MVI subset.
5) Resample from MODIS sinusoidal projection to WGS84.
6) Export to Google Drive folder.
7) Download locally and resample to exact 0.25-degree domain grid using [MODIS_MVI_Prep.R](https://github.com/CharlieDevine/ABoVE_AGB_VOD/blob/main/Code/Preprocessing/MODIS_MVI_Prep.R)

## GEE links (requires sign-in with Google account)
<a href="https://code.earthengine.google.com/2d509d999fddb8a76535cf48cd2f6d7d?noload=true"> Visible Atmospherically Adjusted Index (VARI) from MCD43A4.006

<a href="https://code.earthengine.google.com/e476da967d6960571d6a2fcc93cd2b75?noload=true"> Normalized Difference Vegetation Index (NDVI) and Near-Infrared Reflectance of Vegetation (NIRv) from MCD43A4.006

<a href="https://code.earthengine.google.com/19d15a211059ae6553f9d1e4bf260862?noload=true"> Leaf Area Index (LAI) from MOD15A2H.006

### VARI
```
// ----------------------------------------------------------------------------------------------------------
// Compute and download VARI data for ABoVE region during summer growing season (Jun-Oct) years 2003 - 2017
// ----------------------------------------------------------------------------------------------------------

var above_geom = ee.FeatureCollection("users/cjdevine/ABoVE/Shapefiles/Extended_Domain/ABoVE_bb");

// Assign domain geometry to variable
var above_bb = above_geom;

// Create function to compute VARI using reflectance bands
var variFun = function(image){
  var vari = image.expression(
    '(g - r) / (g + r - b)', {
      'g' : image.select('Nadir_Reflectance_Band4'), // Green (545-565nm)
      'r' : image.select('Nadir_Reflectance_Band1'), // Red (620-670nm)
      'b' : image.select('Nadir_Reflectance_Band3')  // Blue (459-479nm)
    }).rename('VARI');
  return(vari);
};

// Loop through years 2003 through 2017, compute mean VARI for each year's growing season (June through September)
for (var i = 2003; i < 2017; i++) {
  var dateStart = i+'-06-01';
  var dateEnd = i+'-10-01';
  
  // Filter MCD43A4.006 image collection by each year's date range
  var modis = ee.ImageCollection('MODIS/006/MCD43A4').filterDate(dateStart, dateEnd);
  
  // Compute VARI for each image band using pre-defined function
  var VARI = modis.map(variFun);
  
  // Reduce all image bands using mean to get average VARI for each year's growing season
  var VARImean = VARI.reduce(ee.Reducer.mean());
  print(VARImean, 'Mean VARI '+i);
  
  // Clip to bounding box created for ABoVE region and reproject to WGS84 from MODIS sinusoidal projection
  var VARImean_above = VARImean.clip(above_bb);
  var VARImean_above_wgs = VARImean_above.reproject('EPSG:4326', null, 1000);
  
  // Export to Google Drive folder
  Export.image.toDrive({
  image: VARImean_above_wgs,
  maxPixels: 1.0E13,
  folder: 'ABoVE_MCD43A4_VARI',
  description: 'ABoVE_MCD43A4_VARI_Mean_'+i+'-GS',
  crs: 'EPSG:4326',
  scale: 500,
  region: above_bb,
  fileFormat: 'GeoTIFF'
  });
}
```

### NDVI and NIRv
```
// --------------------------------------------------------------------------------------------------------------------
// Compute and download NDVI and NIRv data for ABoVE region during summer growing season (Jun-Oct) years 2003 - 2017
// --------------------------------------------------------------------------------------------------------------------

var above_geom = ee.FeatureCollection("users/cjdevine/ABoVE/Shapefiles/Extended_Domain/ABoVE_bb");

// Assign domain geometry to variable
var above_bb = above_geom;

// Create function to compute NIRv
var ndviFun = function(image){
  var ndvi = image.expression(
    '(nir - r) / (nir + r)',{
     'nir' : image.select('Nadir_Reflectance_Band2'), // NIR Band (841-876nm)
     'r' : image.select('Nadir_Reflectance_Band1') // Red Band (620-670nm)
    }).rename('NDVI');
  return(ndvi);
};

// Create function to extract NIR
var nirGet = function(image){
  var nir = image.select('Nadir_Reflectance_Band2').rename('NIR'); // NIR Band (841-876nm)
  return(nir);
};

// Loop through years 2003 through 2017, compute mean VARI for each year's growing season (June through September)
for (var i = 2003; i < 2017; i++){
  
  var dateStart = i+'-06-01';
  var dateEnd = i+'-10-01';
  
  // Filter MCD43A4.006 image collection by each year's date range
  var modis = ee.ImageCollection('MODIS/006/MCD43A4').filterDate(dateStart, dateEnd);
  
  // Get NDVI and NIRv for each image band using pre-defined function
  var NIR = modis.map(nirGet);
  var NDVI = modis.map(ndviFun);
  
  // Reduce all image bands using mean to get average NIR and NDVI for each year's growing season
  var NIRmean = NIR.reduce(ee.Reducer.mean());
  var NDVImean = NDVI.reduce(ee.Reducer.mean());
  
  // Clip to bounding box created for ABoVE region and reproject to WGS84 from MODIS sinusoidal projection
  var NIRmean_above = NIRmean.clip(above_bb);
  var NIRmean_above_wgs = NIRmean_above.reproject('EPSG:4326', null, 1000);
  var NDVImean_above = NDVImean.clip(above_bb);
  var NDVImean_above_wgs = NDVImean_above.reproject('EPSG:4326', null, 1000);
  
  // Export NDVI to Google Drive folder
  Export.image.toDrive({
  image: NDVImean_above_wgs,
  maxPixels: 1.0E13,
  folder: 'ABoVE_MCD43A4_NDVI+NIRV',
  description: 'ABoVE_MCD43A4_NDVI_Mean_'+i+'-GS',
  crs: 'EPSG:4326',
  scale: 500,
  region: above_bb,
  fileFormat: 'GeoTIFF'
  });
  
  // Export NIRv to Google Drive folder
  Export.image.toDrive({
  image: NIRmean_above_wgs,
  maxPixels: 1.0E13,
  folder: 'ABoVE_MCD43A4_NDVI+NIRV',
  description: 'ABoVE_MCD43A4_NIR_Mean_'+i+'-GS',
  crs: 'EPSG:4326',
  scale: 500,
  region: above_bb,
  fileFormat: 'GeoTIFF'
  });
}                      
```

### LAI
```
// --------------------------------------------------------------------------------------------------------------------
// Download pre-computed LAI data for ABoVE region during summer growing season (Jun-Oct) years 2003 - 2017
// --------------------------------------------------------------------------------------------------------------------

var above_geom = ee.FeatureCollection("users/cjdevine/ABoVE/Shapefiles/Extended_Domain/ABoVE_bb");

// Assign domain geometry to variable
var above_bb = above_geom;

var laiGet = function(image){
  var lai = image.select('Lai_500m').multiply(0.1).rename('LAI_scaled');
  return(lai);
};

// Loop through years 2003 through 2017, extract mean LAI for each year's growing season (June through September)
for (var i = 2003; i < 2017; i++){
  
  var dateStart = i+'-06-01';
  var dateEnd = i+'-10-01';
  
  // Filter MOD13A1.006 image collection by each year's date range for the growing season
  var modis = ee.ImageCollection('MODIS/006/MOD15A2H').filterDate(dateStart, dateEnd);
  
  // Get LAI for each image band using pre-defined function
  var LAI = modis.map(laiGet);
  
  // Reduce all image bands using mean to get average LAI for each year's growing season
  var LAImean = LAI.reduce(ee.Reducer.mean());
  print(LAImean, 'Mean LAI '+i);
  
  // Clip to bounding box created for ABoVE region and reproject to WGS84 from MODIS sinusoidal projection
  var LAImean_above = LAImean.clip(above_bb);
  var LAImean_above_wgs = LAImean_above.reproject('EPSG:4326', null, 1000);
  
  // Export LAI to Google Drive folder
  Export.image.toDrive({
  image: LAImean_above_wgs,
  maxPixels: 1.0E13,
  folder: 'ABoVE_MOD15A2H_LAI',
  description: 'ABoVE_MOD15A2H_LAI_Mean_'+i+'-GS',
  crs: 'EPSG:4326',
  scale: 500,
  region: above_bb,
  fileFormat: 'GeoTIFF'
  });
}
```
