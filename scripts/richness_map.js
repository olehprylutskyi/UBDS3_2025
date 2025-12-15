/**
 * Oleh Prylutskyi, GPL3
 * Ukrainian Biodiversity Data Science School - 2025
 * 
 * Predict fungal species richness across Ukraine using Random Forest regression in GEE.
 * Training data: "projects/ee-olegpril12/assets/SPUN/samples_w_richness", column "richness_adj"
 * Predictors: Elevation (ALOS), Soil features (OpenLandMap), ESA tree cover assets, 
 * CLIMATE predictors (BIOCLIM), Global Aridity Index, ESA LAND COVER,
 * Tree canopy height, Global Plant Functional Types (PFT).
 * Output: Fungal richness map for Ukraine.
 * Visualization: viridis D palette.
 */

// 1. Define region of interest (Ukraine boundary)
var ukraine = ee.FeatureCollection('FAO/GAUL/2015/level0')
    .filter(ee.Filter.eq('ADM0_NAME', 'Ukraine'));

// 2. Load training points and filter out missing rchnss_dj
var samples = ee.FeatureCollection('projects/ee-olegpril12/assets/SPUN/samples_w_richness')
    .filter(ee.Filter.notNull(['rchnss_dj']));

// 3. Load predictor datasets
// Elevation (ALOS, v4.1)
var elevation = ee.ImageCollection('JAXA/ALOS/AW3D30/V4_1')
    .select('DSM')
    .mosaic();

// Soil features (OpenLandMap)
var soc = ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02').select('b10'); // Soil Organic Carbon (10cm)
var clay = ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02').select('b10'); // Clay Content (10cm)
var ph = ee.Image('OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02').select('b10'); // Soil pH in H2O (10cm)

// ESA tree cover datasets
var deciduous = ee.Image('projects/ee-olegpril12/assets/global_data/ESA_TREES-BD_global');
var needleleaved = ee.Image('projects/ee-olegpril12/assets/global_data/ESA_TREES-NE_global');

// CLIMATE predictors (BIOCLIM)
var bioclim = ee.Image('WORLDCLIM/V1/BIO');
var annualMeanTemperature = bioclim.select('bio01').multiply(0.1).rename('bio01').toFloat();
var temperatureSeasonality = bioclim.select('bio04').multiply(0.01).rename('bio04').toFloat();
var precipitation = bioclim.select('bio12').rename('bio12').toFloat();
var precipitationSeasonality = bioclim.select('bio15').rename('bio15').toFloat();

// Global Aridity Index
var aridity = ee.Image("projects/sat-io/open-datasets/global_ai/global_ai_yearly")
    .multiply(0.0001)
    .rename('aridity')
    .toFloat();

// ESA LAND COVER
var esa = ee.ImageCollection('ESA/WorldCover/v200').first()
    .select('Map').multiply(0.1)
    .rename('landcover')
    .toFloat();

// Tree canopy height
var canopy_ht = ee.ImageCollection("projects/meta-forest-monitoring-okw37/assets/CanopyHeight")
    .mosaic()
    .rename('canopy_ht')
    .toFloat();

// Global Plant Functional Types (PFT)
var trees_bd = ee.Image('projects/ee-olegpril12/assets/global_data/ESA_TREES-BD_global')
    .rename('trees_bd').toFloat();
var trees_ne = ee.Image('projects/ee-olegpril12/assets/global_data/ESA_TREES-NE_global')
    .rename('trees_ne').toFloat();

// 4. Stack predictors into one image
var predictors = elevation.rename('elevation')
    .addBands(soc.rename('soc'))
    .addBands(clay.rename('clay'))
    .addBands(ph.rename('ph'))
    .addBands(deciduous.rename('deciduous'))
    .addBands(needleleaved.rename('needleleaved'))
    .addBands(annualMeanTemperature)
    .addBands(temperatureSeasonality)
    .addBands(precipitation)
    .addBands(precipitationSeasonality)
    .addBands(aridity)
    .addBands(esa)
    .addBands(canopy_ht)
    .addBands(trees_bd)
    .addBands(trees_ne);

// Define scale factor for sampling and export
var scaleFactor = 500; // meters

// 5. Sample predictors at training points
var training = predictors.sampleRegions({
    collection: samples,
    properties: ['rchnss_dj'],
    scale: scaleFactor,
    geometries: true
});

// 6. Train Random Forest regressor
var rf = ee.Classifier.smileRandomForest(100)
    .setOutputMode('REGRESSION')
    .train({
        features: training,
        classProperty: 'rchnss_dj',
        inputProperties: [
            'elevation', 'soc', 'clay', 'ph', 'deciduous', 'needleleaved',
            'bio01', 'bio04', 'bio12', 'bio15', 'aridity',
            'landcover', 'canopy_ht', 'trees_bd', 'trees_ne'
        ]
    });

// 7. Apply model to predictors over Ukraine
var richness_pred = predictors.classify(rf);

// 8. Mask and clip to Ukraine
var richness_map = richness_pred.clip(ukraine);

// 9. Visualization parameters (viridis D palette)
var visParams = {
    min: -20, // Set manually or use a reasonable estimate
    max: 20, // Set manually or use a reasonable estimate
    palette: [
        '#440154', '#482878', '#3E4989', '#31688E',
        '#26828E', '#1F9E89', '#35B779', '#6CCE59',
        '#B4DE2C', '#FDE725'
    ] // viridis D
};

// Set Google Hybrid as default background map
Map.setOptions('HYBRID');

// 10. Add layers to map
Map.centerObject(ukraine, 6);
Map.addLayer(richness_map, visParams, 'Predicted Fungal Richness');

// Add sampled training points to the map
Map.addLayer(training, {color: 'red'}, 'Training Samples');

// Optional: Export result
// Be aware of export limits and adjust region/scale accordingly
// Larger region or finer scale may require more processing time
// or caused GEE's limits to be exceeded.
Export.image.toDrive({
    image: richness_map,
    description: 'Ukraine_Fungal_Richness_RF',
    scale: scaleFactor,
    region: ukraine.geometry(),
    maxPixels: 1e13
});

// 11. Assess prediction quality

// RMSE calculation on training data
// Root Mean Square Error (RMSE) is a common metric to evaluate regression models.
// Lower RMSE values indicate a better fit between predicted and actual data.
// RMSE has the same units as the variable being predicted.
var predicted = training.map(function(feat) {
  var pred = richness_pred.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: feat.geometry(),
    scale: scaleFactor
  }).get('classification');
  return feat.set('predicted', pred);
});

var diff = predicted.map(function(feat) {
  var obs = feat.get('rchnss_dj');
  var pred = feat.get('predicted');
  var error = ee.Number(obs).subtract(ee.Number(pred));
  var sqError = error.pow(2);
  return feat.set('sqError', sqError);
});

var rmse = ee.Number(
  diff.aggregate_sum('sqError')
    .divide(diff.size())
    .sqrt()
);
print('RMSE (training):', rmse);

// Variable importance
print('Variable importance:', rf.explain());

// Sample predictors and prediction at training points
var sampled = predictors.addBands(richness_map.rename('predicted'))
  .reduceRegions({
    collection: samples,
    reducer: ee.Reducer.first(),
    scale: scaleFactor
  });

// Export sampled predictors and predictions to Google Drive
Export.table.toDrive({
  collection: sampled,
  folder: 'GEE_data',
  description: 'Ukraine_Fungal_Richness_Training_Predictions',
  fileFormat: 'CSV'
});