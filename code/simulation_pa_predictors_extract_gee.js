/* 
Title: simulation_pa_predictors_extract_gee

By: Patrick Burns [pb463@nau.edu], Northern Arizona University
Collaborators: Zaneta Kaszta, Sam Cushman, Chris Hakkenberg, Patrick Jantz, Andrew Shirk

About: use Google Earth Engine to extact multi-scale predictor variables for SDMs and (optionally) export predictor variable rasters

*/



// ----------------------
// ----- IMPORTS -----
// ----------------------
// Colors
var palettes = require('users/gena/packages:palettes'); 
var color_pal = palettes.colorbrewer.BrBG[11];


// Function to scale predictors variable values between 0 and 1 (currently uses 1st and 99th percentile for min and max)
var scale_preds_0to1 = function(img, scale, aoi){
  var orig_bands = ee.List(img.bandNames()).getInfo() // seems to be necessary to get this to run
  
  //Compute the min and max values per band with (1) the 1st and 99th percentile and (2) absolute min/max
  var perc_red = ee.Reducer.percentile([0,1,99,100])
  var perc_img = img.reduceRegion({geometry: aoi, 
                                  reducer: perc_red, 
                                  scale: scale,
                                  maxPixels: 1e13})

  // Rescale the original image                                
  var img_resc_list = orig_bands.map(function(b){
    var boi = ee.String(b)
    var p0_b = ee.Number(perc_img.get(ee.String(boi.cat('_p0'))))
    var p1_b = ee.Number(perc_img.get(ee.String(boi.cat('_p1'))))
    var p99_b = ee.Number(perc_img.get(ee.String(boi.cat('_p99'))))
    var p100_b = ee.Number(perc_img.get(ee.String(boi.cat('_p100'))))
    
    var pdiff_b = p99_b.subtract(p1_b)
    // Make sure there is a non-zero difference between min/max percentiles. If not, use absolute min/max
    if (pdiff_b === 0){
      var img_b_resc = img.select([b]).unitScale(p0_b,p100_b)
    } else {
      var img_b_resc = img.select([b]).unitScale(p1_b,p99_b)
    }
    return img_b_resc
  })
  var img_resc = ee.Image.cat(img_resc_list)
  return img_resc
}


// Function to standardize predictors by subtracting the mean and dividing by the standard deviation
var standardize_preds_mnsd = function(img, scale, aoi){
  var orig_bands = ee.List(img.bandNames()).getInfo() // seems to be necessary to get this to run
  
  //Compute the mean and standard deviation
  var mnsd_red = img.reduceRegion({geometry: aoi, 
                                   reducer: ee.Reducer.mean()
                                              .combine(ee.Reducer.stdDev(), null, true), 
                                   scale: scale,
                                   maxPixels: 1e13})

  // Standardize the original image                                
  var img_stdz_list = orig_bands.map(function(b){
    var boi = ee.String(b)
    var mn = ee.Number(mnsd_red.get(ee.String(boi.cat('_mean'))))
    var sd = ee.Number(mnsd_red.get(ee.String(boi.cat('_stdDev'))))
    
    var img_b_stdz = img.select([b]).subtract(mn).divide(sd)
    return img_b_stdz
  })
  var img_stdz = ee.Image.cat(img_stdz_list)
  return img_stdz
}


// Function to standardize predictors by dividing by the standard deviation
var standardize_preds_sd = function(img, scale, aoi){
  var orig_bands = ee.List(img.bandNames()).getInfo() // seems to be necessary to get this to run
  
  //Compute the mean and standard deviation
  var sd_red = img.reduceRegion({geometry: aoi, 
                                 reducer: ee.Reducer.stdDev(), 
                                 scale: scale,
                                 maxPixels: 1e13})

  // Standardize the original image                                
  var img_stdz_list = orig_bands.map(function(b){
    var boi = ee.String(b)
    var sd = ee.Number(sd_red.get(b))
    var img_b_stdz = img.select([b]).divide(sd)
    return img_b_stdz
  })
  var img_stdz = ee.Image.cat(img_stdz_list)
  return img_stdz
}

// geom
var box = ee.Geometry.Polygon(
        [[[108.16168533350222, 7.554367650510863],
          [108.16168533350222, -4.493011818555464],
          [119.46661697412722, -4.493011818555464],
          [119.46661697412722, 7.554367650510863]]], null, false);



// ----------------------
// ----- INPUTS -----
// ----------------------

// Which species?
// simulated
var spec = ee.FeatureCollection('projects/ee-gedibio/assets/sdm/borneo_sim_spec_cushman100_kr')
             .filter(ee.Filter.eq('iter_n',1))
             
// Extract values for each point at multiple scales. Note these are radii
var scale_rads = [150,300,600,1200,2400,4800,9600] // 12800 is too big. "Error: Computation timed out. (Error code: 3)"

// Sentinel-2 composite
var s2 = ee.Image('users/pb463/S2_L2A_compos90m_perc_IQR_Borneo_2019-01-01_2020-12-31')
           .select(["Green_p50", "Red_p50", "RE2_p50", "NIR_p50", "SWIR1_p50", "SWIR2_p50", 
                    "kNDVI_p50", "EVI_p50", "RGVI_p50", "NDMI_p50", "NBR_p50","SVVI_p50"])
           .rename(["Imag_Green", "Imag_Red", "Imag_RE2", "Imag_NIR", "Imag_SWIR1", "Imag_SWIR2", 
                    "Imag_kNDVI", "Imag_EVI", "Imag_RGVI", "Imag_NDMI", "Imag_NBR","Imag_SVVI"])
           .set({'min_scale': 90})

// Canopy height (at 2020) from GLAD lab
var can_ht = ee.Image('projects/glad/GLCLU2020/Forest_height_2020').rename('Struct_canht')
               .set({'min_scale': 90})
               
var can_ht_loss = ee.Image('projects/glad/GLCLU2020/Forest_height_netloss').rename('Struct_canhtloss')
                    .set({'min_scale': 90})

// Global Forest Watch Tree Plantations
// ref: https://data.globalforestwatch.org/datasets/gfw::tree-plantations/about
var tree_plant = ee.Image('projects/ee-gedibio/assets/predictors/ag/GFW_Tree_plantations_gen20211210_SEAsia')
                   .rename('Ag_trees')
                   .set({'min_scale': 90})

// Crops
// ref: https://glad.umd.edu/dataset/croplands
var crops_2003 = ee.ImageCollection("users/potapovpeter/Global_cropland_2003")
var crops_2007 = ee.ImageCollection("users/potapovpeter/Global_cropland_2007")
var crops_2011 = ee.ImageCollection("users/potapovpeter/Global_cropland_2011")
var crops_2015 = ee.ImageCollection("users/potapovpeter/Global_cropland_2015")
var crops_2019 = ee.ImageCollection("users/potapovpeter/Global_cropland_2019")
var crops_merge = crops_2003.merge(crops_2007)
                            .merge(crops_2011)
                            .merge(crops_2015)
                            .merge(crops_2019)
                            .reduce(ee.Reducer.anyNonZero())
                            .rename('Ag_crops')
                            .set({'min_scale': 90})

// Human modification
// ref: https://doi.org/10.5194/essd-12-1953-2020
var hm = ee.ImageCollection("CSP/HM/GlobalHumanModification").mosaic().rename('Human_Mod')
           .set({'min_scale': 1000})

// Hunting pressure from Mairin Deith and Jed Brodie (log10 transform applied as recommended by Mairin)
// ref: https://doi.org/10.1098/rspb.2019.2677
var hp = ee.Image('projects/ee-gedibio/assets/predictors/human/Deith_HuntPress_09_07_2021_AllSESources_Summation_nontransformedCurrent')
           .log10()
           .rename('Human_Hunt')
           .set({'min_scale': 1000})

// CHELSA Bioclimatic variables 
var bioclim = ee.Image('projects/ee-gedibio/assets/predictors/climate/CHELSA-v2p1_bioclim_AsiaSE_1981-2010')
var bioclim_new_names = bioclim.bandNames().map(function(n){
  var band_num = ee.String(n).replace('b','')
  return ee.String('bio').cat(ee.String(band_num))
})
bioclim = bioclim.rename(bioclim_new_names)
                 .select(['bio1','bio3','bio12','bio15'])
                 .rename(['Clim_bio1','Clim_bio3','Clim_bio12','Clim_bio15'])
                 .set({'min_scale': 1000})

// Terraclimate
// compute the mean from 1980 to 2010 to match CHELSA
var terraclim = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
                  .filterDate('1980-01-01', '2010-12-31')
                  .mean()
                  .select(['aet', 'def', 'vpd', 'soil', 'srad'])
                  .rename(['Clim_AET', 'Clim_Def', 'Clim_VPD', 'Clim_SoilM', 'Clim_SRad'])
                  .set({'min_scale': 4600})

// Arid_ET
// ref: https://cgiarcsi.community/2019/01/24/global-aridity-index-and-potential-evapotranspiration-climate-database-v2/
var arid_et = ee.Image('projects/ee-gedibio/assets/predictors/climate/ai_et0_v2').rename('Clim_AI')
                .addBands(ee.Image('projects/ee-gedibio/assets/predictors/climate/et0_yr_v2').rename('Clim_ET0'))
                .set({'min_scale': 1000})
                
// Topography
// ref: https://samapriya.github.io/awesome-gee-community-datasets/projects/geomorpho90/
var elev = ee.Image("NASA/NASADEM_HGT/001").select('elevation')
var slope = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/slope").mosaic().rename('Topo_slp');
var twi = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/cti").mosaic().rename('Topo_TWI');
var tri = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/tri").mosaic().rename('Topo_TRI');
var vrm = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/vrm").mosaic().rename('Topo_VRM');
var rough = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/roughness").mosaic().rename('Topo_rgh');
var tpi = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/tpi").mosaic().rename('Topo_TPI');
var hand = ee.ImageCollection("users/gena/global-hand/hand-100").mosaic().rename('Topo_HAND')

var topo = ee.Image.cat([elev, slope, twi, tri, vrm, rough, tpi, hand])
                   .set({'min_scale': 100})

// Soil from SoilGrids, select highest layer in the profile (first band)
var bdod_mean = ee.Image("projects/soilgrids-isric/bdod_mean").select(0).rename('Soil_BDOD');
var cec = ee.Image("projects/soilgrids-isric/cec_mean").select(0).rename('Soil_CEC');
var clay = ee.Image("projects/soilgrids-isric/clay_mean").select(0).rename('Soil_Clay');
var sand = ee.Image("projects/soilgrids-isric/sand_mean").select(0).rename('Soil_Sand');
var nitrogen = ee.Image("projects/soilgrids-isric/nitrogen_mean").select(0).rename('Soil_Nitr');
var phh20 = ee.Image("projects/soilgrids-isric/phh2o_mean").select(0).rename('Soil_pH');
var soc = ee.Image("projects/soilgrids-isric/soc_mean").select(0).rename('Soil_SOC');

var soil = ee.Image.cat([bdod_mean, cec, clay, sand, nitrogen, phh20, soc])
                   .set({'min_scale': 1000}) // nominal resolution is 250 m but I don't buy it

// Country boundaries for clipping and display
var countries = ee.FeatureCollection("USDOS/LSIB/2017")



// ----------------------
// ----- PROCESSING -----
// ----------------------
Map.setOptions('SATELLITE')
Map.addLayer(box,null,'Region of interest')
Map.centerObject(box,8)
// Add layers to map
Map.addLayer(bioclim, {"bands":["Clim_bio1"],"min":12,"max":27}, "Bioclim - MAT", 0)
Map.addLayer(bioclim, {"bands":["Clim_bio12"],"min":1000,"max":6500}, "Bioclim - MAP", 0)
Map.addLayer(topo, {"bands": ['Topo_slp'], "min":0, "max":40, "palette":['ffffff','000000']}, "Topo - Slopeshade", 1)
Map.addLayer(topo, {"bands":["Topo_TWI"],"min":-3.5,"max":3.5}, "Topo - TWI", 0)
Map.addLayer(s2.updateMask(s2.select('Imag_NIR').gt(700)), {"bands":["Imag_SWIR1","Imag_NIR","Imag_Red"],"min":480.1,"max":3656.9,"opacity":0.7}, 'Sentinel-2 Composite', 1)
Map.addLayer(hp, {"min":1.8,"max":4.5,"palette":["273929","f7ff76","ffaf1d","ff0000"]}, "Hunting Pressure")
Map.addLayer(can_ht.updateMask(can_ht.gt(0)), {"min":0, "max": 50}, "GLAD Canopy Height", 0)
Map.addLayer(countries.style({fillColor: 'ffffff00'}), null, "Country Borders")
Map.addLayer(tree_plant, {'min':0, 'max':1}, "Tree Plantations",0)
Map.addLayer(spec.filterMetadata('pa', 'equals', 0), {color: 'red'}, "Species Absence", 1)
Map.addLayer(spec.filterMetadata('pa', 'equals', 1), {color: 'cyan'}, "Species Presence", 1)


// Basic info about presences and absences
// print("----------Species Info----------")
// var pres = spec.filterMetadata('pa', 'equals', 1)
// var abs = spec.filterMetadata('pa', 'equals', 0)
// print("Number of total presences and absences: ", spec.size())
// print("Number of presences: ", pres.size())
// print("Number of absences: ", abs.size())

// add the coordinates as features
// var spec = spec.map(function(f){
//   var x = ee.Number(f.geometry().coordinates().get(0))
//   var y = ee.Number(f.geometry().coordinates().get(1))
//   return f.set({'x': x, 'y': y})
// })


print("----------Model Info----------")

// Combine predictors into a single image and clip to countries
var pred_img_list = [bioclim, terraclim, arid_et, can_ht, can_ht_loss, tree_plant, crops_merge, s2, hm, hp, topo, soil]
var preds_img = ee.Image.cat(pred_img_list)

// Get the names of predictors
var orig_pred_names = ee.List(preds_img.bandNames())
print("Predictor variables: ", orig_pred_names)
print("Scales (radius in meters)", scale_rads)


// Standardize the predictor variables by dividing by the stdDev
var preds_img = standardize_preds_sd(preds_img, 90, box)
                                 
// Calculate focal mean bands
var focal_preds_img_list = scale_rads.map(function(r){
  //var r_str = ee.String(r)
  // Make a list of new names which includes the radius
  var new_names = orig_pred_names.map(function(n){
    var nn = ee.String(n).cat('_mean_').cat(ee.Algorithms.String(r)).cat('m')
    return nn
  })
  
  // Calculate the focal mean of each band using a gaussian kernel
  // The radius corresponds to 2 standard deviations
  var gss_sigma = ee.Number(r).divide(2)
  var kern_gauss = ee.Kernel.gaussian({radius: r, sigma: gss_sigma, units: 'meters', normalize: false})
  var foc = preds_img.reduceNeighborhood({reducer: ee.Reducer.mean(), 
                                          kernel: kern_gauss
                                          }) 
                     .rename(new_names)
  return foc
})
var focal_preds_img = ee.Image.cat(focal_preds_img_list)

// Get the names of the predictors and scales to keep
// We want to exclude bands which have a minimum scales less than the desired set of focal radii
var keep_names = ee.List(pred_img_list.map(function(img){
    // What is the minimum allowable scale
    var min_scale = ee.Number(img.get('min_scale'))
    
    // Get the names of the scales to keep
    var new_names = img.bandNames().map(function(n){
      // Map over potential scales
      var kr = scale_rads.map(function(r) {
        var focal_scale = ee.Number(r).multiply(2)
        return ee.Algorithms.If(min_scale.lt(focal_scale),ee.String(n).cat('_mean_').cat(ee.Algorithms.String(r)).cat('m'),null)
      })
      return ee.List(kr).map(function(i){
        return i
      }, true)
    })
    return new_names.flatten()
})).flatten()
var focal_preds_img = focal_preds_img.select(keep_names)

var pred_names = focal_preds_img.bandNames()
print("Total number of multi-scale predictors", pred_names.size())

// Extract band values from each Pres/Abs point
var pred_ext = focal_preds_img.sampleRegions({collection: spec,
                                              scale: 90,  // resolution to sample raster --25 m value --whatever res to produce resulting map
                                              tileScale: 16,  //memory issues -- set it higher if running out of memory
                                              geometries: true})
                                              
print("Number of features in extract:", pred_ext.size())



// ----------------------
// ----- EXPORT -----
// ----------------------
print("---------- Export ----------")

// Get date+time and add to export file name
var now_date = ee.String(ee.Date(new Date()).format('YYYYMMdd',  'America/Phoenix')).getInfo();
var now_time = ee.String(ee.Date(new Date()).format('HHmmss',  'America/Phoenix')).getInfo();
var exp_prefix = now_date+'_'+now_time+'_Borneo'
print("Switch to the Tasks tab to start model exports.") 
print("The file export prefix is",exp_prefix)

// Export the PA with extracted 
Export.table.toDrive({collection: ee.FeatureCollection(pred_ext), 
                      description: exp_prefix+'_PA_wPreds', 
                      fileNamePrefix: exp_prefix+'_PA_wPreds', 
                      fileFormat: "CSV"
                      })

Export.image.toDrive({image: focal_preds_img, 
                      description: exp_prefix+'_preds_multiscale', 
                      fileNamePrefix: exp_prefix+'_preds_multiscale', 
                      region: box, 
                      scale: 90, 
                      crs: "EPSG:4326", 
                      maxPixels: 1e13, 
                      skipEmptyTiles: true
                      })
