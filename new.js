var MOD09A1 = ee.ImageCollection("MODIS/006/MOD09A1")
var MOD11A2 = ee.ImageCollection("MODIS/006/MOD11A2")
var MOD13A1 = ee.ImageCollection("MODIS/006/MOD13A1")
var sa = table.geometry()
//去云
function bitwiseExtract(value, fromBit, toBit) {
  if (toBit === undefined) toBit = fromBit
  var maskSize = ee.Number(1).add(toBit).subtract(fromBit)
  var mask = ee.Number(1).leftShift(maskSize).subtract(1)
  return value.rightShift(fromBit).bitwiseAnd(mask)}
  
function cloudfree_mod09a1(image){
  var qa = image.select('StateQA')
  var cloudState = bitwiseExtract(qa, 0, 1) 
  var cloudShadowState = bitwiseExtract(qa, 2)
  var cirrusState = bitwiseExtract(qa, 8, 9)
  var mask = cloudState.eq(0) // Clear
  .and(cloudShadowState.eq(0)) // No cloud shadow
  .and(cirrusState.eq(0)) // No cirrus
  return image.updateMask(mask)}

//定义时间范围
var dr0 = ee.DateRange('2006-06-01','2006-09-30');

var DateRG = ee.List([dr0])

//水体掩膜
function cal_mndwi(image){ return image.normalizedDifference(['sur_refl_b04', 'sur_refl_b06']).rename('mndwi')}
var rencentIMG = MOD09A1.filterDate(dr0).filterBounds(sa).map(cloudfree_mod09a1).mean().clip(sa)
var WaterMask = cal_mndwi(ee.Image(rencentIMG)).lt(0.2)

function GetIMGLST(dr){  //计算LST
  function sts_minmax (image){   //获取输入的最小值和最大值
    var minmax = image.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry:sa,
    scale: 500,
    maxPixels: 1e9}).values();
    return minmax;}
    var rawLST =  MOD11A2.filterDate(dr)
                    .filterBounds(sa).mosaic().clip(sa).select('LST_Day_1km')
                    .multiply(0.02).subtract(273.15).rename('lst').updateMask(WaterMask)
    var minMax = sts_minmax(rawLST)
    var LST = rawLST.unitScale(12,minMax.get(0)) //归一化
    return LST
}
var LSTValue = DateRG.map(GetIMGLST) //计算LST
var LSTbands = ee.ImageCollection(LSTValue).toBands() 

Export.image.toDrive({
  image: LSTbands,
  description: 'LST',
  scale: 500,
  region: sa,
  maxPixels:1e13,
  fileFormat: 'GEOTIFF'
})

function GetIMGNDVI(dr){  //计算NDVI
  function sts_minmax (image){   //获取输入的最小值和最大值
    var minmax = image.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry:sa,
    scale: 500,
    maxPixels: 1e9}).values();
    return minmax;}
    var rawNDVI = MOD13A1.filterDate(dr).filterBounds(sa)
                  .mosaic().clip(sa).select('NDVI')
                  .multiply(0.0001).rename('ndvi').updateMask(WaterMask)
    var minMax = sts_minmax(rawNDVI)  
    var NDVI = rawNDVI.unitScale(minMax.get(1),minMax.get(0)) //归一化
    return NDVI
}
var NDVIValue = DateRG.map(GetIMGNDVI) //计算NDIV
var NDVIbands = ee.ImageCollection(NDVIValue).toBands()

Export.image.toDrive({
  image: NDVIbands,
  description: 'NDVI',
  scale: 500,
  region: sa,
  maxPixels:1e13,
  fileFormat: 'GEOTIFF'
})

function GetIMGWET(dr){  //计算wet
  function sts_minmax (image){   //获取输入的最小值和最大值
    var minmax = image.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry:sa,
    scale: 500,
    maxPixels: 1e9}).values();
    return minmax;}
    //定义研究区去云影像
    var srIMG = MOD09A1.filterDate(dr).filterBounds(sa).map(cloudfree_mod09a1)
                        .mosaic().clip(sa)
    //利用公式计算指标
    var rawWET = srIMG.select(0).multiply(0.1147)
              .add(srIMG.select(1).multiply(0.2489))
              .add(srIMG.select(2).multiply(0.2408))
              .add(srIMG.select(3).multiply(0.3132))
              .add(srIMG.select(4).multiply(-0.3122))
              .add(srIMG.select(5).multiply(-0.6416))
              .add(srIMG.select(6).multiply(-0.5087))
              .multiply(0.0001).rename('wet').updateMask(WaterMask)
    var minMax = sts_minmax(rawWET)
    var WET = rawWET.unitScale(minMax.get(1),minMax.get(0))//归一化
    return WET
}
var WETValue = DateRG.map(GetIMGWET) //计算WET
var WETbands = ee.ImageCollection(WETValue).toBands()

Export.image.toDrive({
  image: WETbands,
  description: 'WET',
  scale: 500,
  region: sa,
  maxPixels:1e13,
  fileFormat: 'GEOTIFF'
})

function GetIMGNDBSI(dr){  //计算NDBSI
  function sts_minmax (image){   //获取输入的最小值和最大值
    var minmax = image.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry:sa,
    scale: 500,
    maxPixels: 1e9}).values();
    return minmax;}
   //定义研究区去云影像
    var srIMG = MOD09A1.filterDate(dr).filterBounds(sa).map(cloudfree_mod09a1)
                        .mosaic().clip(sa) 
    //定义变量
    var swir1 = srIMG.select(5);
    var red = srIMG.select(0); 
    var nir1 = srIMG.select(1);
    var blue = srIMG.select(2);
    var green = srIMG.select(3);
    //利用公式计算指标
    var bi = swir1.add(red).subtract(nir1.add(blue)).divide(swir1.add(red).add(nir1.add(blue)));
    var ibi = swir1.multiply(2).divide(swir1.add(nir1)).subtract(nir1.divide(nir1.add(red)).add(green.divide(green.add(swir1))))
                  .divide(swir1.multiply(2).divide(swir1.add(nir1)).add(nir1.divide(nir1.add(red)).add(green.divide(green.add(swir1)))));
    var rawNDBSI = bi.add(ibi).divide(2).rename('ndbsi').updateMask(WaterMask)
    var minMax = sts_minmax(rawNDBSI)
    var NDBSI = rawNDBSI.unitScale(minMax.get(1),minMax.get(0)) //归一化
    return NDBSI
}
var NDBSIValue = DateRG.map(GetIMGNDBSI) //计算NDBSI
var NDBSIbands = ee.ImageCollection(NDBSIValue).toBands()

Export.image.toDrive({
  image: NDBSIbands,
  description: 'NDBSI',
  scale: 500,
  region: sa,
  maxPixels:1e13,
  fileFormat: 'GEOTIFF'
})

function GetIMG(dr){  //计算四个指标的函数
  var rawIMG = ee.Image()
  function sts_minmax (image){   //获取输入的最小值和最大值
    var minmax = image.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry:sa,
    scale: 500,
    maxPixels: 1e9}).values();
    return minmax;}
  return rawIMG.addBands(GetIMGNDVI(dr))
                .addBands(GetIMGWET(dr))
                .addBands(GetIMGNDBSI(dr))
                .addBands(GetIMGLST(dr)).slice(1,20)
}
var IndexCol = DateRG.map(GetIMG) //计算四个指标
print(IndexCol,'IndexCol')

//pca
function pca_model(image){  
  var scale = 500;
  var bandNames = image.bandNames();
  var region = sa;
  var meanDict = image.reduceRegion({ //均值
    reducer: ee.Reducer.mean(),
    geometry:region,
    scale: scale,
    maxPixels: 1e9});
  var means = ee.Image.constant(meanDict.values(bandNames));
  var centered = image.subtract(means);
  var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  })};
  var arrays = centered.toArray();
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
  });
  var covarArray = ee.Array(covar.get('array'));
  var eigens = covarArray.eigen();
  var eigenValues = eigens.slice(1, 0, 1);
  var eigenVectors = eigens.slice(1, 1);
  var arrayImage = arrays.toArray(1);
  var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);
  var sdImage = ee.Image(eigenValues.sqrt())
    .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);
  return principalComponents
    .arrayProject([0])
    .arrayFlatten([getNewBandNames('pc')])
    .divide(sdImage);
}
var PCA1_result = ee.ImageCollection(IndexCol).map(pca_model).select(0) //主成分分析得PCA1
var PCA2_result = ee.ImageCollection(IndexCol).map(pca_model).select(1) //主成分分析得PCA2
var PCA3_result = ee.ImageCollection(IndexCol).map(pca_model).select(2) //主成分分析得PCA3

print(PCA1_result)

//percentageVariance
function pca_percentageVariance(k){
  var image =  ee.Image(k)
  var scale = 500;
  var bandNames = image.bandNames();
  var region = sa;
  var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry:region,
    scale: scale,
    maxPixels: 1e9
  });
  var means = ee.Image.constant(meanDict.values(bandNames));
  var centered = image.subtract(means);
  var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  })};
  var arrays = centered.toArray();
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
  });
  var covarArray = ee.Array(covar.get('array'));
  var eigens = covarArray.eigen();
  var eigenValues = eigens.slice(1, 0, 1);
   //计算主成分载荷
  var eigenValuesList = eigenValues.toList().flatten()
  var total = eigenValuesList.reduce(ee.Reducer.sum())
  var percentageVariance = eigenValuesList.map(function(item) {
  return (ee.Number(item).divide(total)).multiply(100).format('%.2f')
    })
  return percentageVariance}
var percentageVariance = IndexCol.map(pca_percentageVariance)
print(percentageVariance,'percentageVariance')

//eigenValues
function pca_eigenValues(k){
  var image =  ee.Image(k)
  var scale = 500;
  var bandNames = image.bandNames();
  var region = sa;
  var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry:region,
    scale: scale,
    maxPixels: 1e9
  });
  var means = ee.Image.constant(meanDict.values(bandNames));
  var centered = image.subtract(means);
  var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  })};
  var arrays = centered.toArray();
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
  });
  var covarArray = ee.Array(covar.get('array'));
  var eigens = covarArray.eigen();
  var eigenValues = eigens.slice(1, 0, 1);
  return eigenValues}
  var eigenValues = IndexCol.map(pca_eigenValues)
print(eigenValues,'eigenValues')

//eigenVectors
function pca_eigenVectors(k){
  var image =  ee.Image(k)
  var scale = 500;
  var bandNames = image.bandNames();
  var region = sa;
  var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry:region,
    scale: scale,
    maxPixels: 1e9
  });
  var means = ee.Image.constant(meanDict.values(bandNames));
  var centered = image.subtract(means);
  var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  })};
  var arrays = centered.toArray();
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
  });
  var covarArray = ee.Array(covar.get('array'));
  var eigens = covarArray.eigen();
  var eigenVectors = eigens.slice(1, 1);
  return eigenVectors}
  var eigenVectors = IndexCol.map(pca_eigenVectors)
print(eigenVectors,'eigenVectors')

function normlizedPCA123(image){
  function sts_minmax (image){
    var minmax = image.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry:sa,
    scale: 500,
    maxPixels: 1e9}).values();
    return minmax;}
  var minMax = sts_minmax(image)
  return image.unitScale(minMax.get(1),minMax.get(0)).rename('rsei')}
var RSEI = PCA1_result.map(normlizedPCA123) //归一化得RSEI

var DateRG_start = DateRG.map(function (dr){return ee.DateRange(dr).start()})
var zipRSEI = RSEI.toList(DateRG.length()).zip(DateRG_start)
function setSystime(image){
  var systime = ee.Date(ee.List(image).get(1))
  return ee.Image(ee.List(image).get(0)).set('system:time_start',systime)}

var stimeRSEI = zipRSEI.map(setSystime)
print(stimeRSEI,'stimeRSEI')
var options = {
  width: 400,
  height: 240,
  legend: {position: 'top', textStyle: {color: 'blue', fontSize: 16}},
  lineWidth: 1,
  pointSize: 5, 
  vAxis:{
    title: RSEI
  },
  trendlines: {
    0: {
      type: 'linear',
      color: 'red',
      lineWidth: 1,
      opacity: 0.8,
      showR2: true,
      visibleInLegend: true
    }},
  title: 'Timeserise RSEI',

};

var RSEImean = ui.Chart.image.seriesByRegion({
  imageCollection:stimeRSEI,
  regions:sa,
  reducer: ee.Reducer.mean(),
  scale: 500,
  seriesProperty: "RSEI"
}).setOptions(options)
print(RSEImean,'RSEImean')

var stimeRSEIbands = ee.ImageCollection(stimeRSEI).toBands()
print(stimeRSEIbands,'stimeRSEIbands')

Export.image.toDrive({
  image: stimeRSEIbands,
  description: 'stimeRSEI',
  scale: 500,
  region: sa,
  maxPixels:1e13,
  fileFormat: 'GEOTIFF'
})
