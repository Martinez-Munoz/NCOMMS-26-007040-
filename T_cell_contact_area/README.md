# Contact Area of Human Cells Workflow

In order to measure the Human Cell Contact Area, we followed the following steps

## 1. Train a Pixel Clasification Machine Learning Model using Labkit Plugin in FIJI
The investigator anotated one image per experiment using the Labkit drawing tool, setting 2 categories: _cell contact area and background_. A _pixel clasification model_ is then generated and saved.

## 2. Segmentation of the dataset 
The ImageJ macro `2. Apply classifier to folder.ijm` performs segmentation of all images in a folder using a given Labkit classifier. It is applied to every folder of our image dataset and save3 the results for the subsequent steps.

## 3. Manual Curation
The investigator checked manually all the segmented images, correcting segmentation defects manually, via drawing new segmented cells or deleting false positive structures. Although the segmentation was overall accurate and this step was not very time consuming.

## 4. Measuring of features from segmented images
The ImageJ macro `4. Extract Features from segmentations` read all curated images and using the plugin `Analyze Particles` measures the features of the contact area of the segmented cells. Results are saved as csv files

## 5. Processing of Data and Statistical Analysis
FIJI outputs the measuring results as a csv file containing the measured features for each cell in a given image. All this files must be grouped into a Data Frame following an appropiate structure in order to perform an ulterior statistical analysis. The Rmarkdown file `5. Processing.Rmd` read all csv files generated in step 4, concatenate the measurements, generated the categorical variables neccesary for the final statistical analysis, and summarised for the reported area measurements in the study. The investigator finally plotted the performed the statistical analysis using Graphpad Prism.

# Publication Materials and Methods

Bright-field images must be segmented to assess the contact area of human cells. The experimenter manually annotated a representative image, serving as the training data for a Machine Learning Pixel Classifier. This classifier was implemented using the Labkit tool [1] within the FIJI software and applied consistently to all images. Segmented images were inspected and any segmentation defects were manually corrected. Quality of segmentation can be checked in suplementary data. Contact sites smaller than 20 microns² were filtered out, as they were presumed to be debris. Finally, segmented binary images were measured and data was exported in CSV format. ImageJ scripts to automatize this workflow were used and are available upon request. We employed Prism software to generate the plots and conduct statistical tests.

[1] Arzt, M., Deschamps, J., Schmied, C., Pietzsch, T., Schmidt, D., Tomancak, P., … Jug, F. (2022). LABKIT: Labeling and Segmentation Toolkit for Big Image Data. _Frontiers in Computer Science_, _4_. [doi:10.3389/fcomp.2022.777728](https://doi.org/10.3389/fcomp.2022.777728)