# NCOMMS-26-007040-

# Contact Area of Human Cells

Pipeline and Code in the Contact Area sub-folder.

### Publication Materials and Methods

Bright-field images must be segmented to assess the contact area of human cells. The experimenter manually annotated a representative image, serving as the training data for a Machine Learning Pixel Classifier. This classifier was implemented using the Labkit tool [1] within the FIJI software and applied consistently to all images. Segmented images were inspected and any segmentation defects were manually corrected. Quality of segmentation can be checked in suplementary data. Contact sites smaller than 20 microns² were filtered out, as they were presumed to be debris. Finally, segmented binary images were measured and data was exported in CSV format. ImageJ scripts to automatize this workflow were used and are available upon request. We employed Prism software to generate the plots and conduct statistical tests.

[1] Arzt, M., Deschamps, J., Schmied, C., Pietzsch, T., Schmidt, D., Tomancak, P., … Jug, F. (2022). LABKIT: Labeling and Segmentation Toolkit for Big Image Data. _Frontiers in Computer Science_, _4_. [doi:10.3389/fcomp.2022.777728](https://doi.org/10.3389/fcomp.2022.777728)