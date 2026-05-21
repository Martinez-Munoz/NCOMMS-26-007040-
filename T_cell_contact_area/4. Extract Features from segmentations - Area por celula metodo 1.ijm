//
// Measuring of features from segmented images
// Images are the results of macro 2. with 2 channels
//

// 0. Reset FIJI
run("Close All");
run("Clear Results");
print("\\Clear");
counts = roiManager("count");
if(counts !=0) {roiManager("delete");}
run("Options...", "black");
run("Input/Output...", "file=.tsv copy_column copy_row save_column save_row");
run("Set Measurements...", "area centroid center perimeter fit shape feret's area_fraction stack display redirect=None decimal=2");


// 1. Open Dialog

Dialog.create("Select:");
Dialog.addDirectory("Folder with segmented images", "C:");
Dialog.addDirectory("Output folder", "C:");
Dialog.show();
dir = Dialog.getString();
result_folder = Dialog.getString();
folder = split(dir, File.separator);
result_name = folder[folder.length-1]+"_seg_result";

// 2. Read segmentation channel of images and apply "analyze particles"

list = getFileList(dir);
imagen = 0;
for (i=0; i<list.length; i++){
	if (endsWith(toLowerCase(list[i]), "seg_result.tif")){
		imagen = imagen+1;
		// Open and get data
		title=list[i];
		//print(title);
		open(dir+title);
		run("Split Channels");
		selectWindow("C1-"+title);
		run("Make Binary", "method=Default background=Dark");
		run("Analyze Particles...", "display clear exclude stack");
		selectWindow("Results");
		saveAs("Results", result_folder+title+".csv");
		close("*");
		run("Clear Results");
	}
}


print("Terminado");