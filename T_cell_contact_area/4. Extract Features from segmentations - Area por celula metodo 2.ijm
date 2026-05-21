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

setBatchMode(true);
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
		selectWindow("C2-"+title);
		rename("BF_original");
		run("Duplicate...", "title=BF_process duplicate");
		BF = getImageID();
		getDimensions(width, height, channels, slices, frames);
		run("Pseudo flat field correction", "blurring=20 hide stack");
		run("32-bit");
		for (j = 0; j < frames; j++) {
			selectImage(BF);
			Stack.setFrame(j+1);
			run("PHANTAST", "sigma=17 epsilon=0.051 new slice");
			run("Invert");
			run("Fill Holes");
			run("Maximum...", "radius=4 stack");
			run("Watershed");
			rename("temp_stack_"+j);
		
		}
		// generación de los roi y el outline para visualizar
		run("Images to Stack", "name=bin title=temp_stack use");
		Stack.setDimensions(1, 1, frames);
		run("Analyze Particles...", "size=9-Infinity show=Masks exclude add stack");
		rename("outline");
		run("Yellow");
		run("Outline", "stack");
	

		//  medida solo donde el threshold
		selectWindow("C1-"+title);
		rename("contact_bin");
		run("Make Binary", "method=Default background=Dark");
		run("Set Measurements...", "area centroid center perimeter fit shape feret's area_fraction stack limit display redirect=None decimal=2");
		roiManager("measure");
		selectWindow("Results");
		saveAs("Results", result_folder+title+"_metodo2.csv");
	
	
		// pintar los resultados
		run("Merge Channels...", "c4=BF_original c6=contact_bin c7=outline create ignore");
		saveAs("Tiff", result_folder+title+"_metodo2.tif");
		close("*");
		run("Clear Results");
		
	}
}
setBatchMode(false);


print("Terminado");