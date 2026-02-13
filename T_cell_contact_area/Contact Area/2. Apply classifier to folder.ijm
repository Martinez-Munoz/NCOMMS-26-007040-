//
// Segments files in a folder using the classifier
//

// 0. Reset FIJI
run("Close All");
run("Clear Results");
print("\\Clear");
counts = roiManager("count");
if(counts !=0) {roiManager("delete");}
run("Options...", "black");
run("Input/Output...", "file=.tsv copy_column copy_row save_column save_row");
run("Set Measurements...", "area centroid center perimeter fit shape feret's area_fraction stack display redirect=None decimal=3");


// 1. Open a dialog to select the pixel classifier and the images folder

Dialog.create("Select:");
Dialog.addDirectory("Folder with raw images", "C:");
Dialog.addFile("Classifier file ", "C:");
Dialog.addDirectory("Saving Folder", "C:");
Dialog.addNumber("Minimun Area Filter in micron^2", 20);
Dialog.show();
dir = Dialog.getString();
clasificador = Dialog.getString();
result_folder = Dialog.getString();
filtro_area = Dialog.getNumber(); 
// Nombre de la carpeta 
folder = split(dir, File.separator);
result_name = folder[folder.length-1]+"_seg_result";


// 2. Segment every image in the folder using the classifier

list = getFileList(dir); // Lista de archivos en la carpeta destino
imagen = 0;
for (i=0; i<list.length; i++){
	if (endsWith(toLowerCase(list[i]), ".tif")){
		imagen = imagen+1;
		// Open and get data
		title=list[i];
		//print(title);
		open(dir+title);
		original = getImageID();
		//run("Properties...", "pixel_width=0.1612 pixel_height=0.1612 voxel_depth=1");
		setVoxelSize(0.1612, 0.1612, 1, "um");
		run("Segment Image With Labkit", "segmenter_file=\'"+clasificador+"\' use_gpu=false");
		setThreshold(0, 0);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Fill Holes"); // rellena los 'holes' si los hubiera
		rename("pre_mask");
		run("Analyze Particles...", "size="+filtro_area+"-Infinity pixel show=Masks exclude clear");
		run("Grays");
		rename("binaria");
		selectImage(original);
		run("8-bit");
		rename("original");
		run("Merge Channels...", "c4=original c6=binaria create ignore");
		rename(title);
		run("Arrange Channels...", "new=21");
		Stack.setChannel(1);
		//run("Set Label...", "label="+title);
		setMetadata("Label", title);
		close("pre_mask");
	}
}

// 3. Concatenate all images

run("Concatenate...", "all_open title="+result_name+" open");
saveAs("Tiff", result_folder+result_name+".tif");
close("*");
print("Terminado");
// getInfo("slice.label")


