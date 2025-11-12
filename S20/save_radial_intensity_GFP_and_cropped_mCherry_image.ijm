while (nImages>0) { 
  selectImage(nImages); 
  close(); 
} 
setBatchMode(false);
print("\\Clear");

directory = "xxx/"; // location of the image to open

N_objects = 50; // maximum number of ojbects expected in a replicate

position = 0; // the number after 's' for each replicate






condition = newArray("3000", "3000", "3000", "3000", "300", "300", "300", "300", "30", "30", "30", "30", "3", "3", "3", "3");

directory = directory + condition[position-1] + "/s" + position + "/";
mC_image = "a_w3TexRed Camera_s" + toString(position) + "_t60.TIF";
open(directory + mC_image);
selectWindow(mC_image);
rename("mC");
GFP_image = "a_w2GFP Camera_s" + toString(position) + "_t60.TIF";
open(directory + GFP_image);
selectWindow(GFP_image);
rename("GFP");

//setSlice(Math.round(nSlices/2));
for (i=0; i<N_objects; i++){
	object_name = "_s" + toString(position) +"_object_" + toString(i+1);
	
	run("Plots...", "width=1000 height=340 font=14 draw_ticks list minimum=0 maximum=0 interpolate");
	
	selectWindow("GFP");
	resetMinAndMax();
	title = "Select cropping area of object " + toString(i+1);
	msg = "Use the \"Oval\" tool to outline the liposome,\nthen run \"Radial Profile\" and click \"OK\".";
	waitForUser(title, msg);

	data_name = "Values" + object_name + ".csv";
	selectWindow("Plot Values");
	saveAs("Results", directory + data_name);
	close("Radial Profile Plot");
	close(data_name);
	
	selectWindow("mC");
	run("Duplicate...", "title=mC_crop");
	
	selectWindow("mC_crop");
	run("Restore Selection");
	run("Crop");
	image_name = "mC_crop" + object_name + ".tif";
	saveAs("Tiff", directory + image_name);
	close(image_name);
	
	print("Object " + toString(i+1));
}
//print("\n");