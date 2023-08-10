import ij.*;
import ij.plugin.*;
import ij.plugin.LutLoader;
import ij.process.*;
import ij.gui.*;
import ij.measure.*;
import java.util.Arrays;
import java.awt.*;
import java.lang.Math;
import java.awt.image.*;
import java.awt.event.*;
import java.applet.*;
import java.awt.geom.*;
import java.awt.font.*;
import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.FileSystems;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.awt.event.ActionListener;

/**

	The equivalent of Zhou Wang's SSIM MatLab code as a Java plugin inside ImageJ.
	from http://www.cns.nyu.edu/~zwang/files/research/ssim/index.html November 27th  2008.

	Main reference:
	Zhou Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, 
	�Image quality assessment: From error visibility to structural similarity�, 
	IEEE Trans. Image Processing, vol. 13, pp. 600�612, Apr. 2004.

	ImageJ by W. Rasband, U. S. National Institutes of Health, Bethesda, Maryland, USA, 
	http://rsb.info.nih.gov/ij/.  1997-2007. November 27th  2008.

	Java Code by Gabriel Prieto, Margarita Chevalier, Eduardo Guibelalde 26/11/2008.gprietor@med.ucm.es

	Permission to use, copy, or modify this software and its documentation 	for educational and research purposes only and without fee is hereby
	granted, provided that this copyright notice and the original authors' names appear on all copies and supporting documentation. This program
	shall not be used, rewritten, or adapted as the basis of a commercial software or hardware product without first obtaining permission of the
	authors. The authors make no representations about the suitability of this software for any purpose. It is provided "as is" without express
	or implied warranty.

	Please, refer to this version as:

	Gabriel Prieto, Margarita Chevalier, Eduardo Guibelalde. "SSIM Index as a Java plugin for ImageJ"
	Department of Radiology, Faculty of Medicine. Universidad Complutense. Madrid. SPAIN.
	http://www.ucm.es/info/fismed/SSIM_family/SSIM.htm


THIS IS AN AN IMPLEMENTATION AS A PLUGIN FOR IMAGEJ DEVELOPED IN JAVA. IT CALCULATES THE SSIM INDEX (SINGLE-SCALE) AND IT USES THE ZHOU WANG'S ALGORITHM DEVELOPED
IN MATLAB. IT WORKS WITH 8, 16 AND 32 BITS GRAY LEVELS. 
THE DIFERENCE BETWEEN THIS VERSION ANS V_1 IS THE SIMULATION OF VIEWING DISTANCE. ACCORDING TO ZHOU WANG'S HOME PAGE, literal citation:

" The precisely �right� scale depends on both the image resolution and the viewing distance and is usually difficult to be obtained.  In practice, we suggest to use the following empirical formula 
to determine the scale for images viewed from a typical distance (say 3~5 times of the image height): 1) Let F = max(1, round(N/256)), where N is the number of pixels in image height; 
2) Average local F by F pixels and then downsample the image by a factor of F; and 3) apply the ssim_index.m program. For example, for an 512 by 512 image, F = max(1, round(512/256)) = 2, 
so the image should be average within a 2 by 2 window and downsampled by a factor of 2 before applying ssim_index.m "

THIS VERSION INTRODUCES A CONTROL TO GET THE ACTUAL VIEWING DISTANCE. YOU CAN CHANGE THE PROPOSED VALUE. IT IS ONLY A DOWNSAMPLING FACTOR.

*/

public class DSSIM_index implements PlugIn {

	//protected ImagePlus image_1_imp, image_2_imp;  //THIS PLUGIN WORKS WITH TWO (AND ONLY TWO) IMAGES OPEN IN IMAGEJ
	//protected ImageProcessor image_1_p, image_2_p;
	//put the ImagePlus of the stack here or something
	ImagePlus image_imp;

	//ui variables
	double sigma_gauss = 1.5;
	int filter_width = 11;//((int)(2*Math.ceil(sigma_gauss*3)))+1;
	int filter_scale = 10;
	double K1 = 0.01; 
	double K2 = 0.03;
	double downsampled = 1; //downsampling dont know if needed
	double downsampled_backup = downsampled;
	boolean gaussian_window = true;
	String[] window_type = {"Gaussian","Same weight"};  // WE CAN WEIGHTS THE WINDOW WITH A GAUSSIAN WEIGHTING FUNCTION OR GIVING THE SAME WEIGHT TO ALL THE PIXELS IN THE WINDOW
	String window_selection = window_type[0];
	boolean out=false;
	int frame_offset = 1;
	double alpha = 1.0;
	double beta = 1.0;
	double gamma = 1.0;
	boolean show_downsampled_images= false;
	final File f = new File(DSSIM_index.class.getProtectionDomain().getCodeSource().getLocation().getPath());
	String directory= f.toString();
	//String lut_path="C:\\Users\\kpiye\\Desktop\\ImageJ\\luts\\Viridis.lut";
	boolean save_values = false;
	boolean show_ssim_map= true;
	boolean use_color_map= true;
	boolean show_gaussian_filter= false;
	boolean per_frame = true;
	GenericDialog gd = new GenericDialog ("DSSIM");
	
	/**public class FilterSettings implements ActionListener{

		public void actionPerformed(ActionEvent e){
			GenericDialog fst= new GenericDialog ("Filter Settings");

			AdvancedButton advanced = new AdvancedButton();
			AutoFilter autof = new AutoFilter();
			EditSD esd = new EditSD();
		
			//gd.addToSameRow();
			fst.addButton("Standard Deviation", esd);
			fst.addButton("Edit width", advanced);
			fst.addToSameRow();
			fst.addButton("Auto width", autof);

			fst.showDialog();
			if (fst.wasCanceled()) return;
			
			fst.dispose();

		}

	}**/

	public class EditSD implements ActionListener{

		public void actionPerformed(ActionEvent e){
			GenericDialog sdgd= new GenericDialog ("Standard Deviation");
			sdgd.addNumericField ("Standard deviation:", sigma_gauss, 1);
			sdgd.showDialog();
			if (sdgd.wasCanceled()) return;
			sigma_gauss = sdgd.getNextNumber();
			sdgd.dispose();

		}

	}

	public class AutoFilter implements ActionListener{

		public void actionPerformed(ActionEvent e){
			filter_width=((int)(2*Math.ceil(sigma_gauss*3)))+1;
			IJ.error(Integer.toString(filter_width));
			
			
		}

	}

	//Advanced button action
	public class AdvancedButton implements ActionListener{
		boolean first_time = true;
		
		public void actionPerformed(ActionEvent e){
		boolean outt = false;
			while (!outt){
				GenericDialog gdi= new GenericDialog ("Edit width");
				outt=true;
				//GenericDialog gdi= new GenericDialog ("Advanced");
				
				if (first_time){
					first_time = true;
					//gdi.addNumericField ("Filter viewscale:", filter_scale, 0);
					
					gdi.addNumericField ("Filter width:",  filter_width, 0);
					
					
				}
				
				gdi.showDialog();
				if (gdi.wasCanceled()) return;
				//filter_scale = (int) (gdi.getNextNumber());
				
				filter_width = (int) (gdi.getNextNumber());
				
				show_ssim_map = gdi.getNextBoolean();
				per_frame = gdi.getNextBoolean();
				gdi.dispose();
			}
		}
	}

	public class WindowButton implements ActionListener{
		public void actionPerformed(ActionEvent e){  
			
			int fl = filter_width*filter_width;
			float window_weights [] = new float [fl];
			double [] array_gauss_window = new double [fl];

			if (gaussian_window) {
				double value, distance = 0;
				int center = (filter_width/2);
				double total = 0;
				int pointer1;
				double sigma_sq=sigma_gauss*sigma_gauss;
				
					for (int y = 0; y < filter_width; y++){
					for (int x = 0; x < filter_width; x++){
								distance = Math.abs(x-center)*Math.abs(x-center)+Math.abs(y-center)*Math.abs(y-center);
						pointer1 = y*filter_width + x;
									array_gauss_window[pointer1] = Math.exp(-0.5*distance/sigma_sq);
						total = total + array_gauss_window[pointer1];
					}
					}
				for (pointer1=0; pointer1 < fl; pointer1++) {	
					array_gauss_window[pointer1] = array_gauss_window[pointer1] / total;
					window_weights [pointer1] = (float) array_gauss_window[pointer1];
				}
			}
			ColorModel cm=null;
			ImageProcessor gauss_window_ip = new FloatProcessor (filter_width, filter_width, window_weights, cm);
			gauss_window_ip = gauss_window_ip.resize (filter_width*filter_scale);
			String title_filtro_1 = "Sigma: " + sigma_gauss + " Width: "+ filter_width + " pixels"; 	
			ImagePlus gauss_window_imp = new ImagePlus (title_filtro_1, gauss_window_ip);
			//gauss_window_imp.show();
			//gauss_window_imp.updateAndDraw();
			boolean outt = false;
			while (!outt){
				outt = true;
				GenericDialog gdii= new GenericDialog("Gaussian Filter");
				//gdii.setOKLabel("refresh");
				gdii.addImage(gauss_window_imp);
				//gdii.addNumericField ("Filter viewscale:", filter_scale, 0);
				gdii.showDialog();
				if (gdii.wasCanceled()) return;
				//filter_scale = (int) (gdii.getNextNumber());
				//gdii.dispose();
			}
		}
	}
	// public class WindowButton implements ActionListener{
	// 	public void actionPerformed(ActionEvent e){  
	// 		ColorModel cm=null;
	// 		ImageProcessor gauss_window_ip = new FloatProcessor (filter_width, filter_width, window_weights, cm);
	// 		gauss_window_ip = gauss_window_ip.resize (filter_width*filter_scale);
	// 		String title_filtro_1 = "Sigma: " + sigma_gauss + " Width: "+ filter_width + " pixeles"; 	
	// 		ImagePlus gauss_window_imp = new ImagePlus (title_filtro_1, gauss_window_ip);
	// 		gauss_window_imp.show();
	// 		gauss_window_imp.updateAndDraw();
	// 	}
	// }


public void run (String arg) {

	String title_1, title_2;//for downsampling
	int  pointer, filter_length, image_height, image_width, image_dimension, bits_per_pixel_1, bits_per_pixel_2, a, b, c;
	float filter_weights [];
	float af, bf;
	double [] ssim_map;
	double ssim_index;
//
// ERROR CONTROLS. TWO IMAGES SHOULD BE OPENED AND BOTH WITH THE SAME DIMENSIONS
//
	int[] wList = WindowManager.getIDList();
	if (wList==null) {
		IJ.error("There is no image open");
		return;
	}
	image_imp = WindowManager.getCurrentImage();
	//a = WindowManager.getImageCount();

	int stackSize = image_imp.getImageStackSize(); //gonna just check if the stack size is >1
	if (stackSize < 2) {
		IJ.error(String.valueOf(stackSize)+" Stack Size");
		return;
	}

	bits_per_pixel_1=image_imp.getBitDepth();
	/**
	bits_per_pixel_1=image_1_imp.getBitDepth();
	bits_per_pixel_2=image_2_imp.getBitDepth();
	if (bits_per_pixel_1 != bits_per_pixel_2){
		IJ.error("Both images must have the same number of bits per pixel");
		return;
	}**/
	if (bits_per_pixel_1 == 24){
		IJ.error("RGB images are not supportedl");
		return;
	}
//
// END OF CONTROL ERRORS
//	
//
// THIS DIALOG BOX SHOWS DIFFERENT OPTIONS TO CREATE THE WINDOW WE ARE GOING TO USE TO EVALUATE SSIM INDEX OVER THE ENTIRE IMAGES
//	
	image_height = image_imp.getHeight();
	image_width = image_imp.getWidth();
	


	while (!out){
		out=true;
	
		WindowButton windowb  = new WindowButton();

		//FilterSettings filst = new FilterSettings();
		
		AdvancedButton advanced = new AdvancedButton();
		AutoFilter autof = new AutoFilter();
		EditSD esd = new EditSD();
		//gd.addToSameRow();
		gd.addButton("Standard Deviation", esd);
		gd.addButton("Edit width", advanced);
		//gd.addToSameRow();
		gd.addButton("Auto width", autof);

		gd.addButton("Show filter", windowb);
		
		gd.addNumericField ("Frame offset:", frame_offset, 1);
		gd.addNumericField ("Alpha:", alpha, 1);
		gd.addToSameRow();
		gd.addNumericField ("Beta:", beta, 1);
		gd.addToSameRow();
		gd.addNumericField ("Gamma:", gamma, 1);
		gd.addCheckbox("Per-frame contrast", per_frame);
		gd.addCheckbox("Save DSSIM values", save_values);

		gd.addDirectoryField("", directory);
		
	
		//gd.addToSameRow();
		
		
		//gd.addHelp("www.google.com");
		gd.showDialog();
		

		if (gd.wasCanceled()) return;
		
		frame_offset = (int) (gd.getNextNumber());
		alpha = (double) (gd.getNextNumber());
		beta = (double) (gd.getNextNumber());
		gamma = (double) (gd.getNextNumber());
		per_frame = gd.getNextBoolean();
		save_values = gd.getNextBoolean();
		directory = gd.getNextString();
		double d;
		a = filter_width/2;
		d = filter_width -a*2;
		if (window_selection != "Gaussian") gaussian_window = false;
		if (d==0) {
			IJ.error("Filter width and heigth must be odd"); 
			gd = new GenericDialog ("DSSIM");
			out = false;
		}
		if (gaussian_window & sigma_gauss <= 0) {
			IJ.error("Standard deviation must be greater than 0");
			gd = new GenericDialog ("DSSIM");
			out = false;
		}
		if (gaussian_window & filter_scale <= 0) {
			IJ.error("Filter scale must be greater than 0");
			gd = new GenericDialog ("DSSIM");
			out = false;
		}
		if (frame_offset >= stackSize) {
			IJ.error("Frame offset must be less than total frames in stack");
			gd = new GenericDialog ("DSSIM");
			out = false;
		}
		if (frame_offset <= 0) {
			IJ.error("Frame offset must be greater than 0");
			gd = new GenericDialog ("DSSIM");
			frame_offset = 1;
			out = false;
		}
		if (downsampled > downsampled_backup) {
			IJ.error("Miminum height must be 256 pixels (review Viewing scale)");
			gd = new GenericDialog ("DSSIM");
			out = false;
		}
		if (downsampled < 1) {
			IJ.error("Minimun value of Viewing scale must be 1");
			gd = new GenericDialog ("DSSIM");
			out = false;
		}
		gd.dispose();
	}
	
//
// NOW, WE CREATE THE FILTER, GAUSSIAN OR MEDIA FILTER, ACCORDING TO THE VALUE OF boolean "gaussian_window"
//For loop starts around here?
	filter_length = filter_width*filter_width;
	float window_weights [] = new float [filter_length];
	double [] array_gauss_window = new double [filter_length];

	if (gaussian_window) {
		double value, distance = 0;
		int center = (filter_width/2);
  		double total = 0;
		double sigma_sq=sigma_gauss*sigma_gauss;
		
      	  	for (int y = 0; y < filter_width; y++){
			for (int x = 0; x < filter_width; x++){
         				distance = Math.abs(x-center)*Math.abs(x-center)+Math.abs(y-center)*Math.abs(y-center);
				pointer = y*filter_width + x;
                			array_gauss_window[pointer] = Math.exp(-0.5*distance/sigma_sq);
				total = total + array_gauss_window[pointer];
  			}
    		}
		for (pointer=0; pointer < filter_length; pointer++) {	
			array_gauss_window[pointer] = array_gauss_window[pointer] / total;
			window_weights [pointer] = (float) array_gauss_window[pointer];
		}
	}
	else { 								// NO WEIGHTS. ALL THE PIXELS IN THE EVALUATION WINDOW HAVE THE SAME WEIGHT
		for (pointer=0; pointer < filter_length; pointer++) {
			array_gauss_window[pointer]= (double) 1.0/ filter_length;
			window_weights [pointer] = (float) array_gauss_window[pointer];
		}
	}

	
	/**if (show_gaussian_filter) {					// IN CASE OF A GAUSSIAN FILTER, YOU CAN SHOW IT IF YOU WANT	
		ColorModel cm=null;
		ImageProcessor gauss_window_ip = new FloatProcessor (filter_width, filter_width, window_weights, cm);
		gauss_window_ip = gauss_window_ip.resize (filter_width*filter_scale);
		String title_filtro_1 = "Sigma: " + sigma_gauss + " Width: "+ filter_width + " p�xeles"; 	
		ImagePlus gauss_window_imp = new ImagePlus (title_filtro_1, gauss_window_ip);
		gauss_window_imp.show();
		gauss_window_imp.updateAndDraw(); 		
	}**/
//
// END OF FILTER SELECTION							
//
//
// MAIN ALGORITHM
	
	ImageProcessor image_1_p, image_2_p;
	ImageStack stack = image_imp.getImageStack();
	//String message_1;
	//ssim_map = new double [image_dimension];
	ImageStack ssim_stack = new ImageStack((int) (image_width/downsampled), (int) (image_height/downsampled));
	double index_values[] = new double[stackSize-frame_offset];
	boolean swi = true;

	ImageProcessor image_1_original_p;
	ImageProcessor image_2_original_p;
	ImageProcessor mu1_ip;
	ImageProcessor mu2_ip;
	float [] array_mu1_ip;
	float [] array_mu2_ip;
	for (int k = 1; k <= stackSize-frame_offset; k++) { //for loop start

	image_1_original_p = stack.getProcessor(k);
	image_2_original_p = stack.getProcessor(k+frame_offset);
	//ImageProcessor image_1_original_p = image_imp.getProcessor();//stackslices instead returns ip
	//ImageProcessor image_2_original_p = image_imp.getProcessor();	
	//IJ.error("n "+String.valueOf(image_1_original_p.getf(0)));
	//IJ.error("b "+String.valueOf(image_2_original_p.getf(0)));

	image_width = image_1_original_p.getWidth();
	image_width = (int) (image_width/downsampled);//pesky downsampling
	image_1_original_p.setInterpolate(true);
	image_2_original_p.setInterpolate(true);
	image_1_p = image_1_original_p.resize (image_width);
	image_2_p = image_2_original_p.resize (image_width);

	image_height = image_1_p.getHeight();// stack set with width &height but i guess they can take th ip from the stack
	image_width = image_1_p.getWidth();
	image_dimension = image_width*image_height;

	mu1_ip = new FloatProcessor (image_width, image_height);
	mu2_ip = new FloatProcessor (image_width, image_height);
	array_mu1_ip = (float []) mu1_ip.getPixels();
	array_mu2_ip = (float []) mu2_ip.getPixels();

	float [] decheck = new float [image_dimension];
	double desum=0;
	float [] array_mu1_ip_copy = new float [image_dimension];
	float [] array_mu2_ip_copy = new float [image_dimension];
	//IJ.error(String.valueOf(array_mu1_ip[2]));
	a=b=0;
	for (pointer =0; pointer<image_dimension; pointer++) {	

		if (bits_per_pixel_1 == 8) {
			a = (0xff & image_1_p.get (pointer));
			b = (0xff & image_2_p.get(pointer));
			array_mu1_ip [pointer] = array_mu1_ip_copy [pointer] = a; // Float.intBitsToFloat(a);
			array_mu2_ip [pointer] = array_mu2_ip_copy [pointer] = b;
		}
		if (bits_per_pixel_1 == 16) {
			a = (0xffff & image_1_p.get(pointer));
			b = (0xffff & image_2_p.get(pointer));
			array_mu1_ip [pointer] = array_mu1_ip_copy [pointer] = a; // Float.intBitsToFloat(a);
			array_mu2_ip [pointer] = array_mu2_ip_copy [pointer] = b; // Float.intBitsToFloat(b);	
		}
		if (bits_per_pixel_1 == 32) {
			af = (image_1_p.getf(pointer));
			//IJ.error(String.valueOf(a));
			bf = (image_2_p.getf(pointer));
			//IJ.error("b "+String.valueOf(b));
			//break;
			/**String message_1= " ";
			String message_2 = "its 32b"; 
			IJ.showMessage (message_1, message_2);**/
			array_mu1_ip [pointer] = array_mu1_ip_copy [pointer] = af;
			array_mu2_ip [pointer] = array_mu2_ip_copy [pointer] = bf;
			if (swi==true) {
				//decheck[pointer]=af-floor(af);
				desum+=af-Math.floor(af);
			}
		}
		
	}
	swi=false;

	//IJ.error("desum "+String.valueOf(desum));
	//IJ.error("Val "+String.valueOf(image_2_p.getf(0)));
	double C1 = K1;
	C1= C1*C1;
	double C2;
	if (bits_per_pixel_1 == 32) {
		if (desum==0){
			//IJ.error("desum=0 ");
			C2 = (Math.pow(2, bits_per_pixel_1) - 1)*K2;
		} else {
			C2 = K2;
		}
	} else {
		C2 = (Math.pow(2, bits_per_pixel_1) - 1)*K2;
	}

	C2=C2*C2;
	double C3 = C2/2;

	//float[] arr = (float[])image_1_p.getPixels();
	//IJ.error("greyscale: "+String.valueOf(arr[0]));
	//IJ.error("n "+String.valueOf(array_mu1_ip[0]));
	//IJ.error(String.valueOf(array_mu1_ip[0]));
	//IJ.error("x "+String.valueOf(mu1_ip.getf(0)));
	mu1_ip.convolve (window_weights, filter_width, filter_width);
	mu2_ip.convolve (window_weights, filter_width, filter_width);
	//IJ.error("z "+String.valueOf(mu1_ip.getf(0)));
	double [] mu1_sq = new double [image_dimension];
	double [] mu2_sq = new double [image_dimension];
	double [] mu1_mu2 = new double [image_dimension];

	for (pointer =0; pointer<image_dimension; pointer++) {
		mu1_sq[pointer] = (double) (array_mu1_ip [pointer]*array_mu1_ip [pointer]);
		mu2_sq[pointer] = (double) (array_mu2_ip[pointer]*array_mu2_ip[pointer]);
		mu1_mu2 [pointer]= (double) (array_mu1_ip [pointer]*array_mu2_ip[pointer]);
	}

	double [] sigma1_sq = new double [image_dimension];
	double [] sigma2_sq = new double [image_dimension];
	double [] sigma12 = new double [image_dimension];
	double [] sigma1_sigma2 = new double [image_dimension];

	for (pointer =0; pointer<image_dimension; pointer++) {
			
		sigma1_sq[pointer] =(double) (array_mu1_ip_copy [pointer]*array_mu1_ip_copy [pointer]);
		sigma2_sq[pointer] =(double) (array_mu2_ip_copy [pointer]*array_mu2_ip_copy [pointer]);
		sigma12 [pointer] =(double) (array_mu1_ip_copy [pointer]*array_mu2_ip_copy [pointer]);
		sigma1_sigma2 [pointer] =(double) (array_mu1_ip_copy [pointer]*array_mu2_ip_copy [pointer]);
	}
//	
//THERE IS A METHOD IN IMAGEJ THAT CONVOLVES ANY ARRAY, BUT IT ONLY WORKS WITH IMAGE PROCESSORS. THIS IS THE REASON BECAUSE I CREATE THE FOLLOWING PROCESSORS
//
	ImageProcessor soporte_1_ip = new FloatProcessor (image_width, image_height);
	ImageProcessor soporte_2_ip = new FloatProcessor (image_width, image_height);
	ImageProcessor soporte_3_ip = new FloatProcessor (image_width, image_height);
	float [] array_soporte_1 =  (float []) soporte_1_ip.getPixels();
	float [] array_soporte_2 =  (float []) soporte_2_ip.getPixels();
	float [] array_soporte_3 =  (float []) soporte_3_ip.getPixels();

	for (pointer =0; pointer<image_dimension; pointer++) {
		array_soporte_1[pointer] = (float) sigma1_sq[pointer];
		array_soporte_2[pointer] = (float) sigma2_sq[pointer];
		array_soporte_3[pointer] = (float) sigma12[pointer];
	}

	soporte_1_ip.convolve (window_weights, filter_width,  filter_width);
	soporte_2_ip.convolve (window_weights, filter_width,  filter_width); 
	soporte_3_ip.convolve (window_weights, filter_width,  filter_width);
	//IJ.error("sopo "+String.valueOf(soporte_3_ip.getf(0)));

	for (pointer =0; pointer<image_dimension; pointer++) {
		sigma1_sq[pointer] =  array_soporte_1[pointer] - mu1_sq[pointer];
		sigma2_sq[pointer] =  array_soporte_2[pointer ]- mu2_sq[pointer];
		sigma12[pointer] =  array_soporte_3[pointer] - mu1_mu2[pointer];
	}
	ssim_map = new double [image_dimension];
	double suma=0;
	
	/**
	IJ.error("mu1_mu2"+String.valueOf(mu1_mu2[0]));
	IJ.error("mu1_sq"+String.valueOf(mu1_sq[0]));
	IJ.error("mu2_sq"+String.valueOf(mu2_sq[0]));
	IJ.error("c1 "+String.valueOf(C1));
	IJ.error("c2 "+String.valueOf(C2));
	IJ.error("c3 "+String.valueOf(C3));
	**/
	for (pointer =0; pointer<image_dimension; pointer++) {
		//ssim_map[pointer] = (double) (( 2*mu1_mu2[pointer] + C1)* (2*sigma12[pointer] + C2)) / ((mu1_sq[pointer]+mu2_sq[pointer] + C1) * (sigma1_sq[pointer] + sigma2_sq[pointer] + C2));
		ssim_map[pointer] = (double) (Math.pow((2*mu1_mu2[pointer]+C1)/(mu1_sq[pointer]+mu2_sq[pointer]+C1), alpha)*
				Math.pow((2*sigma1_sigma2[pointer]+C2)/(sigma1_sq[pointer] + sigma2_sq[pointer] + C2), beta)*
				Math.pow((sigma12[pointer]+C3)/(sigma1_sigma2[pointer] + C3), gamma));
		//IJ.error("ssim "+String.valueOf(ssim_map[0]));
				ssim_map[pointer] = (1-ssim_map[pointer])/2;
		suma = suma + ssim_map[pointer];
		
	}
	//put if statement
	if (per_frame==true){
		double [] sorted_dssim = ssim_map.clone();
		Arrays.sort(sorted_dssim);
		int perc = (int)(image_dimension*0.01);
		double ssi_low=sorted_dssim[perc];
		double ssi_high=sorted_dssim[image_dimension-perc-1];
		for (int z = 0; z<image_dimension; z++){
			if (ssim_map[z]<ssi_low){
				ssim_map[z]=ssi_low;
			} else if (ssim_map[z]>ssi_high){
				ssim_map[z]=ssi_high;
			}
			ssim_map[z]=(ssim_map[z]-ssi_low)/(ssi_high-ssi_low);
		}
	}
	//IJ.error("s"+String.valueOf(ssim_map[1]));
	ssim_index = (double) suma / image_dimension;
	//message_1= " ";
	index_values[k-1] = ssim_index;
	if (show_ssim_map) {
		ImageProcessor ssim_map_ip = new FloatProcessor (image_width, image_height, ssim_map);
		ImagePlus ssim_map_imp = new ImagePlus (" ", ssim_map_ip);
		ssim_stack.addSlice(ssim_map_imp.getProcessor());
	}
	/**String message_2 = "ssim_index:  " + ssim_index; 
	message_1= " ";
	IJ.showMessage (message_1, message_2);**/
	} //end of for loop
	if (save_values == true){
		try { //print to file
			PrintWriter writer = new PrintWriter(new FileWriter(directory+"values.csv", true));
			//writer.println(index_values.toString());
			for (int i =0; i<index_values.length;i++) {
				if (i>0) {
					writer.print("\n"+String.valueOf(index_values[i]));
				} else {
					writer.print(String.valueOf(index_values[i]));
				}
			}
			writer.close();
			
		}catch (Exception e) {
			IJ.error("No Stack");
		}
	}

	if (show_ssim_map) {
		ImagePlus ssim_draw = new ImagePlus("DSSIM_Stack", ssim_stack);
		ssim_draw.show();
		ssim_draw.updateAndDraw();
		/**if (use_color_map == true){
			LutLoader loader = new LutLoader();
			loader.run(lut_path);
		}**/
	}



	//message_1= " ";
	//String message_2 = "ssim_index:  " + ssim_index; 
	IJ.showProgress(1.0);
	//IJ.showMessage (message_1, message_2);
}
}
