/*RingIntensityDistribution.ijm
Fiji Lifeline 30 May 2017
E.Rebollo, Molecular Imaging Platform IBMB

This macro divides the cell into rings of identical area that converge towards the
nucleus center. The intensity density of the interest signal is then measured per 
ring and normalized to the total intensity density of the cell. 

Open 1) an image containing the signal of interest, and 2) A RoiSet.zip containing 
the two regions of interest, the cell and the nucleus boundaries, in this order. 
Then hit Run...*/

//RETRIEVE IMAGE ID
image = getImageID();

//CREATE DISTANCE MAP
distanceMap();
//distanceMap() creates a distance map within the cell ROI (first ROI in the ROI Manager)
function distanceMap() {
	newImage("distance map", "8-bit black", 2048, 2048, 1);
	roiManager("select",0);
	setForegroundColor(255,255,255);
	run("Fill");
	run("Make Binary");
	run("Distance Map");
	roiManager("Show None");

}

//FILL ARRAYS WITH ROI COORDINATES
//Get nucleus (ROI02) centroid and fill into array
roiManager("Show None");
roiManager("select",1);
run("Set Measurements...", "area centroid redirect=None decimal=2");
run("Measure");
XC=getResult("X");
YC=getResult("Y");
//Get cytoplasmic (ROI01) coordinates and fill into array
roiManager("select",0);
Roi.getCoordinates(X,Y);
kX=newArray(X.length);
kY=newArray(Y.length);

//CORRECT DISTANCE MAP TO SWIFT CENTER TOWARDS THE NUCLEUS CENTROID
selectWindow("distance map");
newMap(kX,kY,XC,YC);
/*newMap() shifts the distance map values towards the centroid of the nucleus (second ROI
in the ROI Manager); 
it first calculates the max distance from the nucleus centroid to the periferal ROI; 
then, using Roi.getContainedPoints, the X and Y coordinates of the 
pixels inside the selection are obtained in two arrays X1 and X2; the distance map is 
converted to 16bits to allow for future calculations > 255; a For loop is used to modify 
the pixel values and correct them according to the maximum distance (dmax)/d;  */

function newMap(Xout,Yout,XC,YC) {
	dmax=0;
	for (j=0;j<Xout.length;j++) {
		dX=Xout[j]-XC;
		dY=Yout[j]-YC;
		d=sqrt(pow(dX,2)+pow(dY,2));
		if (d>dmax) {
			dmax=d;
			}
		}
	roiManager("select",0);
	Roi.getContainedPoints(X1,Y1);
	run("16-bit");
	selectWindow("distance map");
	for (i=0;i<X1.length;i++) {
	dX=X1[i]-XC;
	dY=Y1[i]-YC;
	d=sqrt(pow(dX,2)+pow(dY,2));
	v=getPixel(X1[i],Y1[i]);
	newValue=v*(dmax)/d; 
	setPixel(X1[i],Y1[i],newValue);
}
}

//DELETE NUCELAR ROI (ROI02)
roiManager("select", 1);
roiManager("delete");

//CREATE CELL RINGS
//Threshold the distance map and create selections according to the desired area %
thresholds=newArray(0.75,0.5,0.25); 
n=1;//this is the counter for the function shrinkROI()
for (i=0;i<(thresholds.length);i++){
shrinkROI(i,thresholds);
}

/* shrinkROI() thresholds the distance map image using a second function iterativeDistanceMapThresholding();
It measures the area of the current ROI and applies the second function until an area value is reached*/
function shrinkROI(i,thresholds) {
	setBatchMode(true);
	run("Set Measurements...", "area perimeter shape redirect=None decimal=2");
	roiManager("select",0);
	roiManager("Measure");
	area=getResult("Area");
	roiManager("add");
	roiManager("select",1+i);
	areaThr=thresholds[i]*area;
	while (area>areaThr) {
		iterativeDistanceMapThresholding(n,i);
		roiManager("select",1+i);
		roiManager("Measure");
		area=getResult("Area");
		n++;
	}
	setBatchMode(false);
}

/* iterativeDistanceMapThresholding() creates a threshold on the distance map image (the one actually selected)
 from n+1 (value 2) to 65535. Then converts it into a selection that replaces the previous selection in the 
 ROI Manager */
function iterativeDistanceMapThresholding(n,i){ 
	setBatchMode(true);
	roiManager("Show None");
    setThreshold(n+1, 65535);
    run("Create Selection");
    roiManager("add"); 
    roiManager("select",1+i);
    roiManager("delete");
	setBatchMode(false);
}

//Close distance map
selectWindow("distance map");
run("Close");
run("Clear Results");

//INTENSITY CALCULATIONS
//measure intensity density in each ROI
run("Set Measurements...", "area mean integrated redirect=None decimal=2");
roiManager("deselect");
roiManager("measure");
//Create arrays to store measurements
noRois = roiManager("count"); //hay que actualizar el número de ROIS
myDensities=newArray(noRois);
myAreas=newArray(noRois);
myIntensityPercentages=newArray(noRois);
myMeans=newArray(noRois);
//Fill arrays
for(j=0; j<noRois; j++){
	myDensities[j]=getResult("IntDen",j);
	myAreas[j]=getResult("Area",j);
}	

//Store Total cell intensity Density into a variable for later calculations
totalIntDen=myDensities[0];

//Calculate ring IntDen and Areas by subtraction of subsequent ROIs
ringsIntDen=newArray(noRois-1);
ringsAreas=newArray(noRois-1);
for(i=0; i<noRois-1; i++){
	ringsIntDen[i]=myDensities[i]-myDensities[i+1];
	ringsAreas[i]=myAreas[i]-myAreas[i+1];
}

//measure values from central ring
noRois = roiManager("count"); 
roiManager("select", noRois-1);
centralIntDen=getResult("IntDen");
centralArea=getResult("Area");

//Complete final arrays 
allIntDen=Array.concat(ringsIntDen, centralIntDen);
allAreas=Array.concat(ringsAreas, centralArea);

//Fill final arrays with mean intensity values and intensity % per Ring 
for(i=0; i<noRois; i++){
	myMeans[i]=allIntDen[i]/allAreas[i];
	myIntensityPercentages[i]=100*allIntDen[i]/totalIntDen;
}
//Array.show(myMeans);
Array.show(myMeans);
Array.show(myIntensityPercentages);

//create verification image
selectImage(image);
run("RGB Color");
setForegroundColor(255, 255, 0);
roiManager("deselect");
roiManager("Set Color", "yellow");
roiManager("Set Line Width", 2);
roiManager("draw");

//Close windows
selectWindow("ROI Manager");
run("Close");
selectWindow("Results");
run("Close");



////////////// 
