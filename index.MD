# fMRI Chronic Pain Meta-analysis Documentation

Written by Anna Xu | Replicator: Alina Henn

Last Updated: May 12, 2020

## Setup
For these analyses, you will need the following software:
* [MATLAB](https://www.mathworks.com/products/matlab.html)
* [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
* [BrainNet Viewer](https://www.nitrc.org/projects/bnv/)
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL)

To begin, download [**this GitHub repo**](https://github.com/PennBBL/Xu_fMRIChronicPain). Next, move the **spm12** and **BrainNet Viewer** folders into the `dependencies` folder. After you're done with these steps, open the script [**runALE.m**](https://github.com/PennBBL/Xu_fMRIChronicPain/blob/master/runALE.m) on MATLAB. This will be the script to run analyses. In the **ale** folder, create two subdirectories: **DataMatlab** and **DataRaw**. This is where inputted data will be saved later on.

The first line, `curPath=genpath(pwd); addpath(curPath);` will add all dependencies into your main directory. The next line, `cd ale` will put you in the [**ale**](https://github.com/PennLINC/Xu_fMRIChronicPain/tree/master/ale) folder so that you can begin your analyses.

### For Windows OS

Please note that these scripts were ran using Mac OS (Mojave 10.14.1). Prior to running these scripts, if you are using Windows OS, you will need to do the following:

1. Download Visual Studio for C++: https://visualstudio.microsoft.com/vs/features/cplusplus/

2. Download additional documents specific to Windows in the [**Windows_Specific**]() folder. Move these files into your [**ale**] folder.

3. Run the following commands in MATLAB beforehand:
```
mex -setup
mex 'tfceMex.c'
mex 'tfceMex_pthread.c'
```
You are now all set for running the ALE scripts! One other caveat with Windows is you will need a way to unzip **.gz** files. We recommend the program [**WinRAR**](https://www.chip.de/downloads/WinRAR-64-Bit_37666786.html) for this.

## Conduct Meta-analysis

### Main Analyses

#### Read data

To run the analyses, first make sure the coordinate data is read in to the spreadsheet by running the line `ale_inputCoords('data/aberrant_20200113.xls')`. This reads the coordinate data for the main analyses (aberrant activity) located in [**data/aberrant_20200113.xlsx**](https://github.com/PennLINC/Xu_fMRIChronicPain/blob/master/ale/data/aberrant_20200113.xls) and converts it into a **.mat** format for analyses. Next, run the line `ale_inputCoords('data/painCoords_20200113.xls');` to read the coordinate data for subset analyses of signed differences. This coordinate data is located in [**data/painCoord_20200113.xls**](https://github.com/PennLINC/Xu_fMRIChronicPain/blob/master/ale/data/painCoords_20200113.xls). These coordinate data will now be read in **.mat** format for analyses and are located in **ale/DataMatlab**.

#### ROI-based analyses set up

To conduct ROI-based analyses, fist download the Xu et al., 2020 meta-analytic map from Neurovault as Volume: https://neurovault.org/collections/AJKXTUXX/images/359250/. On your terminal, binarize this image using the line `fslmaths "pain_thresh_cFWE05.nii.gz" -bin "healthyPainBin.nii.gz"`. This outputs a binarized nifti named **healthyPainBin.nii.gz**. Unzip this file and move **healthyPainBin.nii** into the [**ale**](https://github.com/PennBBL/Xu_fMRIChronicPain/tree/master/ale) folder.

#### Run analyses

After, run the line `ale_inputCoords(data/painCoords_20200113.xls)` This reads in the data for signed differences between patients and controls (i.e., separate sub-analyses of experiments using "patients > controls" contrast and experiments using "controls > patients" contrast), located in [**data/painCoords_20200113.xlsx**](https://github.com/PennBBL/Xu_fMRIChronicPain/blob/master/data/painCoords_20200113.xls).

Next, run the line `ale_estimateALE('contrasts/pain_20200113.xlsx');
` to run all analyses.  This line reads the spreadsheet [**contrasts/pain_20200113.xlsx**](https://github.com/PennBBL/Xu_fMRIChronicPain/blob/master/contrasts/pain_20200113.xlsx) and runs the contrasts outlined there. ROI-based analyses will also be conducted with this line.

The total runtime on the computer that ran the analyses locally (Macbook Pro 2018) is around 7 hours. That is the typical expected time on newer laptops.

### Perceptual Ratings Subset

To run sub-analyses looking only at experiments matching pain stimuli between patients and controls with similar painfulness ratings, run the line `ale_inputCoords('data/perceptualData_2020-01-21.xls');` to read the relevant coordinate data.

Next, run the line `ale_estimateALE('contrasts/perceptual_20200121.xlsx');` to run the relevant [**contrasts**](https://github.com/PennBBL/Xu_fMRIChronicPain/blob/master/contrasts/perceptual_20200121.xlsx).

## Interpret Results

The results are now located in a newly generated folder **ale/ALE**. Inside this folder, there is a folder named **Results**. This folder contains the thresholded meta-analytic maps. There is another folder called **Contribution**, which contains any significant coordinate results and % contribution from each experiment. Z-scored unthresholded are located in **VolumesZ** while foci images are located in **Foci**. Results from ROI-based analyses are located in **ROIresults**. These results are saved as **.pdf** files.

### Generate figures

All surface figures were generated using BrainNet Viewer. To use BrainNet, simply type `BrainNet` into the command line in MATLAB.

Unthresholded maps were loaded as the mapping image while **BrainNet_ICBM152.nv** was loaded as the surface file. Images were thresholded to positive values only for display. This is done by going to BrainNet_Option > Volume > Display: Positive only.

Axial view of images were generated using fslview. Images were thresholded to positive values only.

### Distribution of foci
Surface figures of foci distribution were generated using BrainNet. Before creating these figures, you need a **.node** file for BrainNet to read in. This can be generated by creating a **.csv** file in the main folder and running the [**tal2mni_coordinates_20200131.m**](https://github.com/PennBBL/Xu_fMRIChronicPain/blob/master/tal2mni_coordinates_20200131.m) script. This script reads the file [**coordinates_20200131.csv**](https://github.com/PennBBL/Xu_fMRIChronicPain/blob/master/coordinates_20200131.csv) file, which is just the coordinates used in the meta-analysis copied as a separate spreadsheet. Next, it converts any Tal coordinates to MNI and then outputs a **.csv** file with all coordinates in MNI.

Once you have the outputted spreadsheet, create three new columns next to the MNI coordinates. In the first new column, type `1` for all coordinates that are "patient > control" contrast and type `0` for all coordinates that are "control > patient" contrast (this information is under the header "contrast"). In the second new column, type 1. In the third new column, type `-`.

After you're done editing the spreadsheet, copy and paste the MNI coordinates under the columns `mni_x`, `mni_y`, and `mni_z` into a text editor as well as the three new columns generated. If you use TextEditor in MacOS, enter `command+shift+t` first to make it unformatted text editing. Delete the column header (i.e., the first row that contains the column names `mni_x`, `mni_y`, and `mni_z`), and then save this file as **coordinates.node**.

Open `BrainNet` in MATLAB. Next, choose **BrainNet_ICBM152.nv** as the surface figure and choose **coordinates.node** for Data file (node). Under "Node" in Options, go to "Color" and choose the option "Modular". In "Modular", click "More" and select the color red as Module1 and blue as Module2. Click "OK" and generate the figure.
