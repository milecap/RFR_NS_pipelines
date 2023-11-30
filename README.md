# RESL-SIRS postprocessing analysis

This code generates the figures and results presented in the publication:
Stimulus-Induced Rotary Saturation imaging of visually evoked neuroelectric response: preliminary results and data analysis (currently under review)

## Visual-stimulation Dataset
The dataset used for the study can be accessed through the link: https://mb-neuro.medical-blocks.ch/shared/folder/4947b990-8f5a-11ee-aeab-37c16cf6fc6b

Note: After registering, the data repository will allow you to download the dataset. In reason for download please state:
"Dataset for the paper: Stimulus-Induced Rotary Saturation imaging of visually evoked neuroelectric response: preliminary results and data analysis"
Due to our data sharing policies, the high-resolution images cannot be made available, instead, the low-resolution coregistered data is shared.

## Requirements
- MATLAB R2021a or later

All the functions used to run the pipeline are included in the shared folder or are available Matlab toolboxes. 

The violin_2 function is a slightly modified version of the originally published by 
% Hoffmann H, 2015: violin_2.m - Simple violin_2 plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany. 


## Installation
1. Clone the repository to your local machine:
   ```bash
   git clone https://github.com/milecap/RFR_NS_pipelines.git


## Usage
To run the RFR pipeline run the RFR_prepare.m
To run the NS pipeline run the NS_prepare.m
For both cases, you might change the cutoff of the highpass filter (default 0.1 Hz).


## Output
Both the RFR and NS file will create an output folder inside each Sub_n folder. 
Additionaly, running the visualization or statistical analysis will create an output folder in the data folder that sumarises the statistical outputs for all subjects.

## Contributing
Please feel free to let me know if you find a bug or you have suggestions to improve the code/analysis!



