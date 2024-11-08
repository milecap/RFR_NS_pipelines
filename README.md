# RESL-SIRS postprocessing analysis

This code generates the figures and results presented in the publication:
Capiglioni M, Beisteiner R, Cardoso PL, et al. Stimulus-induced rotary saturation imaging of visually evoked response: A pilot study. NMR in Biomedicine. 2024;e5280. doi:10.1002/nbm.5280

## Visual-stimulation Dataset
The dataset used for the study is available upon reasonable request due to privacy/ethical restrictions.

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



