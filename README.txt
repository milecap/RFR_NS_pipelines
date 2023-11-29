# RESL-SIRS postprocessing analysis

The shared code can be used to generate the figures of the publication:
Stimulus-Induced Rotary Saturation imaging of visually evoked neuroelectric response: preliminary results and data analysis



## Table of Contents
- [Requirements]
- [Installation]
- [Usage]
- [Folder Structure]
- [Contributing]
- [License]

## Requirements
- MATLAB R2021a or later

All the functions used to run the pipeline are included in the shared folder or are available matlab toolboxes. 
The violin_2 function is a slightly modified version of the originally published by 
% Hoffmann H, 2015: violin_2.m - Simple violin_2 plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany. 


## Installation
Provide step-by-step instructions for installing and setting up the project.

1. Clone the repository to your local machine:
   ```bash
   git clone https://github.com/milecap/RFR_NS_pipelines


## Usage
To run the RFR pipeline run the RFR_prepare.m
To run the NS pipeline run the NS_prepare.m
For both cases, you might change the cutoff of the highpass filter (default 0.1 Hz).

## Folder Structure
The shared data structure has the form
Subjects_folder
-> 	Sub_n
	->	data
		->	noStim
			->	SLoff
				->	anat
				->	func
			->	SLon
				->	anat
				->	func
		->	VisStim
			->	SLoff
				->	func
			->	SLon
				->	func


## Contributing
Please feel free to let me know if you find a bug or suggestions to improve the code/analysis!

## License


