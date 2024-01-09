# Ants integrate proprioception, visual context and efference copies to make robust predictions

## Description

This repository contains the data and codes necessary to obtain the results from the article "Ants integrate proprioception, visual context and efference copies to make robust predictions" by O. Dauzere-Peres & A. Wystrach (2023).
The raw data files contained in this repository were collected by a virtual reality set-up running with Unity 2020.1 envrironment controlling the rotation of the scene projected on an LED screen around an ant mounted on a trackball. Raw data consists of the movements of the ants on the ball, which were obtained by two sensors quantifying the movement of the ball along the three axis.

The python code files are used to process the raw data in order to visualize them in the way it is shown in the articles figures, and to create the csv files with the variables used for statistical analysis. The cvs files produced are also uploaded in this repository.

## Installation

To reproduce the results of the article you will need to download the files containing the raw data corresponding to the experiments you want to study and to install python (https://www.python.org/downloads/) to be able to execute the codes in this repository. 

Python packages necessary to run the python codes in this repository are specified at the beginning of each file.

## Usage

The file named "extract_raw_data.py" should be executed first. This file contains 4 sub-sections each corresponding to an experiment of the study (3 from the main study and the last one in supplementary). Each section should be executed after you downloaded the raw data corresponding to the associated experiment and you write the path to those files in the line instructed.
This will compute the variables you will need to visualize the ants movements with the other python codes files as well as produce the csv files used for statistical analysis. If you want to check and compare with the csv files used in the study for statistical analysis, they are also uploaded in the repository.
This also means that the sub-section corresponding to the experiment you want the visualize always need to be the last one you executed so that the global variables correspond to the one of this experiment. 

To look at the angular velocity signals or trajectories of individual ants in differents experimental conditions you need the "trajectories_and_signal.py" file and to execute the corresponding function. 
To obtain the mean oscillation cycles for an individual ant or at the population level in given experimental conditions you need the "mean_cycle.py" file and execute the corresponding function. Examples to reproduce the mean oscillation cycles shown in the article are given at the end of the file.

## License

The MIT License (MIT)
Copyright (c) 2024 Océane Dauzere-Peres

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
