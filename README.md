# Sun_Ranking
This repository contains MATLAB code for Conjugate Bayesian Online Ranking for the Thurstone-Mosteller Model using Unified Skew-Normal Distributions. The paper associated with this code is provided in the 'pdf' folder. The code was developed and tested on an AMD Ryzen 3700X processor using MATLAB R2020A.

# How to reproduce results in our paper:
1. Navigate to the 'codes/sun_ranking_simu' folder. Execute the 'sun_ranking_active_simu.m' function to obtain results from the simulation of the SUN ranking active algorithm.
2. Navigate to the 'codes/sun_ranking_mfvi_simu' folder. Run the 'sun_ranking_active_mfvi.m' function to generate results from the simulation of the SUN ranking active algorithm with MFVI.
3. Go to the 'codes/reading_difficulty' folder. Use the 'sun_ranking_active_mfvi_reading.m' function to produce results of the SUN ranking active MFVI algorithm applied to the Reading Difficulty dataset.
4. Enter the 'codes/vqa' folder. Execute the 'sun_ranking_active_mfvi_vqa.m' function to generate results for the SUN ranking active MFVI algorithm applied to the Visual Question Answering (VQA) dataset.
5. Navigate to the 'codes/iqa' folder. Run the 'sun_ranking_active_mfvi_iqa.m' function to obtain results for the SUN ranking active MFVI algorithm applied to the Image Quality Assessment (IQA) dataset.

