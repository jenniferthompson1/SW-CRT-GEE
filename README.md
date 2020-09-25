# Comparison of small-sample standard-error corrections for generalised estimating equations in stepped wedge cluster randomised trials with a binary outcome: A simulation study

This repository contains all programs used to conduct the analysis for our paper "Comparison of small-sample standard-error corrections for generalised estimating equations in stepped wedge cluster randomised trials with a binary outcome: A simulation study". The paper can be read here: https://doi.org/10.1177/0962280220958735

## Description

Details are contained in our paper found here: https://doi.org/10.1177/0962280220958735. 

In brief, the files in DataGeneration sample from a multivariate normal distribution to simulate data from stepped wedge cluster randomised trials with particular set characteristics. For each each simulated trial dataset, several  generalised estimating equations analysis are conducted with small sample corrections. In phase one these 

- Standard error corrections
  - Fay and Graubard 
  - Kauermann and Carroll
  - MacKinnon and White
  - Mancl and DeRouen
  - Morel, Bokossa, and Neerchal
- Degree of freedom corrections
  - Clusters minus parameters
  - Cluster periods minus parameters
  - Cluster periods minus parameters including clusters as parameters
  - Fay and Graubard estimate
  - Pan and Wall estimate

Analysis files then summarise the performance of these analysis methods for each set of characteristics. Performance measures used include

- Intervention effect estimate bias
- Standard error bias
- Mean and standard error of degree of freedom estimators
- Confidence interval coverage
- Power

## Getting Started

### Dependencies

Data generation was run on a Scientific Linux 7 server with R version 3.5.3. Th analysis programs were run on Windows 10 with R version 3.6.1

### Installing

These files have been modified from the versions used to conduct the original analysis. Users may need to modify file paths.

### Executing programs

Master.R contains a suggested order to run the analysis files

## Help

Users should contact jennifer.thompson@lshtm.ac.uk if there are problems downloading and running the programs. 

## Authors

Jennifer Thompson jennifer.thompson@lshtm.ac.uk

## License

This project is licensed under the CC BY 4.0 license.
