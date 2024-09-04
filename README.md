# Altered Brain Network Topology in Children with Auditory Processing Disorder

## Introduction
This study explores brain network organization in children with Auditory Processing Disorder (APD) using resting-state fMRI data and graph theory approaches. Key findings include differences in brain hub architecture and functional connectivity between children with APD and healthy controls, particularly in regions related to auditory processing.

## Background
**Problem Statement:** Auditory Processing Disorder (APD) leads to difficulties in understanding speech despite normal hearing. The origins of APD symptoms are debated, with limited knowledge on the role of altered brain network topology.

**Objectives:** To investigate the functional brain network organization in children with APD and compare it with healthy controls using advanced neuroimaging techniques and network science approaches.

## Methods

### Population for collecting multi-modal data
- 66 children (57 included) aged 8-14 years old (28 with diagnosis of APD and 29 healthy controls).

### Data type
- functional MRI (rs-fMRI) data acquired by multi-echo multi-band imaging sequences

### Analysis Techniques
- Complex network analysis using graph theory, focusing on global and nodal network properties, and brain hub architecture

<!-- ![Alt text](img/1.jpg) -->

<p align="center">
<!-- <img src="img/1.jpg" alt="Alt text" width="700" height="600"/> -->

<!-- for compatibility with screen sizes -->
<img src="img/1.jpg" alt="Description" style="max-width:100%; height:auto;"> 
</p>


<p align="center">
Figure 1: The pipeline used in this project for analyzing the imaging data.
</p>

## Results

### Global Network Properties
- Both APD and control groups showed similar global network properties, but differences emerged in hub architecture

<p align="center">
  <img src="img/2.jpg" alt="Image 1" style="max-width:100%; height:auto;"/>
  <img src="img/3.jpg" alt="Image 2" style="max-width:100%; height:auto;"/>
</p>
Figure 2: Similar brain hub architecture and <a href="https://en.wikipedia.org/wiki/Large-scale_brain_network">intrinsic network model</a> in both APD and HC (Left) as shown by their global network metrics (right).


### Key Findings
- Decreased participation coefficient in auditory cortical regions (<a href="https://www.sciencedirect.com/topics/neuroscience/superior-temporal-gyrus#:~:text=The%20superior%20temporal%20gyrus%20is,short%2Dterm%20auditory%20sensory%20memory.">bilateral superior temporal gyrus</a>) in children with APD, suggesting altered functional connectivity in specific brain networks (<a href="#target-image1">Figure 3: Left</a>).
- Positive correlation between left  <a href="https://www.sciencedirect.com/topics/neuroscience/parahippocampal-gyrus">parahippocampal gyrus </a> connectivity and auditory perception tasks in children with APD (<a href="#target-image2">Figure 3: Right</a>).

<p align="center">
  <img src="img/5.jpg" alt="Image 1" width="400" id="target-image1"/>
  <img src="img/4.jpg" alt="Image 2" width="400" id="target-image2" />
</p>

Figure 3: Major differences in brain regional network in APD and HC based on network metrics.

## Discussion
The findings suggest that children with APD have distinct alterations in brain network topology, particularly in regions associated with auditory processing. These results provide new insights into the neural mechanisms underlying APD and highlight the potential for using advanced neuroimaging and network science approaches to further our understanding of auditory processing disorders.

## Conclusion
This study contributes to the understanding of brain network alterations in APD, with potential implications for developing targeted interventions.

## Publication and Citation
**Read the article:** [here](https://www.sciencedirect.com/science/article/pii/S2213158222002042)

**Citation:** Alvand, A., et al. (2022). Altered brain network topology in children with auditory processing disorder: A resting-state multi-echo fMRI study. *NeuroImage: Clinical, 35*, 103139.

## Funding
This study was funded by [Eisdell Moore Centre](https://www.emcentre.ac.nz/) and Faculty of Science's Research fund from [the University of Auckland](https://www.auckland.ac.nz/en.html).

# Data science behind this study

## Data Organization steps

### Data conversion

### Data restructuring

## Data quality assessment

### Data initial pre-processing (minimal cleaning)

### Selection of optimal de-noising pipelines (benchmarking strategies)

## Modeling the data 

### Constructing network model

#### Defining data points

#### Defining connections

## Processing the data (Application of network algorithms)

### Global and local topology evaluation

### Backbone consistency tests (Hub model)

### Connection consistency test

## Statistics 

### Multivariate tests

### Correlation tests





