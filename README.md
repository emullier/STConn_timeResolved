# STConn_timeResolved


This repository gathers the codes related to the manuscript: </p>
<b> 'Functional brain dynamics are shaped by connectome n-simplicial organization' </b>  
<i> Emeline Mullier, Jakub Vohryzek, Alessandra Griffa, Yasser Alemàn-Gómez, Celia Hacker, Kathryn Hess, Patric Hagmann </i>  (Submitted)    

The dataset used for the analysis in this manuscript is publicly available under a Zenodo repository: 
<a href="url">10.5281/zenodo.3576326</a>. On Zenodo, you can find the preprocessed functional MRI timeseries of 113 healthy controls, a mean structural connectome and the atlas properties (More information can be found in the corresponding repository). Here we provide 3 examples of these subjects processed with the spatio-temporal connectome (STConn) [Griffa2015] <a href="url">https://github.com/agriffa/STConn</a>  


The codes are matlab codes organized into 3 main folders of analysis and 3 adjunct folders.  </p>
The 3 main folders are:  </p>
- GenerationStates: script and functions to generate the states
- DynamicMeasures: script and functions to estimate the dynamic measures and plot the circular graph
- CompStates: script and function to compare the states the seven yeo functional systems. </p>

And the 3 adjunct folders are:  </p>
- DataExample: outputs of the STConn for three example subjects from the dataset.
- BCT: Brain Connectivity Toolbox [Rubinov2010] <a href="url"> https://sites.google.com/site/bctnet/</a>  
- PlotKit: set of codes for brain plotting and results visualization </p>



#### References
<font size="1"> [Griffa2015] 'Transient networks of spatio-temporal connectivity map communication pathways in brain functional systems' - Alessandra Griffa, Benjamin Ricaud, Kirell Benzi, Xavier Bresson, Alessandro Daducci Pierre Vandergheynst, Jean-Philippe Thiran, Patric Hagmann - NeuroImage. Volume 155, 15 July 2017, Pages 490-502 </font>


<font size="1"> [Rubinov2010] 'Complex network measures of brain connectivity: Uses and interpretations.' - Mikail Rubinov, Olaf Sporns - NeuroImage 52:1059-69 </font>


Copyright (c) 2019 Lausanne University Hospital (UNIL-CHUV)
