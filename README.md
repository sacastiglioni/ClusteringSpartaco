# ClusteringSpartaco
by Sara Agavnì Castiglioni and Andrea Sottosanti.

In this package, the base and penalized clustering versions of the Spartaco model by A. Sottosanti and D. Risso are implemented.
For more details, please refer to the article on the Spartaco method, Co-clustering of Spatially Resolved Transcriptomic Data, by A. Sottosanti and D. Risso (2023), available at https://arxiv.org/abs/2110.04872, and to the Master's Thesis by Sara Agavnì Castiglioni in Statistical Sciences at the University of Padua.

# Installation
```
remotes::install_github("sacastiglioni/ClusteringSpartaco")
```
# Model Estimation
Let ```x``` be the matrix containing the expression of ```nrow(x)``` genes measured on ```ncol(x)``` spots, with genes in rows and spots in columns. The spatial coordinates of the spots are contained in the matrix ```coordinates```, while the spot labels are specified in the vector ```column.labels```. The function uses the parameters ```lambda.mu``` and ```lambda.tau``` to adjust the penalization of the model estimation, which are set to 0 by default. To perform the search for ```K``` gene clusters, execute the ```RCspartaco``` (Row-Clustering-spartaco) function with the following code:



```
library(ClusteringSpartaco)
RCspartaco(x = x, coordinates = coordinates, K = K, column.labels = column.labels)
```

# Convergence
By default, the estimation procedure is performed for a maximum of ```max.iter``` iterations, but is terminated earlier if a certain convergence criterion (```conv.criterion```) is reached. If ```conv.criterion = NULL```, there is no termination condition, and the procedure is executed for ```max.iter``` iterations. Otherwise, it is terminated when the increment of the classification log-likelihood is below a certain threshold ```conv.criterion$epsilon``` for ```conv.criterion$iterations``` consecutive times.

# Vignette
We show how to use the functions of the R package `ClusteringSpartaco` for fitting the clustering version of the SpaRTaCo model ([*Sottosanti and Risso, 2023*](https://arxiv.org/pdf/2110.04872.pdf)). The data and results can be retrieved at the links below:
- for the ```PathologistAnnotation.csv``` https://drive.google.com/file/d/1xunCdaatsB_Pzb9uIz99ar4l4ja8hFYE/view?usp=sharing
- for the ```spe_raw.RData``` data https://drive.google.com/file/d/1HAhzrN2bkN3Z4BtHZqU78Wck6GGKiBzh/view?usp=sharing or https://www.10xgenomics.com/resources/datasets/human-prostate-cancer-adenocarcinoma-with-invasive-carcinoma-ffpe-1-standard-1-3-0
- for the ```sign.symbol.RData``` https://drive.google.com/file/d/1Z0q2r8ftJQSlR8gb3LuAsC2XWdXCaFBH/view?usp=sharing
- for the model ```results``` https://drive.google.com/file/d/1VRGBoFdkCH3XOBbOM_KXft3hng_IYs5S/view?usp=sharing
