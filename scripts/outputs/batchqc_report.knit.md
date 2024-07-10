---
title: "BatchQC Report"
date: "2024-06-28"
output: 
  html_vignette:
    toc: true
    toc_depth: 2
    template: batchQC.html
    self_contained: no
    lib_dir: libs
---


Summary
=======
## Confounding
### Number of samples in each Batch and Condition

-----------------------------------------------------
        &nbsp;           Batch 1   Batch 2   Batch 3 
----------------------- --------- --------- ---------
 **Condition control**     72        14        12    
-----------------------------------------------------

### Measures of confounding between Batch and Condition

----------------------------------------------------------------------
            &nbsp;                Standardized Pearson     Cramer's V 
                                 Correlation Coefficient              
------------------------------- ------------------------- ------------
  **Confounding Coefficients               NA                  NA     
 (0=no confounding, 1=complete                                        
        confounding)**                                                
----------------------------------------------------------------------

## Variation Analysis
### Variation explained by Batch and Condition
![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


----------------------------------------------------------
   &nbsp;      Full (Condition+Batch)   Condition   Batch 
------------- ------------------------ ----------- -------
  **Min.**             0.001                0       0.001 

 **1st Qu.**           12.25                0       12.25 

 **Median**            35.22                0       35.22 

  **Mean**             38.59                0       38.59 

 **3rd Qu.**           64.14                0       64.14 

  **Max.**              91.1                0       91.1  
----------------------------------------------------------

## P-value Analysis
### Distribution of Batch and Condition Effect p-values Across Genes

--------------------------------------------------------------------------------------------
         &nbsp;           Min.   1st Qu.    Median     Mean     3rd Qu.     Max.    Ps<0.05 
------------------------ ------ --------- ---------- --------- ---------- -------- ---------
   **Batch P-values**      0        0      1.17e-09   0.06358   0.002014   0.9997   0.8592  

 **Condition P-values**    1        1         1          1         1         1         0    
--------------------------------------------------------------------------------------------

![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


Differential Expression
=======================
## Expression Plot
Boxplots for all values for each of the samples and are colored by batch membership.

![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

## LIMMA



Median Correlations
===================
This plot helps identify outlying samples.
![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


Heatmaps
========
## Heatmap
This is a heatmap of the given data matrix showing the batch effects and variations with different conditions.
![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

## Sample Correlations
This is a heatmap of the correlation between samples.
![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-16-1.png)<!-- -->


Circular Dendrogram
===================
This is a Circular Dendrogram of the given data matrix colored by batch to show the batch effects.
![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-18-1.png)<!-- -->


PCA: Principal Component Analysis
=================================
## PCA
This is a plot of the top two principal components colored by batch to show the batch effects.
![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

## Explained Variation

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  &nbsp;    Proportion of Variance (%)   Cumulative Proportion of   Percent Variation Explained by   Percent Variation Explained by   Condition Significance   Percent Variation Explained by   Batch Significance (p-value) 
                                               Variance (%)           Either Condition or Batch                Condition                    (p-value)                      Batch                                             
---------- ---------------------------- -------------------------- -------------------------------- -------------------------------- ------------------------ -------------------------------- ------------------------------
 **PC1**              37.49                       37.49                          88.4                              0                            1                           88.4                             0               

 **PC2**              7.752                       45.24                          33.3                              0                            1                           33.3                             0               

 **PC3**               6.32                       51.56                          24.8                              0                            1                           24.8                             0               

 **PC4**              4.357                       55.92                          0.1                               0                            1                           0.1                            0.9332            

 **PC5**              3.839                       59.76                          17.5                              0                            1                           17.5                          0.00011            

 **PC6**              3.359                       63.11                          5.7                               0                            1                           5.7                           0.06164            

 **PC7**              2.677                       65.79                          2.4                               0                            1                           2.4                            0.3204            

 **PC8**               2.41                        68.2                          3.6                               0                            1                           3.6                            0.1713            

 **PC9**              2.105                       70.31                          6.6                               0                            1                           6.6                           0.03931            

 **PC10**             1.822                       72.13                          1.2                               0                            1                           1.2                            0.5705            

 **PC11**             1.462                       73.59                          0.2                               0                            1                           0.2                            0.9169            

 **PC12**             1.312                        74.9                           1                                0                            1                            1                             0.6218            

 **PC13**             1.243                       76.15                           2                                0                            1                            2                             0.3906            

 **PC14**             1.198                       77.34                           1                                0                            1                            1                             0.6231            

 **PC15**             1.128                       78.47                          1.8                               0                            1                           1.8                            0.4254            

 **PC16**             1.097                       79.57                          0.1                               0                            1                           0.1                            0.939             

 **PC17**             0.9689                      80.54                          0.4                               0                            1                           0.4                            0.8265            

 **PC18**             0.9022                      81.44                          0.6                               0                            1                           0.6                            0.747             

 **PC19**             0.8939                      82.33                          1.3                               0                            1                           1.3                            0.5485            

 **PC20**             0.8664                       83.2                          0.5                               0                            1                           0.5                            0.8062            

 **PC21**             0.8113                      84.01                          0.3                               0                            1                           0.3                            0.8659            

 **PC22**             0.7637                      84.78                          0.2                               0                            1                           0.2                            0.9024            

 **PC23**             0.7254                       85.5                          0.1                               0                            1                           0.1                            0.9415            

 **PC24**             0.6757                      86.18                          0.2                               0                            1                           0.2                            0.9239            

 **PC25**             0.6518                      86.83                          0.8                               0                            1                           0.8                            0.6672            

 **PC26**             0.6252                      87.45                          0.2                               0                            1                           0.2                            0.8965            

 **PC27**             0.5659                      88.02                          0.2                               0                            1                           0.2                            0.9212            

 **PC28**             0.5528                      88.57                          0.1                               0                            1                           0.1                            0.9462            

 **PC29**             0.5487                      89.12                          0.1                               0                            1                           0.1                            0.975             

 **PC30**             0.5268                      89.65                          0.1                               0                            1                           0.1                            0.9737            

 **PC31**             0.4987                      90.15                          0.2                               0                            1                           0.2                            0.9098            

 **PC32**             0.4832                      90.63                          0.8                               0                            1                           0.8                            0.667             

 **PC33**             0.4711                       91.1                          0.4                               0                            1                           0.4                            0.826             

 **PC34**             0.4574                      91.56                          0.1                               0                            1                           0.1                            0.9317            

 **PC35**             0.4444                        92                           0.1                               0                            1                           0.1                            0.9731            

 **PC36**             0.4268                      92.43                          0.2                               0                            1                           0.2                            0.9116            

 **PC37**             0.4036                      92.83                           0                                0                            1                            0                             0.9997            

 **PC38**             0.3756                      93.21                           0                                0                            1                            0                             0.9991            

 **PC39**             0.3386                      93.55                          0.2                               0                            1                           0.2                            0.9171            

 **PC40**             0.3333                      93.88                           0                                0                            1                            0                             0.9801            

 **PC41**             0.3171                       94.2                           0                                0                            1                            0                             0.9981            

 **PC42**             0.2983                       94.5                          0.1                               0                            1                           0.1                            0.9724            

 **PC43**             0.2901                      94.79                          0.1                               0                            1                           0.1                            0.9549            

 **PC44**             0.2814                      95.07                           0                                0                            1                            0                             0.9887            

 **PC45**             0.2703                      95.34                          0.1                               0                            1                           0.1                            0.9628            

 **PC46**             0.2528                      95.59                          0.2                               0                            1                           0.2                            0.9088            

 **PC47**             0.2403                      95.83                          0.1                               0                            1                           0.1                            0.9592            

 **PC48**             0.2255                      96.06                          0.2                               0                            1                           0.2                             0.91             

 **PC49**             0.2222                      96.28                          0.1                               0                            1                           0.1                            0.9671            

 **PC50**             0.2122                      96.49                          0.3                               0                            1                           0.3                            0.8708            

 **PC51**             0.1931                      96.68                           0                                0                            1                            0                             0.9819            

 **PC52**             0.1842                      96.87                          0.3                               0                            1                           0.3                            0.877             

 **PC53**             0.1816                      97.05                           0                                0                            1                            0                             0.9904            

 **PC54**             0.1618                      97.21                          0.1                               0                            1                           0.1                            0.9548            

 **PC55**             0.1569                      97.37                           0                                0                            1                            0                             0.9818            

 **PC56**             0.1548                      97.52                           0                                0                            1                            0                             0.9926            

 **PC57**             0.152                       97.68                           0                                0                            1                            0                             0.9979            

 **PC58**             0.1442                      97.82                          0.1                               0                            1                           0.1                            0.9705            

 **PC59**             0.129                       97.95                          0.1                               0                            1                           0.1                            0.9745            

 **PC60**             0.1272                      98.08                           0                                0                            1                            0                             0.9873            

 **PC61**             0.122                        98.2                          0.1                               0                            1                           0.1                            0.9502            

 **PC62**             0.1177                      98.32                          0.2                               0                            1                           0.2                            0.8954            

 **PC63**             0.1015                      98.42                           0                                0                            1                            0                             0.9928            

 **PC64**            0.09906                      98.52                           0                                0                            1                            0                             0.9975            

 **PC65**             0.0949                      98.61                           0                                0                            1                            0                             0.9962            

 **PC66**             0.0921                       98.7                          0.1                               0                            1                           0.1                            0.9555            

 **PC67**            0.08657                      98.79                          0.2                               0                            1                           0.2                            0.9308            

 **PC68**            0.08169                      98.87                           0                                0                            1                            0                             0.9842            

 **PC69**            0.07945                      98.95                          0.1                               0                            1                           0.1                            0.9457            

 **PC70**            0.07505                      99.03                           0                                0                            1                            0                             0.9945            

 **PC71**            0.07453                       99.1                          0.2                               0                            1                           0.2                            0.9302            

 **PC72**            0.06996                      99.17                           0                                0                            1                            0                             0.9984            

 **PC73**            0.06435                      99.23                           0                                0                            1                            0                             0.993             

 **PC74**            0.05953                      99.29                           0                                0                            1                            0                             0.9914            

 **PC75**            0.05783                      99.35                           0                                0                            1                            0                             0.9988            

 **PC76**            0.05478                      99.41                          0.2                               0                            1                           0.2                            0.8937            

 **PC77**            0.04936                      99.46                           0                                0                            1                            0                             0.9921            

 **PC78**            0.04794                       99.5                           0                                0                            1                            0                             0.9859            

 **PC79**             0.0458                      99.55                           0                                0                            1                            0                             0.9981            

 **PC80**            0.04344                      99.59                           0                                0                            1                            0                             0.9995            

 **PC81**            0.04176                      99.63                           0                                0                            1                            0                             0.9806            

 **PC82**            0.03753                      99.67                           0                                0                            1                            0                             0.9936            

 **PC83**            0.03579                      99.71                           0                                0                            1                            0                             0.9795            

 **PC84**            0.03454                      99.74                           0                                0                            1                            0                             0.983             

 **PC85**            0.03098                      99.77                           0                                0                            1                            0                             0.9986            

 **PC86**            0.02966                       99.8                           0                                0                            1                            0                             0.9918            

 **PC87**            0.02691                      99.83                           0                                0                            1                            0                             0.997             

 **PC88**            0.02476                      99.86                           0                                0                            1                            0                             0.9974            

 **PC89**            0.02271                      99.88                           0                                0                            1                            0                             0.9992            

 **PC90**            0.02142                       99.9                           0                                0                            1                            0                             0.9907            

 **PC91**            0.02032                      99.92                           0                                0                            1                            0                             0.9994            

 **PC92**            0.01765                      99.94                           0                                0                            1                            0                             0.9992            

 **PC93**            0.01669                      99.95                           0                                0                            1                            0                             0.9988            

 **PC94**            0.01479                      99.97                           0                                0                            1                            0                             0.9937            

 **PC95**            0.01257                      99.98                           0                                0                            1                            0                             0.9992            

 **PC96**            0.01059                      99.99                           0                                0                            1                            0                             0.9993            

 **PC97**            0.008152                      100                            0                                0                            1                            0                             0.998             

 **PC98**           7.393e-30                      100                           35.9                              0                            1                           35.9                             0               
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Shape
=====
This is a heatmap plot showing the variation of gene expression mean, variance, skewness and kurtosis between samples grouped by batch to see the batch effects variation
![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```
## Note: Sample-wise p-value is calculated for the variation across samples on the measure across genes. Gene-wise p-value is calculated for the variation of each gene between batches on the measure across each batch. If the data is quantum normalized, then the Sample-wise measure across genes is same for all samples and Gene-wise p-value is a good measure.
```


Combat Plots
============
This is a plot showing whether parametric or non-parameteric prior is appropriate for this data. It also shows the Kolmogorov-Smirnov test comparing the parametric and non-parameteric prior distribution.

```
## Found 3 batches
## Adjusting for 0 covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
```

![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-25-1.png)<!-- -->![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-25-2.png)<!-- -->![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-25-3.png)<!-- -->![](C:/Users/sarah.dandou/Documents/MELANOMODELE/melano-db/melano/scripts/batchqc_report_files/figure-html/unnamed-chunk-25-4.png)<!-- -->

```
## Batch mean distribution across genes: Normal vs Empirical distribution
## Two-sided Kolmogorov-Smirnov test
## Selected Batch: 1
## Statistic D = 0.06166
## p-value = 0.08715
## 
## 
## Batch Variance distribution across genes: Inverse Gamma vs Empirical distribution
## Two-sided Kolmogorov-Smirnov test
## Selected Batch: 1
## Statistic D = 0.3811
## p-value = 0Note: The non-parametric version of ComBat takes much longer time to run and we recommend it only when the shape of the non-parametric curve widely differs such as a bimodal or highly skewed distribution. Otherwise, the difference in batch adjustment is very negligible and parametric version is recommended even if p-value of KS test above is significant.
```


SVA
===
## Summary

```
## Number of Surrogate Variables found in the given data: 1
```
