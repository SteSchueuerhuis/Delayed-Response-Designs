# A Two-Stage Group-Sequential Design for Delayed Treatment Responses with the Possibility of Trial Restart

Authors: Stephen Schüürhuis, Frank Konietschke, Cornelia Ursula Kunz

The code to reproduce the plots and results of the article is demonstrated in the Quarto document reproduce_results.qmd / reproduce_results.html in the folder application. 

The file init.R needs to be run to import packages and read in the code files.

# Information on R

The code was produced with the following versions of R and packages:

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252   
[3] LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.1252    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] kableExtra_1.3.4    tibble_3.1.7        ggpubr_0.4.0        tictoc_1.1         
 [5] ggplot2_3.4.1       nleqslv_3.3.2       crayon_1.5.1        numDeriv_2016.8-1.1
 [9] greekLetters_0.0.7  mvtnorm_1.1-3       rpact_3.4.0         MASS_7.3-57        

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.8.3      svglite_2.1.1     tidyr_1.2.0       assertthat_0.2.1  digest_0.6.31    
 [6] utf8_1.2.2        R6_2.5.1          backports_1.4.1   evaluate_0.22     httr_1.4.3       
[11] highr_0.9         pillar_1.9.0      rlang_1.1.3       rstudioapi_0.13   jquerylib_0.1.4  
[16] car_3.0-13        R.utils_2.11.0    R.oo_1.24.0       styler_1.10.2     rmarkdown_2.14   
[21] textshaping_0.3.6 webshot_0.5.4     stringr_1.4.0     munsell_0.5.0     broom_1.0.3      
[26] compiler_4.1.2    xfun_0.31         pkgconfig_2.0.3   systemfonts_1.0.4 htmltools_0.5.5  
[31] tidyselect_1.2.0  gridExtra_2.3     viridisLite_0.4.0 fansi_1.0.3       dplyr_1.1.0      
[36] withr_2.5.0       R.methodsS3_1.8.1 grid_4.1.2        jsonlite_1.8.4    gtable_0.3.0     
[41] lifecycle_1.0.3   pacman_0.5.1      magrittr_2.0.3    scales_1.2.1      cli_3.3.0        
[46] stringi_1.7.6     carData_3.0-5     farver_2.1.0      ggsignif_0.6.3    bslib_0.3.1      
[51] xml2_1.3.3        ellipsis_0.3.2    ragg_1.2.5        generics_0.1.2    vctrs_0.6.1      
[56] cowplot_1.1.1     tools_4.1.2       R.cache_0.15.0    glue_1.6.2        purrr_0.3.4      
[61] abind_1.4-5       fastmap_1.1.0     colorspace_2.0-3  rstatix_0.7.0     rvest_1.0.2      
[66] knitr_1.39        sass_0.4.5  

