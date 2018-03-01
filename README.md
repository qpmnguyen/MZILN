# ZILN
## A multivariate zero-inflated logistic model for microbiome relative abundance data

### Things needed to be developed   
* Functionality to clean the data (remove all zero rows and all zero columns)  
* Functionality for input data 
	* Includes which covariates to be put in the model
	* Includes the number of taxa  
	* Includes the number of samples  
	* Data formatting is correct      
* Functionality of the code  
	* Proper cross-validation.  
* Functionality of output  
	* Includes the lambda curve  
	* Includes the taxa names divided by the number of possible covariates   
	* Includes the option of seeking for possible covariates and

### Stretch goals  
* CRAN support.  
* Implement Rcpp in order to improve performance.  
