library(testthat)
library(mixedsubjects)

test_check("mixedsubjects")

library(devtools)                  
load_all()                                                                               
test() 