install.packages(setdiff(c("dtw"), rownames(installed.packages())))
# Warning: this step pattern was found empirically. It is NOT actually normalizable because the weights
# are path-dependent.  See https://dynamictimewarping.github.io/faq/

asymmetricP1.1 <- dtw:::stepPattern(c(
  1, 1 , 2 , -1 ,
  1, 0 , 1 , 1 ,
  1, 0 , 0 , 1 ,
  2, 1 , 1 , -1 ,
  2, 0 , 0 ,  1 ,
  3, 2 , 1 , -1 ,
  3, 1 , 0 ,  2 ,
  3, 0 , 0 ,  2
),"N")
