#root -l -b -q Likelihood.C | tee testLikelihood.log 
root -l -b -q '../Scripts/Likelihood.C("Data")' | tee log_Likelihood_Data.log 
