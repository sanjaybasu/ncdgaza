
# Extract numeric columns 
cvd_incidence <- as.numeric(results_summary_1_4$Incidence[results_summary_1_4$Condition == "CVD"])
cvd_dalys <- as.numeric(results_summary_1_4$DALYs[results_summary_1_4$Condition == "CVD"]) 
cvd_mortality <- as.numeric(results_summary_1_4$Mortality[results_summary_1_4$Condition == "CVD"])

# Extract numeric columns 
cvd_incidence <- as.numeric(results_summary_1_4$Incidence[results_summary_1_4$Condition == "CVD"])
cvd_dalys <- as.numeric(results_summary_1_4$DALYs[results_summary_1_4$Condition == "CVD"])  
cvd_mortality <- as.numeric(results_summary_1_4$Mortality[results_summary_1_4$Condition == "CVD"])

diabetes_incidence <- as.numeric(results_summary_1_4$Incidence[results_summary_1_4$Condition == "Diabetes"])  
diabetes_dalys <- as.numeric(results_summary_1_4$DALYs[results_summary_1_4$Condition == "Diabetes"])  
diabetes_mortality <- as.numeric(results_summary_1_4$Mortality[results_summary_1_4$Condition == "Diabetes"])  

copd_asthma_incidence <- as.numeric(results_summary_1_4$Incidence[results_summary_1_4$Condition == "COPD/Asthma"])  
copd_asthma_dalys <- as.numeric(results_summary_1_4$DALYs[results_summary_1_4$Condition == "COPD/Asthma"])    
copd_asthma_mortality <- as.numeric(results_summary_1_4$Mortality[results_summary_1_4$Condition == "COPD/Asthma"])

breast_ca_incidence <- as.numeric(results_summary_1_4$Incidence[results_summary_1_4$Condition == "Breast Cancer"])  
breast_ca_dalys <- as.numeric(results_summary_1_4$DALYs[results_summary_1_4$Condition == "Breast Cancer"])    
breast_ca_mortality <- as.numeric(results_summary_1_4$Mortality[results_summary_1_4$Condition == "Breast Cancer"])   

colorectal_ca_incidence <- as.numeric(results_summary_1_4$Incidence[results_summary_1_4$Condition == "Colorectal Cancer"])  
colorectal_ca_dalys <- as.numeric(results_summary_1_4$DALYs[results_summary_1_4$Condition == "Colorectal Cancer"])    
colorectal_ca_mortality <- as.numeric(results_summary_1_4$Mortality[results_summary_1_4$Condition == "Colorectal Cancer"])   



mean_cvd <- c(Condition = "CVD", Incidence = mean(cvd_incidence), DALYs = mean(cvd_dalys), Mortality = mean(cvd_mortality))
quantile_cvd <- cbind(
  Condition = "CVD",
  Incidence = quantile(cvd_incidence, probs = c(0.025, 0.975)), 
  DALYs = quantile(cvd_dalys, probs = c(0.025, 0.975)), 
  Mortality = quantile(cvd_mortality, probs = c(0.025, 0.975))  
)

mean_diabetes <- c(Condition = "Diabetes", Incidence = mean(diabetes_incidence), DALYs = mean(diabetes_dalys), Mortality = mean(diabetes_mortality))
quantile_diabetes <- cbind(
  Condition = "Diabetes",
  Incidence = quantile(diabetes_incidence, probs = c(0.025, 0.975)), 
  DALYs = quantile(diabetes_dalys, probs = c(0.025, 0.975)), 
  Mortality = quantile(diabetes_mortality, probs = c(0.025, 0.975))  
) 

mean_copd_asthma <- c(Condition = "COPD/Asthma", Incidence = mean(copd_asthma_incidence), DALYs = mean(copd_asthma_dalys), Mortality = mean(copd_asthma_mortality))
quantile_copd_asthma <- cbind(
  Condition = "COPD/Asthma",
  Incidence = quantile(copd_asthma_incidence, probs = c(0.025, 0.975)),  
  DALYs = quantile(copd_asthma_dalys, probs = c(0.025, 0.975)),  
  Mortality = quantile(copd_asthma_mortality, probs = c(0.025, 0.975))    
)   

mean_breast_ca <- c(Condition = "Breast Cancer", Incidence = mean(breast_ca_incidence), DALYs = mean(breast_ca_dalys), Mortality = mean(breast_ca_mortality))  
quantile_breast_ca <- cbind(
  Condition = "Breast Cancer",
  Incidence = quantile(breast_ca_incidence, probs = c(0.025, 0.975)),  
  DALYs = quantile(breast_ca_dalys, probs = c(0.025, 0.975)),  
  Mortality = quantile(breast_ca_mortality, probs = c(0.025, 0.975))   
)

mean_colorectal_ca <- c(Condition = "Colorectal Cancer", Incidence = mean(colorectal_ca_incidence), DALYs = mean(colorectal_ca_dalys), Mortality = mean(colorectal_ca_mortality))  
quantile_colorectal_ca <- cbind(
  Condition = "Colorectal Cancer",
  Incidence = quantile(colorectal_ca_incidence, probs = c(0.025, 0.975)), 
  DALYs = quantile(colorectal_ca_dalys, probs = c(0.025, 0.975)), 
  Mortality = quantile(colorectal_ca_mortality, probs = c(0.025, 0.975))   
)


mean_results <- rbind(mean_cvd, mean_diabetes, mean_copd_asthma, mean_breast_ca, mean_colorectal_ca)
quantile_results <- rbind(quantile_cvd, quantile_diabetes, quantile_copd_asthma, quantile_breast_ca, quantile_colorectal_ca)



