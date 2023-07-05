# NCD risk reduction in Gaza
# Sanjay Basu, et al., 2023

library(haven)
library(tidyverse)
library(collapse)
library(globorisk)
library(labelled)
library(mice)
library(knitr)
library(kableExtra)
library(tableone)
library(ggplot2)

df0 <- read_dta("4. Dataset labelled clean.dta")

df1 = df0 %>%
  select(unique_id,DEM06,DEM04,DM02_2,DM03,BPS02,BPS03,BPS01,HBP03,TOBAC02,NARG02,BMI,WC01,PA01_time,LSTL04,LSTL05,CRD01,TOBAC01,NARG01) %>%
  remove_val_labels() 

df1$smk = df1$TOBAC02==1 | df1$TOBAC02==2 | df1$NARG02==1 | df1$NARG02==2 

df1 = df1 %>%
  select(-TOBAC02, -NARG02)

lab = read_dta("5. UNRWA_data_imported_noID.dta") %>%
  pivot_wider(names_from = TestName, values_from = TestResult) %>% 
  select(unique_id, "CHOLESTEROL - NCD", "CREATININE-NCD" ,   "TRIGLYCERIDES" ,    "LDL CHOLESTEROL" , "HDL CHOLESTEROL",   "HB FOR OTHERS"    ) %>%
  mutate_at(c("CHOLESTEROL - NCD", "CREATININE-NCD" ,   "TRIGLYCERIDES" ,    "LDL CHOLESTEROL" , "HDL CHOLESTEROL",   "HB FOR OTHERS"), as.numeric) %>%
  group_by(unique_id) %>% 
  summarise_each(~(flast(., na.rm = TRUE))) %>%
  remove_val_labels() %>%
  rename(tchol = "CHOLESTEROL - NCD", cr = "CREATININE-NCD" ,   tri = "TRIGLYCERIDES" ,  ldl=  "LDL CHOLESTEROL" , hdl ="HDL CHOLESTEROL",  hba1c= "HB FOR OTHERS" )

df2 = df1 %>%
  left_join(lab,by="unique_id") %>%
  select(-unique_id)

# Impute missing values using mice
imputed_data <- mice(df2, m = 1, maxit = 10, seed = 123)

# Generate the complete dataset
df <- complete(imputed_data)


df$age = df$DEM06
df$sex = df$DEM04==2
df$diabetes = df$DM02_2==1 | df$DM03!='NA' | (df$hba1c>=6.5)
df$systolic_bp <- NA

# Use only BPS01 if BPS02 and BPS03 are missing
df$systolic_bp[is.na(df$BPS02) & is.na(df$BPS03)] <- df$BPS01[is.na(df$BPS02) & is.na(df$BPS03)]

# Use average of BPS01 and BPS02 if BPS03 is missing
df$systolic_bp[is.na(df$BPS03) & !is.na(df$BPS02)] <- 
  (df$BPS01[is.na(df$BPS03) & !is.na(df$BPS02)] + df$BPS02[is.na(df$BPS03) & !is.na(df$BPS02)])/2

# Use average of all three if none are missing
df$systolic_bp[!is.na(df$BPS01) & !is.na(df$BPS02) & !is.na(df$BPS03)] <- 
  (df$BPS01[!is.na(df$BPS01) & !is.na(df$BPS02) & !is.na(df$BPS03)] + 
     df$BPS02[!is.na(df$BPS01) & !is.na(df$BPS02) & !is.na(df$BPS03)] +
     df$BPS03[!is.na(df$BPS01) & !is.na(df$BPS02) & !is.na(df$BPS03)])/3

df$hypertension_treatment_yes=(df$HBP03==1)

# convert from mg/dl to mmol/L
df$TCHOL =  df$tchol*0.02586
df$HDL =  df$hdl*0.02586
df$LDL = df$ldl*0.02586
  
# This section does the following:
#   - Calculates the Globorisk score for cardiovascular disease risk using the specified equations 
# - Calculates the FINDRISC score for diabetes risk, modified for Middle Eastern ethnicity 
# - Calculates a relative risk score for asthma/COPD using self-reported history and smoking exposure
# - Calculates the Gail model breast cancer risk score 
# - Calculates the CRC-PRO colorectal cancer risk score
# The scores are calculated for each individual in the df data frame using their provided characteristics.

# Globorisk score 
df$globage = df$age 
df$globage[df$age>80] = 80

globorisk_lbn = globorisk(
  sex = df$sex,
  age = df$globage,
  sbp = df$systolic_bp,
  tc = df$TCHOL,
  dm = df$diabetes,
  smk = df$smk,
  iso = rep("LBN",dim(df)[1]),
  year = rep(2020,dim(df)[1]),
  version = "lab",
  type = "risk"
)

globorisk_syr = globorisk(
  sex = df$sex,
  age = df$globage,
  sbp = df$systolic_bp,
  tc = df$TCHOL,
  dm = df$diabetes,
  smk = df$smk,
  iso = rep("SYR",dim(df)[1]),
  year = rep(2020,dim(df)[1]),
  version = "lab",
  type = "risk"
)

globorisk_jor = globorisk(
  sex = df$sex,
  age = df$globage,
  sbp = df$systolic_bp,
  tc = df$TCHOL,
  dm = df$diabetes,
  smk = df$smk,
  iso = rep("JOR",dim(df)[1]),
  year = rep(2020,dim(df)[1]),
  version = "lab",
  type = "risk"
)


globorisk <- rowMeans(cbind(globorisk_jor,
                   globorisk_lbn,
                   globorisk_syr))
  

# FINDRISC score 
# Function to calculate FINDRISC score
findrisc <- function(age, bmi, waist_circumference, physical_activity, fruit_vegetable_intake,
                     hypertension_treatment, high_blood_glucose, family_history) {
  score <- 0
  
  # Age
  if (age <= 54) {
    score <- score + 2
  } else if (age >= 55 & age <= 64) {
    score <- score + 3
  } else if (age > 64) {
    score <- score + 4
  }
  
  # BMI
  if (bmi >= 25 & bmi <= 30) {
    score <- score + 1
  } else if (bmi > 30) {
    score <- score + 3
  }
  
  # Waist circumference
  if (waist_circumference >= 94 & waist_circumference <= 102) {
    score <- score + 1
  } else if (waist_circumference > 102) {
    score <- score + 3
  }
  
  # Physical activity
  if (physical_activity == 0) {
    score <- score + 2
  }
  
  # Fruit and vegetable intake
  if (fruit_vegetable_intake == 0) {
    score <- score + 1
  }
  
  # Hypertension treatment
  if (hypertension_treatment == 1) {
    score <- score + 2
  }
  
  # High blood glucose
  if (high_blood_glucose == 1) {
    score <- score + 5
  }
  
  # Family history of diabetes
  if (family_history == 1) {
    score <- score + 5
  } else if (family_history == 2) {
    score <- score + 3
  }
  
  return(score)
}

age <- df$age
bmi <- df$BMI
waist_circumference <- df$WC01
physical_activity <- (df$PA01_time>30) # 1 for at least 30 min of physical activity daily, 0 for less
fruit_vegetable_intake <- (df$LSTL04==1 | df$LSTL05==1) # 1 for daily fruit/vegetable consumption, 0 for less
hypertension_treatment <- (df$hypertension_treatment_yes==1) # 1 for currently on blood pressure medication, 0 for not
high_blood_glucose <- (df$DM02==1) # 1 for high blood glucose, 0 for not
family_history <- rep(1,dim(df)[1]) # 0 for no family history, 1 for first-degree relative, 2 for second-degree relative

findrisc_score = rep(0,dim(df)[1])
findrisc_score_pre = rep(0,dim(df)[1])

# Define the function outside of the loop
risk_category <- function(findrisc_score) {
  if (findrisc_score >= 0 & findrisc_score <= 6) {
    return(0.005) #"Low risk (10-year risk < 1%)"
  } else if (findrisc_score >= 7 & findrisc_score <= 11) {
    return(0.025) # "Slightly elevated risk (10-year risk 1-5%)"
  } else if (findrisc_score >= 12 & findrisc_score <= 14) {
    return(0.10) #"Moderate risk (10-year risk 5-15%)"
  } else if (findrisc_score >= 15 & findrisc_score <= 20) {
    return(0.225) #"High risk (10-year risk 15-30%)"
  } else {
    return(0.3) # "Very high risk (10-year risk > 30%)"
  }
}

for (i in 1:length(age)){
  
  findrisc_score_pre[i] <- findrisc(age[i], bmi[i], waist_circumference[i], physical_activity[i], fruit_vegetable_intake[i],
                                    hypertension_treatment[i], high_blood_glucose[i], family_history[i])
  
  # Call the function inside the loop
  findrisc_score[i] <- risk_category(findrisc_score_pre[i])
  
}




df$copd_asthma_dx = (df$CRD01==1)
df$ever_smoker = (df$TOBAC01==1 | df$NARG01==1)

  
# Asthma/COPD risk 
#   This equation now calculates COPD/asthma risk as the sum of:
#   (Increased risk due to COPD/asthma diagnosis * risk due to smoking * risk from higher BMI * risk from lower BMI) 
# It uses the specified relative risks for:
# - Smoking (df$ever_smoker==1): Relative risk of 3.51 for ever smokers
# - Higher BMI (df$bmi - 25): Relative risk increases by 1.21 BMI over 25
# - Lower BMI (22 - df$bmi / 1.5): Relative risk of 1.5 for BMI below 22
# https://academic.oup.com/ije/article/47/6/1865/5113268
# https://bmcpulmmed.biomedcentral.com/articles/10.1186/1471-2466-11-36
# https://erj.ersjournals.com/content/35/6/1235
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7810869/

df$copd_asthma_risk_by_age = (df$sex==0)*(0.00002722538182*df$age-0.00054017418)+
  (df$sex==T)*(0.00001799354545*df$age+0.0002806929545)
copd_asthma_risk = max(1,3.51*(df$ever_smoker==1))*
  max(1,1.21*(df$BMI>25))*
  max(1,1.5*(df$BMI<22))*
  df$copd_asthma_risk_by_age




# Gail model
risk_factors <- c('age', 'age_at_menarche', 'first_degree_relatives', 
                  'number_of_biopsies', 'age_at_first_live_birth')

df$age_at_menarche = 13
df$first_degree_relatives = 0
df$number_of_biopsies = 0
df$age_at_first_live_birth = 24

# Calculate relative risks 
rr_age <- ifelse(df$age <= 50, 1, df$age * 0.025 - 0.5)
rr_menarche <- ifelse(df$age_at_menarche <= 12, 1.05, 1)  
rr_family <- ifelse(df$first_degree_relatives >= 1, 1.7, 1)
rr_biopsies <- ifelse(df$number_of_biopsies == 1, 1.3,  
                      ifelse(df$number_of_biopsies == 2, 1.6,  
                             ifelse(df$number_of_biopsies >= 3, 2.45, 1)))
rr_birth <- ifelse(df$age_at_first_live_birth <= 24,  
                   1.05 * (35 - df$age_at_first_live_birth) * 0.01, 1)

# Calculate risk 
base_risk <- (57.14/100000)*(df$sex)+(0.48/100000*(1-df$sex))  # age standardized incidence in Palestine
df$age_risk <- base_risk * rr_age 
df$menarche_risk <- df$age_risk * rr_menarche  
df$family_risk <- df$menarche_risk * rr_family
df$biopsy_risk <- df$family_risk * rr_biopsies
df$birth_risk <- df$biopsy_risk * rr_birth
gail_model <- round(df$birth_risk * (1 - exp(-0.25 * 10)), 4) 

# Function to calculate CRC-Pro model risk
crc_pro <- function(age, sex, family_history, smoking, bmi) {
  # Coefficients for the logistic regression model
  intercept <- -8.09
  age_coef <- 0.07
  sex_coef <- 0.31
  family_history_coef <- 0.69
  smoking_coef <- 0.18
  bmi_coef <- 0.03
  
  # Calculate the linear predictor
  linear_predictor <- intercept + age_coef * age + sex_coef * sex + family_history_coef * family_history + smoking_coef * smoking + bmi_coef * bmi
  
  # Calculate the risk using logistic function
  risk <- exp(linear_predictor) / (1 + exp(linear_predictor))
  
  return(risk)
}

df$family_history <- 0 # 0 for no family history, 1 for family history
crcrisk <- crc_pro(df$age, (1-df$sex), family_history, df$smk, df$BMI)


# This section:
# - Generates individual risk estimates for each condition using the risk calculators (Globorisk, FINDRISC, etc.)
# - Calibrates those individual risks to match Palestine incidence rates from GBD
# - Sums the calibrated individual risks for each condition to determine the total 10-year incidence
# - Rounds the incidence to whole numbers to estimate cases for each condition over 10 years

# Set population size and time horizon
pop_size <- 10000  
time_horizon <- 10
cv_inc <- (2679.58+756.49)/100000   # CV inc in Palestine (IHD + Stroke) 
diabetes_inc <- 1395.5/100000 #  diabetes inc in Palestine
copd_asthma_inc <- (389.25+316.94)/100000   #  COPD/asthma inc in Palestine         
breast_inc <-  98.8/100000     #  breast ca inc in Palestine
colorectal_ca_inc <- 112.04/100000 #  CRC inc in Palestine
cv_mortality_prop <- (839.33+475.41)/(21755.91+5738.55)   # CV mortality in Palestine  
diabetes_mortality_prop <- 294.79/1395.5 #  diabetes mortality in Palestine
copd_asthma_mortality_prop <- (66.63+20.33)/(389.25+316.94)   # COPD/asthma mortality in Palestine         
breast_ca_mortality_prop <-  52/98.8     # breast ca mortality in Palestine
colorectal_ca_mortality_prop <- 82.69/112.04 # CRC mortality in Palestine



generate_summaries <- function(cv_inc,diabetes_inc,copd_asthma_inc,breast_inc,colorectal_ca_inc,
                               cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop, colorectal_ca_mortality_prop) {

# Globorisk CV risk  
cv_risk <- globorisk   
# Scale risk to match annual  incidence in Palestine
cv_scale_adj <- cv_inc*time_horizon / mean(cv_risk,na.rm=T)  
cv_risk_caled_adj <- cv_risk * cv_scale_adj
cv_risk_caled_adj[cv_risk_caled_adj>1] = 1

# Diabetes risk 
diabetes_risk <- findrisc_score   
# Scale risk to match annual  incidence in Palestine
diabetes_scale_adj <- diabetes_inc*time_horizon / mean(diabetes_risk,na.rm=T)  
diabetes_risk_caled_adj <- diabetes_risk * diabetes_scale_adj
diabetes_risk_caled_adj[diabetes_risk_caled_adj>1]=1

# COPD/asthma risk  
copd_asthma_risk <- copd_asthma_risk   
# Scale risk to match annual  incidence in Palestine
copd_asthma_scale_adj <- copd_asthma_inc*time_horizon / mean(copd_asthma_risk,na.rm=T)  
copd_asthma_risk_caled_adj <- copd_asthma_risk * copd_asthma_scale_adj
copd_asthma_risk_caled_adj[copd_asthma_risk_caled_adj>1]=1

# Breast cancer risk 
breast_ca_risk <- gail_model
# Scale risk to match annual  incidence in Palestine
breast_ca_scale_adj <- breast_inc*time_horizon / mean(breast_ca_risk,na.rm=T)  
breast_ca_risk_caled_adj <- breast_ca_risk * breast_ca_scale_adj
breast_ca_risk_caled_adj[breast_ca_risk_caled_adj>1]=1

# CRC risk   
colorectal_ca_risk <- crcrisk   
# Scale risk to match annual  incidence in Palestine
colorectal_ca_scale_adj <- colorectal_ca_inc*time_horizon / mean(colorectal_ca_risk,na.rm=T)  
colorectal_ca_risk_caled_adj <- colorectal_ca_risk * colorectal_ca_scale_adj
colorectal_ca_risk_caled_adj[colorectal_ca_risk_caled_adj>1]=1


# Discount rate
discount <- 0.03

iters = 1000
results_list <- list()
for (i in 1:iters) {
  
# Simulate cases based on binomial probability 
set.seed(i)
cv_cases <- rbinom(pop_size, 1, cv_risk_caled_adj)  
set.seed(i)
diabetes_cases <- rbinom(pop_size, 1, diabetes_risk_caled_adj)  
set.seed(i)
copd_asthma_cases <- rbinom(pop_size, 1, copd_asthma_risk_caled_adj)  
set.seed(i)
breast_ca_cases <- rbinom(pop_size, 1, breast_ca_risk_caled_adj)  
set.seed(i)
colorectal_ca_cases <- rbinom(pop_size, 1, colorectal_ca_risk_caled_adj)

# Calculate incidence 
cv_incidence <- sum(cv_cases) / pop_size  
diabetes_incidence <- sum(diabetes_cases) / pop_size  
copd_asthma_incidence <- sum(copd_asthma_cases) / pop_size  
breast_ca_incidence <- sum(breast_ca_cases) / pop_size  
colorectal_ca_incidence <- sum(colorectal_ca_cases) / pop_size 

# Cases over time horizon
cv_Cases <- cv_incidence * pop_size * time_horizon /10
diabetes_Cases <- diabetes_incidence * pop_size * time_horizon /10
copd_asthma_Cases <- copd_asthma_incidence * pop_size * time_horizon /10
breast_ca_Cases <- breast_ca_incidence * pop_size * time_horizon /10
colorectal_ca_Cases <- colorectal_ca_incidence * pop_size * time_horizon /10

# Simulate mortality based on binomial probability
set.seed(i)
cv_mortality <- rbinom(cv_Cases, 1, cv_mortality_prop)  
set.seed(i)
diabetes_mortality <- rbinom(diabetes_Cases, 1, diabetes_mortality_prop)
set.seed(i)
copd_asthma_mortality <- rbinom(copd_asthma_Cases, 1, copd_asthma_mortality_prop)  
set.seed(i)
breast_ca_mortality <- rbinom(breast_ca_Cases, 1, breast_ca_mortality_prop)      
set.seed(i)
colorectal_ca_mortality <- rbinom(colorectal_ca_Cases, 1, colorectal_ca_mortality_prop )



# This section: 
# - Calculates DALYs lost for each condition over 10 years based on disability weights from the literature 
# - Sums years of life lost and years of live lifed with disability and sums DALYs at a 3% anual discount rate

# Disability weights https://ghdx.healthdata.org/record/ihme-data/gbd-2019-disability-weights
cv_disability_weight <- 0.236295455  
diabetes_disability_weight <- 0.147875  
copd_disability_weight <- 0.1045  
breast_ca_disability_weight <- 0.241166667
colorectal_ca_disability_weight <- 0.260333333 

# Average age of disease onset  
cv_age_onset <- 51 #   https://pubmed.ncbi.nlm.nih.gov/15364185/
diabetes_age_onset <- 55 #   https://pubmed.ncbi.nlm.nih.gov/21245594/
copd_asthma_age_onset <- 43 # https://www.sciencedirect.com/science/article/pii/S0954611112700124#:~:text=This%20was%20most%20frequently%20chronic,was%2043.0%C2%B111.7%20years.  
breast_ca_age_onset <- 49 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7325818/#:~:text=4.2.&text=Najjar%20and%20Easson%20%5B15%5D%20found,than%2050%20years%20%5B16%5D.  
colorectal_ca_age_onset <- 60  # https://ascopubs.org/doi/full/10.1200/EDBK_390520#:~:text=The%20mean%20age%20at%20CRC,age%2060%20v%2071%20years).&text=In%20addition%2C%20the%20incidence%20of,and%20in%20the%20Middle%20East.

# Age distribution 
ages <- df$age

# Discount rate
discount <- 0.03

# Simulate age at disease onset 
cv_age_at_onset <- rep(cv_age_onset,cv_Cases)
diabetes_age_at_onset <- rep(diabetes_age_onset,diabetes_Cases)
copd_asthma_age_at_onset <- rep(copd_asthma_age_onset,copd_asthma_Cases)
breast_ca_age_at_onset <- rep(breast_ca_age_onset, breast_ca_Cases)  
colorectal_ca_age_at_onset <- rep(colorectal_ca_age_onset, colorectal_ca_Cases)


# Set life expectancy 
life_expectancy <- 74

# Simulate age at death based on disease duration and life expectancy  
# The durations are  sampled as:
# - CV: 0 to 8 years; prevalence/incidence = duration: (21756+5738)/(2679+756)
# - Diabetes: 0 to 24 years : 33526/1396
# - COPD/asthma: 0 to 12: (3587+5204)/(317+389)
# - Breast cancer: 0 to 8 years: 811/99
# - Colorectal cancer: 0 to 4: 438/112

set.seed(i)
cv_age_at_death <- pmin(life_expectancy, cv_age_at_onset + sample(0:8, length(cv_age_at_onset), replace=TRUE))  
set.seed(i)
diabetes_age_at_death <- pmin(life_expectancy, diabetes_age_at_onset + sample(0:24, length(diabetes_age_at_onset), replace=TRUE)) 
set.seed(i)
copd_asthma_age_at_death <- pmin(life_expectancy, copd_asthma_age_at_onset + sample(0:12, length(copd_asthma_age_at_onset), replace=TRUE))    
set.seed(i)
breast_ca_age_at_death <- pmin(life_expectancy, breast_ca_age_at_onset + sample(0:8, length(breast_ca_age_at_onset), replace=TRUE))         
set.seed(i)
colorectal_ca_age_at_death <- pmin(life_expectancy, colorectal_ca_age_at_onset + sample(0:4, length(colorectal_ca_age_at_onset), replace=TRUE))  

# Years of life lost (YLL)
cv_yll <- (life_expectancy - cv_age_at_death) * cv_mortality * (1 - discount)^(cv_age_at_death - cv_age_at_onset)  
diabetes_yll <- (life_expectancy - diabetes_age_at_death) * diabetes_mortality * (1 - discount)^(diabetes_age_at_death - diabetes_age_at_onset)
copd_asthma_yll <- (life_expectancy - copd_asthma_age_at_death) * copd_asthma_mortality * (1 - discount)^(copd_asthma_age_at_death - copd_asthma_age_at_onset)
breast_ca_yll <- (life_expectancy - breast_ca_age_at_death) * breast_ca_mortality * (1 - discount)^(breast_ca_age_at_death - breast_ca_age_at_onset)   
colorectal_ca_yll <- (life_expectancy - colorectal_ca_age_at_death) * colorectal_ca_mortality * (1 - discount)^(colorectal_ca_age_at_death - colorectal_ca_age_at_onset)

# Years lived with disability (YLD)
cv_yld <- cv_incidence * cv_disability_weight * (cv_age_at_death - cv_age_at_onset) * (1 - discount)^(cv_age_at_death - cv_age_at_onset)  
diabetes_yld <- diabetes_incidence * diabetes_disability_weight * (diabetes_age_at_death - diabetes_age_at_onset) * (1 - discount)^(diabetes_age_at_death - diabetes_age_at_onset)
copd_asthma_yld <- copd_asthma_incidence * copd_disability_weight * (copd_asthma_age_at_death - copd_asthma_age_at_onset) * (1 - discount)^(copd_asthma_age_at_death - copd_asthma_age_at_onset)  
breast_ca_yld <- breast_ca_incidence * breast_ca_disability_weight * (breast_ca_age_at_death - breast_ca_age_at_onset) * (1 - discount)^(breast_ca_age_at_death - breast_ca_age_at_onset)     
colorectal_ca_yld <- colorectal_ca_incidence * colorectal_ca_disability_weight * (colorectal_ca_age_at_death - colorectal_ca_age_at_onset) * (1 - discount)^(colorectal_ca_age_at_death - colorectal_ca_age_at_onset)



# DALYs
cv_dalys <- sum(cv_yll + cv_yld )
diabetes_dalys <- sum(diabetes_yll + diabetes_yld)
copd_asthma_dalys <- sum(copd_asthma_yll + copd_asthma_yld  )
breast_ca_dalys <- sum(breast_ca_yll + breast_ca_yld)
colorectal_ca_dalys <- sum(colorectal_ca_yll + colorectal_ca_yld)


# Create data frame with results
results <- data.frame(
  Condition = c("CVD", "Diabetes", "COPD/Asthma", "Breast Cancer", "Colorectal Cancer"),
  Incidence = c(cv_Cases, diabetes_Cases, copd_asthma_Cases, breast_ca_Cases, 
            colorectal_ca_Cases),
  DALYs = c(cv_dalys, diabetes_dalys, copd_asthma_dalys, breast_ca_dalys, 
            colorectal_ca_dalys),
  Mortality = c(sum(cv_mortality), sum(diabetes_mortality), sum(copd_asthma_mortality),
            sum(breast_ca_mortality), sum(colorectal_ca_mortality))
)


results_list[[i]] <- results
}

results_summary <- do.call(rbind, results_list)
}



results_summary = generate_summaries(cv_inc,diabetes_inc,copd_asthma_inc,breast_inc,colorectal_ca_inc,
                                     cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop, colorectal_ca_mortality_prop)


# proportion of population reached by intervention, divided by 2 for assuming 50% coverage on average over 10 years
prop =0.64/2

results_summary_1_4 = generate_summaries(cv_inc*(1-prop + prop* (mean(2.08/2.88,2.71/3.85))),diabetes_inc,copd_asthma_inc*(1-prop + prop* (1-0.049)),breast_inc,colorectal_ca_inc,
                                         cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop, colorectal_ca_mortality_prop)

results_summary_2_11 = generate_summaries(cv_inc*(1-prop + prop* (1-(0.521-0.475)+0.82*((0.521-0.475)))),diabetes_inc*(1-prop + prop* (1-(0.521-0.475)+0.86*((0.521-0.475)))),copd_asthma_inc,breast_inc*(1-prop + prop* (1-(0.521-0.475)+0.89*((0.521-0.475)))),colorectal_ca_inc*(1-prop + prop* (1-(0.521-0.475)+0.91*((0.521-0.475)))),
                                         cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop, colorectal_ca_mortality_prop)

results_summary_2_13 = generate_summaries(cv_inc*(1-prop + prop* (1-(0.07)+0.82*((0.07)))),diabetes_inc*(1-prop + prop* (1-(0.07)+0.86*((0.07)))),copd_asthma_inc,breast_inc*(1-prop + prop* (1-(0.07)+0.89*((0.07)))),colorectal_ca_inc*(1-prop + prop* (1-(0.07)+0.91*((0.07)))),
                                          cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop, colorectal_ca_mortality_prop)

results_summary_3_1 = generate_summaries(cv_inc*(1-prop + prop* (1-(0.042)+0.76*((0.042)))),diabetes_inc*(1-prop + prop* (1-(0.042)+0.69*((0.042)))),copd_asthma_inc*(1-prop + prop* (1-(0.042)+0.72*((0.042)))),breast_inc*(1-prop + prop* (1-(0.042)+0.81*((0.042)))),colorectal_ca_inc*(1-prop + prop* (1-(0.042)+0.76*((0.042)))),
                                          cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop, colorectal_ca_mortality_prop)

results_summary_3_2 = generate_summaries(cv_inc*(1-prop + prop* (1-(0.11)+0.76*((0.11)))),diabetes_inc*(1-prop + prop* (1-(0.11)+0.69*((0.11)))),copd_asthma_inc*(1-prop + prop* (1-(0.11)+0.72*((0.11)))),breast_inc*(1-prop + prop* (1-(0.11)+0.81*((0.11)))),colorectal_ca_inc*(1-prop + prop* (1-(0.11)+0.76*((0.11)))),
                                         cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop, colorectal_ca_mortality_prop)

results_summary_6_2 = generate_summaries(cv_inc,diabetes_inc,copd_asthma_inc,breast_inc,colorectal_ca_inc,
                                         cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop*(1-prop + prop* (1-0.28)), colorectal_ca_mortality_prop)

results_summary_6_3 = generate_summaries(cv_inc,diabetes_inc,copd_asthma_inc,breast_inc,colorectal_ca_inc,
                                         cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop, colorectal_ca_mortality_prop*(1-prop) + prop*(1-0.65))

results_summary_6_4 = generate_summaries(cv_inc,diabetes_inc,copd_asthma_inc,breast_inc,colorectal_ca_inc,
                                         cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop*(1-prop)+prop*(1-0.86), colorectal_ca_mortality_prop)

results_summary_7_3 = generate_summaries(cv_inc,diabetes_inc,copd_asthma_inc*(1-prop + prop* (1-0.15)),breast_inc,colorectal_ca_inc,
                                         cv_mortality_prop,diabetes_mortality_prop, copd_asthma_mortality_prop, breast_ca_mortality_prop, colorectal_ca_mortality_prop)




iters = 1000
discount = 0.03

# Extract DALYs for each set of results
dalys_results_summary <- sum(results_summary[,3])/iters
dalys_results_summary_1_4 <- sum(results_summary_1_4[,3])/iters

# Calculate difference in DALYs
daly_difference_1_4 <- dalys_results_summary - dalys_results_summary_1_4

# Define incremental costs 
incremental_costs_1_4 <- pop_size*rnorm(iters,0.20,0.20/2/1.96)*time_horizon * (1 - discount)^time_horizon
quantile_cost_diff_1_4 <- quantile(incremental_costs_1_4, probs = c(0.025, 0.975))  

# Calculate ICER  
ICER_1_4 <- mean(incremental_costs_1_4 / daly_difference_1_4)

# Define 95% quantiles for DALY difference  
quantile_daly_diff_1_4 <- quantile(daly_difference_1_4, probs = c(0.025, 0.975))  

# Calculate 95% quantile bounds for ICER  
icer_lower_1_4 <- quantile_cost_diff_1_4[2] / quantile_daly_diff_1_4[1]  
icer_upper_1_4 <- quantile_cost_diff_1_4[1] / quantile_daly_diff_1_4[2]  

# Print results  
cat("The ICER between results_summary and results_summary_1_4 is: ", ICER_1_4,   "\n")
cat("Lower 95% bound:", icer_lower_1_4, "; Upper 95% bound:", icer_upper_1_4, "\n")




# Extract DALYs for each set of results
dalys_results_summary <- sum(results_summary[,3])/iters
dalys_results_summary_2_11 <- sum(results_summary_2_11[,3])/iters

# Calculate difference in DALYs
daly_difference_2_11 <- dalys_results_summary - dalys_results_summary_2_11

# Define incremental costs 
incremental_costs_2_11 <- pop_size*rnorm(iters,4,4*7/57)*time_horizon * (1 - discount)^time_horizon
quantile_cost_diff_2_11 <- quantile(incremental_costs_2_11, probs = c(0.025, 0.975))  

# Calculate ICER  
ICER_2_11 <- mean(incremental_costs_2_11 / daly_difference_2_11)

# Define 95% quantiles for DALY difference  
quantile_daly_diff_2_11 <- quantile(daly_difference_2_11, probs = c(0.025, 0.975))  

# Calculate 95% quantile bounds for ICER  
icer_lower_2_11 <- quantile_cost_diff_2_11[2] / quantile_daly_diff_2_11[1]  
icer_upper_2_11 <- quantile_cost_diff_2_11[1] / quantile_daly_diff_2_11[2]  

# Print results  
cat("The ICER between results_summary and results_summary_2_11 is: ", ICER_2_11,   "\n")
cat("Lower 95% bound:", icer_lower_2_11, "; Upper 95% bound:", icer_upper_2_11, "\n")






# Extract DALYs for each set of results
dalys_results_summary <- sum(results_summary[,3])/iters
dalys_results_summary_2_13 <- sum(results_summary_2_13[,3])/iters

# Calculate difference in DALYs
daly_difference_2_13 <- dalys_results_summary - dalys_results_summary_2_13

# Define incremental costs 
incremental_costs_2_13 <- pop_size*rnorm(iters,0.30,0.30/2/1.96)*time_horizon * (1 - discount)^time_horizon
quantile_cost_diff_2_13 <- quantile(incremental_costs_2_13, probs = c(0.025, 0.975))  

# Calculate ICER  
ICER_2_13 <- mean(incremental_costs_2_13 / daly_difference_2_13)

# Define 95% quantiles for DALY difference  
quantile_daly_diff_2_13 <- quantile(daly_difference_2_13, probs = c(0.025, 0.975))  

# Calculate 95% quantile bounds for ICER  
icer_lower_2_13 <- quantile_cost_diff_2_13[2] / quantile_daly_diff_2_13[1]  
icer_upper_2_13 <- quantile_cost_diff_2_13[1] / quantile_daly_diff_2_13[2]  

# Print results  
cat("The ICER between results_summary and results_summary_2_13 is: ", ICER_2_13,   "\n")
cat("Lower 95% bound:", icer_lower_2_13, "; Upper 95% bound:", icer_upper_2_13, "\n")







# Extract DALYs for each set of results
dalys_results_summary <- sum(results_summary[,3])/iters
dalys_results_summary_3_1 <- sum(results_summary_3_1[,3])/iters

# Calculate difference in DALYs
daly_difference_3_1 <- dalys_results_summary - dalys_results_summary_3_1

# Define incremental costs 
incremental_costs_3_1 <- pop_size*rnorm(iters,18.9,18.9/2/1.96)*time_horizon * (1 - discount)^time_horizon
quantile_cost_diff_3_1 <- quantile(incremental_costs_3_1, probs = c(0.025, 0.975))  

# Calculate ICER  
ICER_3_1 <- mean(incremental_costs_3_1 / daly_difference_3_1)

# Define 95% quantiles for DALY difference  
quantile_daly_diff_3_1 <- quantile(daly_difference_3_1, probs = c(0.025, 0.975))  

# Calculate 95% quantile bounds for ICER  
icer_lower_3_1 <- quantile_cost_diff_3_1[2] / quantile_daly_diff_3_1[1]  
icer_upper_3_1 <- quantile_cost_diff_3_1[1] / quantile_daly_diff_3_1[2]  

# Print results  
cat("The ICER between results_summary and results_summary_3_1 is: ", ICER_3_1,   "\n")
cat("Lower 95% bound:", icer_lower_3_1, "; Upper 95% bound:", icer_upper_3_1, "\n")


\






# Extract DALYs for each set of results
dalys_results_summary <- sum(results_summary[,3])/iters
dalys_results_summary_3_2 <- sum(results_summary_3_2[,3])/iters

# Calculate difference in DALYs
daly_difference_3_2 <- dalys_results_summary - dalys_results_summary_3_2

# Define incremental costs 
incremental_costs_3_2 <- pop_size*rnorm(iters,110,110/2/1.96)*time_horizon * (1 - discount)^time_horizon
quantile_cost_diff_3_2 <- quantile(incremental_costs_3_2, probs = c(0.025, 0.975))  

# Calculate ICER  
ICER_3_2 <- mean(incremental_costs_3_2 / daly_difference_3_2)

# Define 95% quantiles for DALY difference  
quantile_daly_diff_3_2 <- quantile(daly_difference_3_2, probs = c(0.025, 0.975))  

# Calculate 95% quantile bounds for ICER  
icer_lower_3_2 <- quantile_cost_diff_3_2[2] / quantile_daly_diff_3_2[1]  
icer_upper_3_2 <- quantile_cost_diff_3_2[1] / quantile_daly_diff_3_2[2]  

# Print results  
cat("The ICER between results_summary and results_summary_3_2 is: ", ICER_3_2,   "\n")
cat("Lower 95% bound:", icer_lower_3_2, "; Upper 95% bound:", icer_upper_3_2, "\n")







# Extract DALYs for each set of results
dalys_results_summary <- sum(results_summary[,3])/iters
dalys_results_summary_6_2 <- sum(results_summary_6_2[,3])/iters

# Calculate difference in DALYs
daly_difference_6_2 <- dalys_results_summary - dalys_results_summary_6_2

# Define incremental costs 
incremental_costs_6_2 <- pop_size*rnorm(iters,19.95,19.95/2/1.96)*time_horizon * (1 - discount)^time_horizon
quantile_cost_diff_6_2 <- quantile(incremental_costs_6_2, probs = c(0.025, 0.975))  

# Calculate ICER  
ICER_6_2 <- mean(incremental_costs_6_2 / daly_difference_6_2)

# Define 95% quantiles for DALY difference  
quantile_daly_diff_6_2 <- quantile(daly_difference_6_2, probs = c(0.025, 0.975))  

# Calculate 95% quantile bounds for ICER  
icer_lower_6_2 <- quantile_cost_diff_6_2[2] / quantile_daly_diff_6_2[1]  
icer_upper_6_2 <- quantile_cost_diff_6_2[1] / quantile_daly_diff_6_2[2]  

# Print results  
cat("The ICER between results_summary and results_summary_6_2 is: ", ICER_6_2,   "\n")
cat("Lower 95% bound:", icer_lower_6_2, "; Upper 95% bound:", icer_upper_6_2, "\n")












# Extract DALYs for each set of results
dalys_results_summary <- sum(results_summary[,3])/iters
dalys_results_summary_6_3 <- sum(results_summary_6_3[,3])/iters

# Calculate difference in DALYs
daly_difference_6_3 <- dalys_results_summary - dalys_results_summary_6_3

# Define incremental costs 
incremental_costs_6_3 <- colorectal_ca_inc*pop_size*rnorm(iters,14736.09392,14736.09392/2/1.96)*time_horizon * (1 - discount)^time_horizon
quantile_cost_diff_6_3 <- quantile(incremental_costs_6_3, probs = c(0.025, 0.975))  

# Calculate ICER  
ICER_6_3 <- mean(incremental_costs_6_3 / daly_difference_6_3)

# Define 95% quantiles for DALY difference  
quantile_daly_diff_6_3 <- quantile(daly_difference_6_3, probs = c(0.025, 0.975))  

# Calculate 95% quantile bounds for ICER  
icer_lower_6_3 <- quantile_cost_diff_6_3[2] / quantile_daly_diff_6_3[1]  
icer_upper_6_3 <- quantile_cost_diff_6_3[1] / quantile_daly_diff_6_3[2]  

# Print results  
cat("The ICER between results_summary and results_summary_6_3 is: ", ICER_6_3,   "\n")
cat("Lower 95% bound:", icer_lower_6_3, "; Upper 95% bound:", icer_upper_6_3, "\n")











# Extract DALYs for each set of results
dalys_results_summary <- sum(results_summary[,3])/iters
dalys_results_summary_6_4 <- sum(results_summary_6_4[,3])/iters

# Calculate difference in DALYs
daly_difference_6_4 <- dalys_results_summary - dalys_results_summary_6_4

# Define incremental costs 
incremental_costs_6_4 <- breast_inc*pop_size*rnorm(iters,2214,2214/2/1.96)*time_horizon * (1 - discount)^time_horizon
quantile_cost_diff_6_4 <- quantile(incremental_costs_6_4, probs = c(0.025, 0.975))  

# Calculate ICER  
ICER_6_4 <- mean(incremental_costs_6_4 / daly_difference_6_4)

# Define 95% quantiles for DALY difference  
quantile_daly_diff_6_4 <- quantile(daly_difference_6_4, probs = c(0.025, 0.975))  

# Calculate 95% quantile bounds for ICER  
icer_lower_6_4 <- quantile_cost_diff_6_4[2] / quantile_daly_diff_6_4[1]  
icer_upper_6_4 <- quantile_cost_diff_6_4[1] / quantile_daly_diff_6_4[2]  

# Print results  
cat("The ICER between results_summary and results_summary_6_4 is: ", ICER_6_4,   "\n")
cat("Lower 95% bound:", icer_lower_6_4, "; Upper 95% bound:", icer_upper_6_4, "\n")









# Extract DALYs for each set of results
dalys_results_summary <- sum(results_summary[,3])/iters
dalys_results_summary_7_3 <- sum(results_summary_7_3[,3])/iters

# Calculate difference in DALYs
daly_difference_7_3 <- dalys_results_summary - dalys_results_summary_7_3

# Define incremental costs 
incremental_costs_7_3 <- copd_asthma_inc*pop_size*rnorm(iters,24,24/2/1.96)*time_horizon * (1 - discount)^time_horizon
quantile_cost_diff_7_3 <- quantile(incremental_costs_7_3, probs = c(0.025, 0.975))  

# Calculate ICER  
ICER_7_3 <- mean(incremental_costs_7_3 / daly_difference_7_3)

# Define 95% quantiles for DALY difference  
quantile_daly_diff_7_3 <- quantile(daly_difference_7_3, probs = c(0.025, 0.975))  

# Calculate 95% quantile bounds for ICER  
icer_lower_7_3 <- quantile_cost_diff_7_3[2] / quantile_daly_diff_7_3[1]  
icer_upper_7_3 <- quantile_cost_diff_7_3[1] / quantile_daly_diff_7_3[2]  

# Print results  
cat("The ICER between results_summary and results_summary_7_3 is: ", ICER_7_3,   "\n")
cat("Lower 95% bound:", icer_lower_7_3, "; Upper 95% bound:", icer_upper_7_3, "\n")








  
# Extract numeric columns 
cvd_incidence <- as.numeric(results_summary$Incidence[results_summary$Condition == "CVD"])
cvd_dalys <- as.numeric(results_summary$DALYs[results_summary$Condition == "CVD"]) 
cvd_mortality <- as.numeric(results_summary$Mortality[results_summary$Condition == "CVD"])

# Extract numeric columns 
cvd_incidence <- as.numeric(results_summary$Incidence[results_summary$Condition == "CVD"])
cvd_dalys <- as.numeric(results_summary$DALYs[results_summary$Condition == "CVD"])  
cvd_mortality <- as.numeric(results_summary$Mortality[results_summary$Condition == "CVD"])

diabetes_incidence <- as.numeric(results_summary$Incidence[results_summary$Condition == "Diabetes"])  
diabetes_dalys <- as.numeric(results_summary$DALYs[results_summary$Condition == "Diabetes"])  
diabetes_mortality <- as.numeric(results_summary$Mortality[results_summary$Condition == "Diabetes"])  

copd_asthma_incidence <- as.numeric(results_summary$Incidence[results_summary$Condition == "COPD/Asthma"])  
copd_asthma_dalys <- as.numeric(results_summary$DALYs[results_summary$Condition == "COPD/Asthma"])    
copd_asthma_mortality <- as.numeric(results_summary$Mortality[results_summary$Condition == "COPD/Asthma"])

breast_ca_incidence <- as.numeric(results_summary$Incidence[results_summary$Condition == "Breast Cancer"])  
breast_ca_dalys <- as.numeric(results_summary$DALYs[results_summary$Condition == "Breast Cancer"])    
breast_ca_mortality <- as.numeric(results_summary$Mortality[results_summary$Condition == "Breast Cancer"])   

colorectal_ca_incidence <- as.numeric(results_summary$Incidence[results_summary$Condition == "Colorectal Cancer"])  
colorectal_ca_dalys <- as.numeric(results_summary$DALYs[results_summary$Condition == "Colorectal Cancer"])    
colorectal_ca_mortality <- as.numeric(results_summary$Mortality[results_summary$Condition == "Colorectal Cancer"])   



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


df$dmcorr = (df0$DM05==1 | df0$DM04==1 | df0$DM03!='NA') 
df$dmcorr[is.na(df$dmcorr)] = 'FALSE'

df$dmmd = (df0$DM05==1 | df0$DM04==1)
df$dmmd[is.na(df$dmmd)] = 'FALSE'


df$htndx = df0$HBP02==1
df$htndx[is.na(df$htndx)] = 'FALSE'

df$htnmdcorr = (df0$HBP03==1)
df$htnmdcorr[is.na(df$htnmdcorr)] = 'FALSE'

df$cig = df0$TOBAC01==1
df$cig[is.na(df$cig)] = 'FALSE'

df$pipe = df0$NARG01==1
df$pipe[is.na(df$pipe)] = 'FALSE'

df$dyslip = df0$CHOL02==1
df$dyslip[is.na(df$dyslip)] = 'FALSE'


df$statin = df0$CVD04==1
df$statin[is.na(df$statin)] = 'FALSE'


# Specify the variables to summarize
vars <- c("age", "dmcorr",  "dmmd", "htndx", "htnmdcorr", "systolic_bp", "dyslip", "statin", "copd_asthma_dx","cig", "pipe", "BMI")

# Create a TableOne object stratified by sex
mytable <- CreateTableOne(vars = vars, strata = "sex", data = df)

# Print the table
print(mytable,quote = T,smd = T,
      catDigits = 1,
      contDigits = 1)


dfsex = df0 %>%
  select(unique_id,DEM04)

labdf = lab %>%
  left_join(dfsex)
# Specify the variables to summarize
vars <- c("tchol", "tri", "ldl", "hdl", "hba1c", "cr")

# Create a TableOne object stratified by sex
mytable2 <- CreateTableOne(vars = vars, strata = "DEM04", data = labdf)

# Print the table
print(mytable2,quote = T,smd = T,
      catDigits = 1,
      contDigits = 1)


interventions <- read.table(text = "
                           Intervention IncrementalCostmean CostSD IncrementalDALYmean DALYSD
                           1_4 14719.28 3665.102 433 2
                           2_11 294799.4 35951.72 19 4.6
                           2_13 22187.1 5681.36 30 6.1
                           3_1 1386611  344407 37 7.7
                           3_2 8186501  1944942 94 0.5
                           6_2 1471433 373880.3 87 19.9
                           6_3 1223437  315883.3 160 7.1 
                           6_4 158439.7  40794.48 217 19.9 
                           7_3 12588.36  2978.863 90 3.1", 
                            header = TRUE)

# Sample from normal distributions 1,000 times  
boot_costs <- lapply(1:nrow(interventions), function(i) { 
  rnorm(1000, interventions$IncrementalCostmean[i], 
        interventions$CostSD[i]) 
})
boot_effects <- lapply(1:nrow(interventions), function(i) {  
  rnorm(1000, interventions$IncrementalDALYmean[i], 
        interventions$DALYSD[i])
})

# Combine results into data frame        
ce_results <- data.frame(
  cost = unlist(boot_costs),
  effects = unlist(boot_effects)
)

ce_results$intervention <- 
  rep(interventions$Intervention, each=1000)

# Plot 
ggplot(ce_results, aes(x = effects, y = cost, 
                       color = intervention, shape = intervention)) +
  geom_point() +
  scale_color_manual(values =  
                       c("#1f78b4","#ff7f00", "#a6cee3", "#b2df8a",  
                         "#33a02c","#fb9a99", "#e31a1c", "#fdbf6f", "#cab2d6")) +  
  scale_shape_manual(values = 1:9)+ 
  xlab("Incremental DALYs") + 
  ylab("Incremental costs ($2023)")


ggplot(ce_results, aes(x = effects, y = cost, 
                       color = intervention, shape = intervention)) +
  geom_point() +
  scale_color_manual(values =  
                       c("#1f78b4","#ff7f00", "#a6cee3", "#b2df8a",  
                         "#33a02c","#fb9a99", "#e31a1c", "#fdbf6f", "#cab2d6")) +  
  scale_shape_manual(values = 1:9)+ 
  xlab("Incremental DALYs") + 
  ylab("Incremental costs ($2023)")+
  ylim(0, 3e6)  



