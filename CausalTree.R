# LIBRARIES
library(tidyverse)
library(lubridate)
library(grf) # Generalized Random Forest for CausalTree estimation of heterogenous treatent effect
library(imputeTS)
library(gghighlight)
library(DiagrammeR)
library(ggthemes)
library(knitr)

# DATA IMPORT
WDIData <- read.csv("WDIData 1960-2019.csv", stringsAsFactors = F)
WDISeries <- read.csv("WDISeries.csv", stringsAsFactors = F)
brent_month <- read.csv("brent-month_csv 1987-2020.csv", stringsAsFactors = F)
WDICountry <- read.csv("WDICountry.csv", stringsAsFactors = F)

inflation_variables <- c("Inflation, GDP deflator (annual %)", "GDP deflator (base year varies by country)", "Inflation, consumer prices (annual %)")
inflation_Indicator.Code <- c("NY.GDP.DEFL.KD.ZG", "NY.GDP.DEFL.ZS", "FP.CPI.TOTL.ZG")

G20_countries <- c("ARG", "AUS", "BRA", "CAN", "CHN", "DEU", "FRA", "IND", "IDN", "ITA", "JPN", "MEX", "RUS", "SAU", "ZAF", "KOR", "TUR", "GBR" ,"USA")


# DATA WRANGLING
WDIData$Country.Name <- NULL
WDIData$Indicator.Name <- NULL
WDIData_1 <- gather(WDIData, key = "Year", value = "Value", -Country.Code, -Indicator.Code)
WDIData_2 <- spread(WDIData_1, key = "Indicator.Code", value = "Value")
WDIData_2$Year <- str_replace(WDIData_2$Year, "X", "") %>% as.numeric()

# removing inflation related variables' columns
WDIData_2A <- WDIData_2 %>% select(-NY.GDP.DEFL.KD.ZG, -NY.GDP.DEFL.ZS, -FP.CPI.TOTL.ZG)

# removing missing columns where more than 1/2 of observations are missing values
WDIData_col_NA <- map_int(WDIData_2A, ~sum(is.na(.x)))
nrow_WDIData_2A <- nrow(WDIData_2A)
WDIData_col_NA_which <- ifelse(WDIData_col_NA <= nrow_WDIData_2A*0.50, T, F) %>% which() # nrow_WDIData_2/2 because only 50% i.e. 1/2 missing values are allowed
WDIData_3 <- WDIData_2A[, WDIData_col_NA_which]

# removing rows where year is missing and truncating it to 1986-2016
WDIData_4 <- WDIData_3[which(!is.na(WDIData_3$Year)),]

# imputing missing values column wise using linear interpolation
WDIData_5 <- WDIData_4 %>% arrange(Country.Code, Year)
for(i in 1:ncol(WDIData_5)){
  WDIData_5[, i] <- na_interpolation(WDIData_5[, i], option = "linear")
}



# left joining outcome variable inflation
inflation_consumer <- WDIData_2 %>% select(Country.Code, Year, Inflation = FP.CPI.TOTL.ZG)
inflation_consumer <- inflation_consumer[which(!is.na(inflation_consumer$Year)),]
inflation_consumer <- inflation_consumer %>% arrange(Country.Code, Year)
inflation_consumer$Inflation <- na_interpolation(inflation_consumer$Inflation, option = "linear")
WDIData_6 <- left_join(WDIData_5, inflation_consumer, by = c("Country.Code", "Year"))

# binary treatment variable for oil price shock
# if increase in price is more than 50 % from 12 month lag then a price shock, only positive shocks included
brent_month$Date <- brent_month$Date %>% as.Date()
brent_month$Year <- year(brent_month$Date)
brent_month$price_change <- (brent_month$Price - lag(brent_month$Price, n = 12))/brent_month$Price
brent_month <- brent_month %>% filter(price_change > 0)
brent_month$Shock <- ifelse(brent_month$price_change >= 0.5, 1, 0)
brent_month$price_change <- NULL

treatment <- brent_month %>% 
  na.omit() %>% 
  group_by(Year) %>% 
  summarise(treatment = sum(Shock)) %>% 
  mutate(treatment = ifelse(treatment > 0, 1, 0))

# we will ad one year in treatment so that feature is one year back and is not affect by treatment
treatment$Year <- treatment$Year + 1

# joining treatment to our final dataset
WDIData_6 <- left_join(WDIData_6, treatment, by = "Year")
final_dataset <- WDIData_6 %>% na.omit()

# proportion of treatment variables
sum(final_dataset$treatment)/nrow(final_dataset) # 0.07142857

# arranging by country and year
final_dataset <- final_dataset %>% arrange(Country.Code, Year)

# DATA ANALYSIS
# estimating heterogenous treatment effect of oil price shock on inflation

# dividing data into train and test
set.seed(1996)
cases <- sample(seq_len(nrow(final_dataset)), round(nrow(final_dataset)*0.6))
train <- final_dataset[cases,]
test <- final_dataset[-cases,]

# estimating the causal forest
CF <- causal_forest(
  X = train[, -c(1, 2, ncol(final_dataset)-1, ncol(final_dataset))],
  Y = train$Inflation,
  W = train$treatment,
  num.trees = 5000,
  seed = 1996
)

preds <- predict(
  object = CF, 
  newdata = test[, -c(1, 2, ncol(final_dataset)-1, ncol(final_dataset))], 
  estimate.variance = TRUE
) %>% abs()

test <- cbind(test, preds)

variable_importance <- CF %>% 
  variable_importance() %>% 
  as.data.frame() %>% 
  mutate(variable = colnames(CF$X.orig)) %>% 
  arrange(desc(V1))

colnames(WDISeries)[1] <- colnames(variable_importance)[2]
variable_importance <- left_join(variable_importance, WDISeries, by = "variable")

# saving variable importance

write.csv(variable_importance, file = "variable_importanceA.csv")

# Estimating individual effect for G20 countries
G20_dataset <- final_dataset %>% filter(Country.Code %in% G20_countries, Year == 2019) 

G20_preds <- predict(
  object = CF, 
  newdata = G20_dataset[, -c(1, 2, ncol(G20_dataset)-1, ncol(G20_dataset))], 
  estimate.variance = TRUE) %>%
  abs()

G20_dataset <- cbind(G20_dataset, G20_preds)

# ATE India

IND_pred <- predict(
  object = CF, 
  newdata = final_dataset[, -c(1, 2, ncol(G20_dataset)-1, ncol(G20_dataset))] %>% filter(Country.Code == "IND"), 
  estimate.variance = TRUE) %>%
  abs()

# Average treatment effect
average_treatment_effect <- average_treatment_effect(CF) %>% abs()
paste("95% CI for the ATE:", round ( average_treatment_effect[1] , 3) ,"+/ -", round ( qnorm (0.975) * average_treatment_effect[2] , 3))


# PLOTS

# plotting relation between top four important varibles and the predictions

p1 <- ggplot(test, aes(x = log(test[,variable_importance$variable[1]]), y = predictions)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1) +
  theme_light() +
  xlab(variable_importance$Indicator.Name[1]) +
  ylab("Predictions")
p2 <- ggplot(test, aes(x = log(test[,variable_importance$variable[2]]), y = predictions)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1) +
  theme_light() +
  xlab(variable_importance$Indicator.Name[2]) +
  ylab("Predictions")
p3 <- ggplot(test, aes(x = log(test[,variable_importance$variable[3]]), y = predictions)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1) +
  theme_light() +
  xlab(variable_importance$Indicator.Name[3]) +
  ylab("Predictions")
p4 <- ggplot(test, aes(x = log(test[,variable_importance$variable[4]]), y = predictions)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1) +
  theme_light() +
  xlab(variable_importance$Indicator.Name[4]) +
  ylab("Predictions")

cowplot::plot_grid(p1, p2, p3, p4)

# plotting the heterogeneous treatment effect

plot_htes <- function(cf_preds, ci = FALSE, z = 1.96) {
  out <- ggplot(
    mapping = aes(
      x = rank(cf_preds$predictions), 
      y = cf_preds$predictions
    )
  ) +
    geom_point() +
    labs(x = "Rank", y = "Estimated Treatment Effect") +
    theme_light()
  
  if (ci == 1) {
    out <- out +
      geom_errorbar(
        mapping = aes(
          ymin = cf_preds$predictions + z * sqrt(cf_preds$variance.estimates),
          ymax = cf_preds$predictions - z * sqrt(cf_preds$variance.estimates)
        )
      )
  }
  
  return(out)
}
hte_plot <- plot_htes(preds[sample(1:nrow(preds), size = nrow(preds)/8),])
hte_plot_CI <- plot_htes(preds[sample(1:nrow(preds), size = nrow(preds)/8),], ci = T)

# G20 countries bar plot
G20_dataset %>%
  ggplot(aes(x = Country.Code, y = predictions)) +
  geom_col() +
  gghighlight(Country.Code == "IND") +
  ggtitle("Estimated Treatment Effect in 2019 for G20 Countries") +
  ylab("Predictions") +
  theme_economist()


# plotting one of the tree in forest
CF %>% get_tree(222) %>% plot()

# accuracy of predictions
rmse <- sqrt(sum((test$Inflation - test$predictions)^2)/nrow(test))

# FOR KNITTING RMD DOCUMENTS IN CURRENT SESSION
rmarkdown::render("causal_tree.Rmd", output_format = "word_document")


# ROUGH
A <- unique(WDIData$Indicator.Name) %>% length()
B <- str_replace("X1555", "X", "")
is.na(inflation_consumer$Inflation) %>% sum() # no NAs
WDIData_NA_omit <- na.omit(WDIData_2) # zero obs
any(is.na(WDIData_5)) # no NAs
sum(brent_month$Shock, na.rm = T) # 45 shocks
sum(treatment$treatment) # 12 shocks in 33 years
final_dataset$treatment %>% sum() # 2362 shock in 6648 observations i.e. approx 35%
temp <- final_dataset[, -c(1, 2, ncol(final_dataset)-1, ncol(final_dataset))]
plot(x = final_dataset$treatment, y = final_dataset$Inflation)
plot(x = final_dataset$Year, final_dataset$Inflation)
write.csv(variable_importance, file = "variable_importance.csv")
sum(is.na(WDIData))/(nrow(WDIData)*ncol(WDIData)) # 65% missing values
sum(is.na(WDIData_4))/(nrow(WDIData_4)*ncol(WDIData_4)) # 30 % missing values
CPIAData <- read.csv("CPIAData 2005-2018.csv")
WGIData <- read.csv("WGIData 1996-2018.csv")
sum(is.na(CPIAData))/(nrow(CPIAData)*ncol(CPIAData)) # 11 % missing values
sum(is.na(WGIData))/(nrow(WGIData)*ncol(WGIData)) # 6 % missing values
corruption_perceptions_index <- read.csv("corruption-perceptions-index data_csv.csv", na.strings = "-")
sum(is.na(corruption_perceptions_index))/(nrow(corruption_perceptions_index)*ncol(corruption_perceptions_index)) # 32 % missing values
max(final_dataset$Year) #2019
min(final_dataset$Year) #1989
any(is.na(WDIData_6)) # No
max(WDIData_6$Year) #2019
max(inflation_consumer$Year, na.rm = T) # 2019
# average treatment effect
mean(test$predictions) 
# -11.04526 abs  1-0
# -11.57539 no abs 1-0
# average treatment effect
mean(test$predictions) # -11.04526  / with abs 11.11803
# proposrtion of treatment in total observations
sum(final_dataset$treatment)/nrow(final_dataset) # 0.3548387
# Plotting relationship between top 4 variables and prediction

p1 <- ggplot(test, aes(x = test[,variable_importance$variable[1]], y = predictions)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1) +
  theme_light() +
  xlab(variable_importance$Indicator.Name[1]) +
  ylab("Predictions")

p2 <- ggplot(test, aes(x = NE.IMP.GNFS.CN, y = predictions)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1) +
  theme_light()

p3 <- ggplot(test, aes(x = NE.EXP.GNFS.CN, y = predictions)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1) +
  theme_light()

p4 <- ggplot(test, aes(x = NY.GDP.MKTP.KD.ZG, y = predictions)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1) +
  theme_light()

cowplot::plot_grid(p1, p2, p3, p4)

# inflation varibles
inf_var <- c("PA.NUS.PPP", "PA.NUS.PRVT.PP", "NY.GDP.DEFL.KD.ZG.AD", "NY.GDP.DEFL.ZS.AD", "PA.NUS.PPPC.RF")
# removing inflation variable
inf_var_which <- which(!colnames(final_dataset) %in% inf_var)
final_dataset <- final_dataset[, inf_var_which]
