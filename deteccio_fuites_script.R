#######################
# DETECCIÓ DE FUITES  #
#######################

# Carregar llibreries
library(caret)
library(car)
library(corrplot)
library(data.table)
library(datawizard)
library(dplyr)
library(forcats)
library(ggplot2)
library(GGally)
library(glmnet)
library(graphics)
library(ggraph)
library(igraph)
library(isotone)    
library(lubridate)
library(MASS)
library(pdp)
library(pROC)
library(PRROC)
library(psych)
library(ranger)
library(readr)
library(showtext)
library(skimr)
library(smotefamily)
library(tibble)
library(tidyr)
library(tidygraph)
library(tidyverse)
library(xgboost)
library(MLmetrics)
library(ROSE)


# CARREGAR DADES

df <- read_csv("dataset_fugas_gas_realista_v4.csv")

# Taula de classes de 'fuga'
table(df$fuga)
# Proporció de cada classe
prop.table(table(df$fuga))
# Estadístiques descriptives detallades
describe(df)
# Estadístiques per classe de 'fuga'
df %>%
  group_by(fuga) %>%
  summarise(across(where(is.numeric), list(media = mean, sd = sd),
                   .names = "{.col}_{.fn}"))


# Seleccionar variables numèriques
numeric_vars <- sapply(df, is.numeric)
hist_data    <- df[, numeric_vars]

# Taula descriptiva amb noms en català
descripcion <- purrr::map_dfr(
  .x = names(hist_data),
  .f = function(varname) {
    vec <- hist_data[[varname]]
    tibble(
      variable        = varname,
      n               = sum(!is.na(vec)),
      mitjana         = mean(vec, na.rm = TRUE),
      màxim           = max(vec, na.rm = TRUE),
      mínim           = min(vec, na.rm = TRUE),
      mediana         = median(vec, na.rm = TRUE),
      sd              = sd(vec, na.rm = TRUE),
      skew            = skewness(vec, na.rm = TRUE)
    )
  }
)

knitr::kable(
  descripcion,
  digits = c(0, 1, 5, 5, 5, 5, 3, 2),
  caption = "Estadístiques descriptives de les variables numèriques: n, mitjana, màxim, mínim, mediana, sd i asimetria."
)

# Matriu de correlacions
cor_matrix <- cor(df)
corrplot(cor_matrix, method = "shade", tl.cex = 0.5)

# Correlació amb la variable 'fuga'
cor_fuga <- cor_matrix[, "fuga"]
print(cor_fuga)
# Correlacions ordenades descendentment
cor_fuga_sorted <- sort(cor_fuga, decreasing = TRUE)
print(cor_fuga_sorted)

# Top 5 variables més correlacionades positivament amb 'fuga'
top5_pos <- head(cor_fuga_sorted[names(cor_fuga_sorted) != "fuga"], 5)
print(top5_pos)
# Top 5 variables més correlacionades negativament

top5_neg <- head(sort(cor_fuga, decreasing = FALSE)[names(cor_fuga) != "fuga"], 5)
print(top5_neg)
# Filtrar correlacions amb |r| > 0.1
cor_fuga_filtrat <- cor_fuga_sorted[abs(cor_fuga_sorted) > 0.1]
print(cor_fuga_filtrat)


# Boxplot: delta_flux_nodo vs fuga
ggplot(df, aes(x = factor(fuga), y = delta_flujo_nodo, fill = factor(fuga))) +
  geom_boxplot() +
  labs(title = "Delta de flux segons classe de fuga", x = "Fuga", y = "Delta flux")

# Boxplot: zscore_presio_local vs fuga
ggplot(df, aes(x = factor(fuga), y = zscore_presion_local, fill = factor(fuga))) +
  geom_boxplot() +
  labs(title = "Zscore de pressió local segons classe de fuga", x = "Fuga", y = "Zscore pressió local")

# Histograma: pressió real per classe
ggplot(df, aes(x = presion_real, fill = factor(fuga))) +
  geom_histogram(position = "identity", bins = 50, alpha = 0.6) +
  labs(title = "Distribució de pressió real segons fuga", fill = "Fuga")

# Histograma: pressió sensor per classe
ggplot(df, aes(x = presion_sensor, fill = factor(fuga))) +
  geom_histogram(position = "identity", bins = 50, alpha = 0.6) +
  labs(title = "Distribució de pressió sensor segons fuga", fill = "Fuga")

# Densitat: centralitat per classe
ggplot(df, aes(x = centralitat, fill = factor(fuga))) +
  geom_density(alpha = 0.6) +
  labs(title = "Distribució de centralitat segons fuga", fill = "Fuga")


#Test de diferències de mitjanes

t.test(df$delta_flujo_nodo ~ df$fuga)
t.test(df$delta_presion_media_vecinos ~ df$fuga)
t.test(df$flujo_entrada_tuberia ~ df$fuga)
t.test(df$zscore_presion_local ~ df$fuga)
t.test(df$centralidad ~ df$fuga)
t.test(df$flujo_salida_tuberia ~ df$fuga)
t.test(df$grado_conectividad ~ df$fuga)


# VIF i selecció de variables

# Predictors seleccionats
predictors <- c(
  "presion_real","presion_sensor","flujo_entrada_tuberia","flujo_salida_tuberia",
  "delta_flujo_nodo","delta_presion_media_vecinos","centralidad","distancia_fuente",
  "grado_conectividad","es_consumo_discreto","hora_dia","log_presion_real",
  "std_presion_vecinos","zscore_presion_local"
)

# Dataframe reduït
df2 <- df %>% select(fuga, all_of(predictors))

# Càlcul del VIF
mod_vif <- lm(fuga ~ ., data = df2)
vif_vals <- vif(mod_vif)
print(vif_vals)
high_vif <- vif_vals[vif_vals > 5]
print(high_vif)

# Anàlisi de Components Principals (PCA)
X_scaled <- scale(df2 %>% dplyr::select(-fuga))
pca <- prcomp(X_scaled, center = TRUE, scale. = FALSE)
summary(pca)
plot(pca, type = "l")


# Llegir el dataset i transformar variables
df <- read_csv("dataset_fugas_gas_realista_v4.csv") %>%
  mutate(
    fuga = factor(fuga, levels = c(0,1)),
    hora_cat = case_when(
      hora_dia %in% 0:5   ~ "madrugada",
      hora_dia %in% 6:11  ~ "matí",
      hora_dia %in% 12:17 ~ "tarda",
      hora_dia %in% 18:23 ~ "nit"
    ) %>% factor(levels = c("madrugada","matí","tarda","nit"))
  )

# Variables predictives
vars_pred <- c(
  "delta_flujo_nodo",
  "flujo_entrada_tuberia",
  "zscore_presion_local",
  "delta_presion_media_vecinos",
  "centralidad",
  "flujo_salida_tuberia",
  "grado_conectividad",
  "distancia_fuente",
  "es_consumo_discreto",
  "std_presion_vecinos",
  "hora_cat"
)

# Divisió Train/Test
set.seed(42)
idx_train <- createDataPartition(df$fuga, p = 0.7, list = FALSE)
df_train <- df[idx_train, ]
df_test  <- df[-idx_train, ]

# Codificació One-Hot de les categories a train i test
ohe <- dummyVars(~ ., data = df_train[vars_pred], fullRank = TRUE)
X_train <- predict(ohe, df_train[vars_pred]) %>% as.data.frame()
X_test  <- predict(ohe, df_test[vars_pred])  %>% as.data.frame()

# Adjuntar la variable resposta
df_train_ohe <- X_train %>% mutate(fuga = df_train$fuga)
y_train_num  <- as.numeric(as.character(df_train_ohe$fuga))
y_test_num   <- as.numeric(as.character(df_test$fuga))

# Reequilibri amb SMOTE sobre el conjunt d'entrenament
sm <- SMOTE(
  X = df_train_ohe %>% dplyr::select(-fuga),
  target = y_train_num,
  K = 5,
  dup_size = 2
)
df_train_bal <- sm$data %>%
  rename(fuga = class) %>%
  mutate(fuga = factor(fuga, levels = c("0","1")))
X_tr_mat <- as.matrix(df_train_bal %>% dplyr::select(-fuga))
y_tr_num <- as.numeric(as.character(df_train_bal$fuga))

# Entrenament models
set.seed(42)
cv_glm <- cv.glmnet(
  x = X_tr_mat,
  y = y_tr_num,
  family = "binomial",
  alpha = 1,
  nfolds = 5
)
mod_glm <- glmnet(
  x = X_tr_mat,
  y = y_tr_num,
  family = "binomial",
  alpha = 1,
  lambda = cv_glm$lambda.min
)

mod_rf <- ranger(
  formula      = fuga ~ .,
  data         = df_train_bal,
  probability  = TRUE,
  num.trees    = 500,
  importance   = "permutation"
)

dtrain <- xgb.DMatrix(data = X_tr_mat, label = y_tr_num)
dtest  <- xgb.DMatrix(data = as.matrix(X_test), label = y_test_num)
set.seed(42)
cv_xgb <- xgb.cv(
  params               = list(objective="binary:logistic", eval_metric="auc"),
  data                 = dtrain,
  nrounds              = 100,
  nfold                = 5,
  early_stopping_rounds= 10,
  verbose              = FALSE
)
mod_xgb <- xgb.train(
  params   = list(objective="binary:logistic", eval_metric="auc"),
  data     = dtrain,
  nrounds  = cv_xgb$best_iteration
)

# Prediccions sobre el conjunt de test
preds <- tibble(
  prob_glm = as.numeric(predict(mod_glm, newx = as.matrix(X_test), type = "response")),
  prob_rf  = predict(mod_rf, data = X_test, type = "response")$predictions[,"1"],
  prob_xgb = predict(mod_xgb, dtest),
  truth    = y_test_num
)

# Predictors seleccionats
predictors <- c(
  "presion_real","presion_sensor","flujo_entrada_tuberia","flujo_salida_tuberia",
  "delta_flujo_nodo","delta_presion_media_vecinos","centralidad","distancia_fuente",
  "grado_conectividad","es_consumo_discreto","hora_dia","log_presion_real",
  "std_presion_vecinos","zscore_presion_local"
)

# Dataframe reduït
df2 <- df %>% select(fuga, all_of(predictors))

# Càlcul del VIF
mod_vif <- lm(fuga ~ ., data = df2)
vif_vals <- vif(mod_vif)
print(vif_vals)
high_vif <- vif_vals[vif_vals > 5]
print(high_vif)

# Anàlisi de Components Principals (PCA)
X_scaled <- scale(df2 %>% dplyr::select(-fuga))
pca <- prcomp(X_scaled, center = TRUE, scale. = FALSE)
summary(pca)
plot(pca, type = "l")


# Llegir el dataset i transformar variables
df <- read_csv("dataset_fugas_gas_realista_v4.csv") %>%
  mutate(
    fuga = factor(fuga, levels = c(0,1)),
    hora_cat = case_when(
      hora_dia %in% 0:5   ~ "madrugada",
      hora_dia %in% 6:11  ~ "matí",
      hora_dia %in% 12:17 ~ "tarda",
      hora_dia %in% 18:23 ~ "nit"
    ) %>% factor(levels = c("madrugada","matí","tarda","nit"))
  )

# Variables predictives
vars_pred <- c(
  "delta_flujo_nodo",
  "flujo_entrada_tuberia",
  "zscore_presion_local",
  "delta_presion_media_vecinos",
  "centralidad",
  "flujo_salida_tuberia",
  "grado_conectividad",
  "distancia_fuente",
  "es_consumo_discreto",
  "std_presion_vecinos",
  "hora_cat"
)

# Divisió Train/Test
set.seed(42)
idx_train <- createDataPartition(df$fuga, p = 0.7, list = FALSE)
df_train <- df[idx_train, ]
df_test  <- df[-idx_train, ]

# Codificació One-Hot de les categories a train i test
ohe <- dummyVars(~ ., data = df_train[vars_pred], fullRank = TRUE)
X_train <- predict(ohe, df_train[vars_pred]) %>% as.data.frame()
X_test  <- predict(ohe, df_test[vars_pred])  %>% as.data.frame()

# Adjuntar la variable resposta
df_train_ohe <- X_train %>% mutate(fuga = df_train$fuga)
y_train_num  <- as.numeric(as.character(df_train_ohe$fuga))
y_test_num   <- as.numeric(as.character(df_test$fuga))

# Reequilibri amb SMOTE sobre el conjunt d'entrenament
sm <- SMOTE(
  X = df_train_ohe %>% dplyr::select(-fuga),
  target = y_train_num,
  K = 5,
  dup_size = 2
)
df_train_bal <- sm$data %>%
  rename(fuga = class) %>%
  mutate(fuga = factor(fuga, levels = c("0","1")))
X_tr_mat <- as.matrix(df_train_bal %>% dplyr::select(-fuga))
y_tr_num <- as.numeric(as.character(df_train_bal$fuga))

# Entrenament models
set.seed(42)
cv_glm <- cv.glmnet(
  x = X_tr_mat,
  y = y_tr_num,
  family = "binomial",
  alpha = 1,
  nfolds = 5
)
mod_glm <- glmnet(
  x = X_tr_mat,
  y = y_tr_num,
  family = "binomial",
  alpha = 1,
  lambda = cv_glm$lambda.min
)

mod_rf <- ranger(
  formula      = fuga ~ .,
  data         = df_train_bal,
  probability  = TRUE,
  num.trees    = 500,
  importance   = "permutation"
)

dtrain <- xgb.DMatrix(data = X_tr_mat, label = y_tr_num)
dtest  <- xgb.DMatrix(data = as.matrix(X_test), label = y_test_num)
set.seed(42)
cv_xgb <- xgb.cv(
  params               = list(objective="binary:logistic", eval_metric="auc"),
  data                 = dtrain,
  nrounds              = 100,
  nfold                = 5,
  early_stopping_rounds= 10,
  verbose              = FALSE
)
mod_xgb <- xgb.train(
  params   = list(objective="binary:logistic", eval_metric="auc"),
  data     = dtrain,
  nrounds  = cv_xgb$best_iteration
)

# Prediccions sobre el conjunt de test
preds <- tibble(
  prob_glm = as.numeric(predict(mod_glm, newx = as.matrix(X_test), type = "response")),
  prob_rf  = predict(mod_rf, data = X_test, type = "response")$predictions[,"1"],
  prob_xgb = predict(mod_xgb, dtest),
  truth    = y_test_num
)

# Càlcul mètriques
metrics <- preds %>%
  summarise(
    AUC_glm   = as.numeric(roc(truth, prob_glm)$auc),
    AUC_rf    = as.numeric(roc(truth, prob_rf)$auc),
    AUC_xgb   = as.numeric(roc(truth, prob_xgb)$auc),
    Brier_glm = mean((prob_glm - truth)^2),
    Brier_rf  = mean((prob_rf  - truth)^2),
    Brier_xgb = mean((prob_xgb - truth)^2),
    Accuracy_glm = mean((prob_glm >= 0.5) == truth),
    Accuracy_rf  = mean((prob_rf  >= 0.5) == truth),
    Accuracy_xgb = mean((prob_xgb >= 0.5) == truth),
    LogLoss_glm = -mean(truth*log(prob_glm) + (1-truth)*log(1-prob_glm)),
    LogLoss_rf  = -mean(truth*log(prob_rf)  + (1-truth)*log(1-prob_rf)),
    LogLoss_xgb = -mean(truth*log(prob_xgb)+ (1-truth)*log(1-prob_xgb)),
    Sensitivity_glm = sum((prob_glm >= 0.5) & (truth == 1)) / sum(truth == 1),
    Sensitivity_rf  = sum((prob_rf  >= 0.5) & (truth == 1)) / sum(truth == 1),
    Sensitivity_xgb = sum((prob_xgb >= 0.5) & (truth == 1)) / sum(truth == 1),
    Specificity_glm = sum((prob_glm <  0.5) & (truth == 0)) / sum(truth == 0),
    Specificity_rf  = sum((prob_rf  <  0.5) & (truth == 0)) / sum(truth == 0),
    Specificity_xgb = sum((prob_xgb <  0.5) & (truth == 0)) / sum(truth == 0),
    Precision_glm   = sum((prob_glm >= 0.5) & (truth == 1)) / sum(prob_glm >= 0.5),
    Precision_rf    = sum((prob_rf  >= 0.5) & (truth == 1)) / sum(prob_rf  >= 0.5),
    Precision_xgb   = sum((prob_xgb >= 0.5) & (truth == 1)) / sum(prob_xgb >= 0.5)
  ) %>%
  mutate(
    F1_glm = 2 * (Precision_glm * Sensitivity_glm) / (Precision_glm + Sensitivity_glm),
    F1_rf  = 2 * (Precision_rf  * Sensitivity_rf)  / (Precision_rf  + Sensitivity_rf),
    F1_xgb = 2 * (Precision_xgb * Sensitivity_xgb) / (Precision_xgb + Sensitivity_xgb)
  ) %>%
  pivot_longer(everything(), names_to = c("metric","model"), names_sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value)

print(metrics)
cal_data <- preds %>%
  pivot_longer(cols = starts_with("prob_"), names_prefix = "prob_", 
               names_to = "model", values_to = "pred") %>%
  mutate(obs = factor(truth, levels = c(0, 1), labels = c("NoFuga", "Fuga")))

# Generar objecte de calibratge amb facetes per model
tab_cal <- caret::calibration(obs ~ pred | model, data = cal_data, 
                              class = "Fuga", cuts = 10)

# Gràfic amb lattice sense utilitzar '+'
lattice::xyplot(tab_cal,
                auto.key = list(columns = 2),
                xlab = "Probabilitat Predita",
                ylab = "Freqüència Observada",
                main = "Corbes de Calibratge per Model")

# Calcular valors SHAP
shap_values <- predict(mod_xgb, dtest, predcontrib = TRUE)

# Calcular importància SHAP (mitjana de valors absoluts, excloent l'intercepte)
shap_importance <- colMeans(abs(shap_values[, -ncol(shap_values)]))
shap_importance_df <- data.frame(
  Importance = shap_importance
) %>% arrange(desc(Importance))

print(shap_importance_df)


# Dividir df_train_bal en train_main / df_cal (20% per a calibratge)
set.seed(42)
idx_cal       <- createDataPartition(df_train_bal$fuga, p = 0.2, list = FALSE)
df_cal        <- df_train_bal[idx_cal, ]
df_train_main <- df_train_bal[-idx_cal, ]

# Construir matriu de predictors per al calibratge
X_cal <- df_cal %>% dplyr::select(-fuga)
X_cal_mat <- as.matrix(X_cal)
y_cal     <- as.numeric(as.character(df_cal$fuga))

# Obtenir probabilitats crues sobre df_cal
prob_glm_crudo <- as.numeric(predict(mod_glm,
                                     newx = X_cal_mat,
                                     type = "response")[,1])
prob_rf_crudo  <- predict(mod_rf,
                          data = X_cal,
                          type = "response")$predictions[,"1"]
prob_xgb_crudo <- predict(mod_xgb,
                          xgb.DMatrix(data = X_cal_mat))

cal_preds <- tibble(
  truth    = y_cal,
  prob_glm = prob_glm_crudo,
  prob_rf  = prob_rf_crudo,
  prob_xgb = prob_xgb_crudo
)

# Ajustar calibradors de Platt (regressió logística)
cal_glm <- glm(truth ~ prob_glm, data = cal_preds, family = binomial)
cal_rf  <- glm(truth ~ prob_rf,  data = cal_preds, family = binomial)
cal_xgb <- glm(truth ~ prob_xgb, data = cal_preds, family = binomial)

# Aplicar calibradors al conjunt de test
test_preds <- preds %>%
  mutate(
    prob_glm_cal = predict(cal_glm, newdata = tibble(prob_glm = prob_glm), type = "response"),
    prob_rf_cal  = predict(cal_rf,  newdata = tibble(prob_rf  = prob_rf),  type = "response"),
    prob_xgb_cal = predict(cal_xgb, newdata = tibble(prob_xgb = prob_xgb), type = "response")
  )

# Funció per calcular mètriques
compute_metrics <- function(truth, prob){
  roc_obj <- roc(truth, prob)
  auc     <- as.numeric(roc_obj$auc)
  brier   <- mean((prob - truth)^2)
  logloss <- -mean(truth * log(prob) + (1 - truth) * log(1 - prob))
  acc     <- mean((prob >= 0.5) == truth)
  sens    <- sum((prob >= 0.5 & truth == 1)) / sum(truth == 1)
  spec    <- sum((prob <  0.5 & truth == 0)) / sum(truth == 0)
  prec    <- sum((prob >= 0.5 & truth == 1)) / sum(prob >= 0.5)
  f1      <- 2 * prec * sens / (prec + sens)
  tibble(AUC = auc, Brier = brier, LogLoss = logloss,
         Accuracy = acc, Sensitivity = sens,
         Specificity = spec, Precision = prec,
         F1 = f1)
}


# Calcular mètriques abans i després
metrics_before <- bind_rows(
  glm  = compute_metrics(preds$truth,     preds$prob_glm),
  rf   = compute_metrics(preds$truth,     preds$prob_rf),
  xgb  = compute_metrics(preds$truth,     preds$prob_xgb),
  .id = "model"
)

metrics_after <- bind_rows(
  glm  = compute_metrics(test_preds$truth, test_preds$prob_glm_cal),
  rf   = compute_metrics(test_preds$truth, test_preds$prob_rf_cal),
  xgb  = compute_metrics(test_preds$truth, test_preds$prob_xgb_cal),
  .id = "model"
)

# Comparativa final de mètriques abans/després de la calibració amb regressió logística
results <- metrics_before %>%
  rename_with(~ paste0(.,"_before"), -model) %>%
  left_join(
    metrics_after %>% rename_with(~ paste0(.,"_after"), -model),
    by = "model"
  )

print(results)



# Dividir df_train_bal en train_main / df_cal (20 % per a calibratge)
set.seed(42)
idx_cal       <- createDataPartition(df_train_bal$fuga, p = 0.20, list = FALSE)
df_cal        <- df_train_bal[idx_cal, ]
df_train_main <- df_train_bal[-idx_cal, ]

# Obtenir probabilitats “crues” sobre df_cal
X_cal     <- df_cal %>% dplyr::select(-fuga)
X_cal_mat <- as.matrix(X_cal)
y_cal     <- as.numeric(as.character(df_cal$fuga))

prob_glm_crudo <- as.numeric(
  predict(mod_glm, newx = X_cal_mat, type = "response")[,1]
)
prob_rf_crudo  <- predict(
  mod_rf, data = X_cal, type = "response"
)$predictions[,"1"]
prob_xgb_crudo <- predict(
  mod_xgb, xgb.DMatrix(data = X_cal_mat)
)

cal_preds <- tibble(
  truth    = y_cal,
  prob_glm = prob_glm_crudo,
  prob_rf  = prob_rf_crudo,
  prob_xgb = prob_xgb_crudo
)

# Ajustar calibradors mitjançant Regressió Isotònica
iso_glm <- stats::isoreg(cal_preds$prob_glm, cal_preds$truth)
iso_rf  <- stats::isoreg(cal_preds$prob_rf,  cal_preds$truth)
iso_xgb <- stats::isoreg(cal_preds$prob_xgb, cal_preds$truth)

# Aplicar calibradors al conjunt de test
preds_iso <- preds %>%
  mutate(
    prob_glm_iso = approx(x = iso_glm$x, y = iso_glm$yf,
                          xout = prob_glm, rule = 2)$y,
    prob_rf_iso  = approx(x = iso_rf$x,  y = iso_rf$yf,
                          xout = prob_rf,  rule = 2)$y,
    prob_xgb_iso = approx(x = iso_xgb$x, y = iso_xgb$yf,
                          xout = prob_xgb, rule = 2)$y
  )

# Funció per calcular mètriques
compute_metrics <- function(truth, prob){
  roc_obj <- roc(truth, prob)
  auc     <- as.numeric(roc_obj$auc)
  brier   <- mean((prob - truth)^2)
  logloss <- -mean(truth * log(prob) + (1 - truth) * log(1 - prob))
  acc     <- mean((prob >= 0.5) == truth)
  sens    <- sum((prob >= 0.5 & truth == 1)) / sum(truth == 1)
  spec    <- sum((prob <  0.5 & truth == 0)) / sum(truth == 0)
  prec    <- sum((prob >= 0.5 & truth == 1)) / sum(prob >= 0.5)
  f1      <- 2 * prec * sens / (prec + sens)
  tibble(AUC = auc,
         Brier = brier,
         LogLoss = logloss,
         Accuracy = acc,
         Sensitivity = sens,
         Specificity = spec,
         Precision = prec,
         F1 = f1)
}

# Mètriques abans i després del calibratge isotònic
metrics_before <- bind_rows(
  glm = compute_metrics(preds$truth,     preds$prob_glm),
  rf  = compute_metrics(preds$truth,     preds$prob_rf),
  xgb = compute_metrics(preds$truth,     preds$prob_xgb),
  .id = "model"
)

metrics_after <- bind_rows(
  glm = compute_metrics(preds_iso$truth, preds_iso$prob_glm_iso),
  rf  = compute_metrics(preds_iso$truth, preds_iso$prob_rf_iso),
  xgb = compute_metrics(preds_iso$truth, preds_iso$prob_xgb_iso),
  .id = "model"
)

# Comparativa final
results <- metrics_before %>%
  rename_with(~ paste0(.,"_before"), -model) %>%
  left_join(
    metrics_after %>% rename_with(~ paste0(.,"_after"), -model),
    by = "model"
  )

print(results)



# 2) Funció robusta: Isotònica (gpava) + spline suavitzat amb fallback
smooth_gpava_spline_safe <- function(raw_train, y_train, raw_test, spar = 0.5) {
  # 2.1) Ajust isotònic pur
  iso_fit <- gpava(z = raw_train, y = y_train)
  xs <- iso_fit$x
  ys <- iso_fit$y
  
  # 2.2) Neteja: valors únics, finits i ordenats
  idx <- which(!duplicated(xs) & is.finite(xs) & is.finite(ys))
  xs <- xs[idx]; ys <- ys[idx]
  o  <- order(xs); xs <- xs[o]; ys <- ys[o]
  
  # 2.3) Si hi ha <3 punts, només interpolar (sense spline)
  if(length(xs) < 3) {
    pred <- approx(x = xs, y = ys, xout = raw_test, rule = 2)$y
    return(pmin(pmax(pred, 0), 1))
  }
  
  # 2.4) Intentar suavitzar amb spline; si falla, utilitzar interpolació
  out <- tryCatch({
    spline_fit <- smooth.spline(x = xs, y = ys, spar = spar)
    pred <- predict(spline_fit, raw_test)$y
    pmin(pmax(pred, 0), 1)
  }, error = function(e) {
    approx(x = xs, y = ys, xout = raw_test, rule = 2)$y %>%
      { pmin(pmax(., 0), 1) }
  })
  
  return(out)
}

# 3) Supòsits: ja tens cal_preds (truth, prob_glm, prob_rf, prob_xgb)
#              i preds      (truth, prob_glm, prob_rf, prob_xgb)

# 4) Calibrar amb la funció segura
prob_glm_iso2 <- smooth_gpava_spline_safe(
  raw_train = cal_preds$prob_glm,
  y_train   = cal_preds$truth,
  raw_test  = preds$prob_glm,
  spar      = 0.5
)
prob_rf_iso2  <- smooth_gpava_spline_safe(
  raw_train = cal_preds$prob_rf,
  y_train   = cal_preds$truth,
  raw_test  = preds$prob_rf,
  spar      = 0.5
)
prob_xgb_iso2 <- smooth_gpava_spline_safe(
  raw_train = cal_preds$prob_xgb,
  y_train   = cal_preds$truth,
  raw_test  = preds$prob_xgb,
  spar      = 0.5
)

# 5) Funció per calcular mètriques
compute_metrics <- function(truth, prob) {
  roc_obj <- roc(truth, prob)
  tibble(
    AUC         = as.numeric(roc_obj$auc),
    Brier       = mean((prob - truth)^2),
    LogLoss     = -mean(truth * log(prob) + (1 - truth) * log(1 - prob)),
    Accuracy    = mean((prob >= 0.5) == truth),
    Sensitivity = sum((prob >= 0.5 & truth == 1)) / sum(truth == 1),
    Specificity = sum((prob <  0.5 & truth == 0)) / sum(truth == 0),
    Precision   = sum((prob >= 0.5 & truth == 1)) / sum(prob >= 0.5),
    F1          = {
      prec <- sum((prob >= 0.5 & truth == 1)) / sum(prob >= 0.5)
      sens <- sum((prob >= 0.5 & truth == 1)) / sum(truth == 1)
      2 * prec * sens / (prec + sens)
    }
  )
}

# 6) Mètriques abans i després del calibratge
metrics_before <- bind_rows(
  glm = compute_metrics(preds$truth,     preds$prob_glm),
  rf  = compute_metrics(preds$truth,     preds$prob_rf),
  xgb = compute_metrics(preds$truth,     preds$prob_xgb),
  .id = "model"
)

metrics_after_iso2 <- bind_rows(
  glm = compute_metrics(preds$truth,     prob_glm_iso2),
  rf  = compute_metrics(preds$truth,     prob_rf_iso2),
  xgb = compute_metrics(preds$truth,     prob_xgb_iso2),
  .id = "model"
)

# 7) Comparativa final de resultats abans i després del calibratge robust
results_iso2 <- metrics_before %>%
  rename_with(~ paste0(.,"_before"), -model) %>%
  left_join(
    metrics_after_iso2 %>% rename_with(~ paste0(.,"_after"), -model),
    by = "model"
  )

print(results_iso2)

library(caret)
library(PRROC)
library(dplyr)
library(ggplot2)

# Funció per calcular corba Precision–Recall + llindar òptim (F1)
pr_with_thresh <- function(truth, prob, model_name){
  pr <- pr.curve(
    scores.class0 = prob[truth == 1],
    scores.class1 = prob[truth == 0],
    curve = TRUE
  )
  df_pr <- data.frame(
    Recall    = pr$curve[,1],
    Precision = pr$curve[,2],
    Threshold = pr$curve[,3]
  )
  # Llindar que maximitza F1
  df_pr <- df_pr %>%
    mutate(F1 = 2*Precision*Recall/(Precision+Recall)) 
  best <- df_pr[which.max(df_pr$F1), ]
  best$Model <- model_name
  
  list(curve = df_pr, best = best)
}

# Calcular per a cada model
pr_glm <- pr_with_thresh(preds$truth, preds$prob_glm, "GLM")
pr_rf  <- pr_with_thresh(preds$truth, preds$prob_rf,  "RF")
pr_xgb <- pr_with_thresh(preds$truth, preds$prob_xgb, "XGB")

# Unir corbes i generar gràfic
all_pr <- bind_rows(
  pr_glm$curve %>% mutate(Model="GLM"),
  pr_rf$curve  %>% mutate(Model="RF"),
  pr_xgb$curve %>% mutate(Model="XGB")
)

ggplot(all_pr, aes(x = Recall, y = Precision, color = Model)) +
  geom_line(linewidth = 1) +
  labs(title="Corbes Precision–Recall", x="Recall", y="Precisió") +
  theme_minimal()

# Combinar les files dels millors llindars i extreure columnes de forma robusta
best_thresholds <- bind_rows(
  pr_glm$best,
  pr_rf$best,
  pr_xgb$best
) %>% 
  as_tibble()

thresholds_tbl <- best_thresholds %>%
  dplyr::select(Model, Threshold, Precision, Recall, F1)

print(thresholds_tbl)


# Reconstruir un objecte xgb.Booster amb noms de variables
mod_xgb$feature_names <- colnames(X_train)

# Funció ràpida de PDP (Partial Dependence Plot)
plot_pdp <- function(model, data_train, var){
  pd <- partial(
    object    = model, 
    pred.var  = var, 
    train     = data_train,
    prob      = TRUE,
    grid.resolution = 30
  )
  autoplot(pd) + 
    labs(title = paste("PDP per", var), y="P(p[fuga]=1)") +
    theme_minimal()
}

# Prepara un dataframe amb les mateixes columnes usades per XGB
dtrain_df <- as.data.frame(X_train)

# Gràfics de dependència parcial
plot_pdp(mod_xgb, dtrain_df, "delta_flujo_nodo")
plot_pdp(mod_xgb, dtrain_df, "flujo_entrada_tuberia")
plot_pdp(mod_xgb, dtrain_df, "zscore_presion_local")

library(dplyr)
library(cluster)
library(factoextra)

# 1) Agregació per node
df_node <- df %>%
  group_by(nodo) %>%
  summarise(
    delta_flujo  = mean(delta_flujo_nodo, na.rm=TRUE),
    flujo_ent    = mean(flujo_entrada_tuberia, na.rm=TRUE),
    zscore_loc   = mean(zscore_presion_local, na.rm=TRUE),
    delta_pres   = mean(delta_presion_media_vecinos, na.rm=TRUE),
    centralidad  = mean(centralidad, na.rm=TRUE),
    flujo_sal    = mean(flujo_salida_tuberia, na.rm=TRUE),
    grado        = mean(grado_conectividad, na.rm=TRUE),
    distancia    = mean(distancia_fuente, na.rm=TRUE),
    consumo_disc = mean(es_consumo_discreto, na.rm=TRUE),
    std_pres     = mean(std_presion_vecinos, na.rm=TRUE),
    p_fuga       = mean(fuga, na.rm=TRUE)
  ) %>%
  ungroup()

# 2) Construir data.frame només numèric i eliminar columnes no finites
num_df <- df_node %>% dplyr::select(-nodo)
num_df <- num_df[, sapply(num_df, function(col) all(is.finite(col)))]

# 3) Escalar
scaled_df <- scale(num_df)

# 4) K-means
set.seed(42)
km <- kmeans(scaled_df, centers = 3, nstart = 50)

# 5) Visualitzar els clústers
fviz_cluster(km, data = scaled_df,
             ellipse.type = "convex",
             repel = TRUE) +
  labs(title = "Clustering de nodes (k = 3)")

# 6) Perfil mitjà de cada clúster
df_node %>%
  mutate(cluster = factor(km$cluster)) %>%
  group_by(cluster) %>%
  summarise(across(-nodo, mean)) %>%
  print()


# Error original
df_err <- df %>%
  mutate(err_orig = presion_sensor - presion_real)

# Mostra la distribució de l'error original
summary(df_err$err_orig)

# Factors d’escala per simular menys/més soroll
scales <- c(0.5, 1, 2)

# Inicialitza tibble de resultats
res <- tibble(scale = numeric(), RMSE = numeric(), MAE = numeric(), Bias = numeric())

# Bucle que aplica cada factor i calcula mètriques
for(s in scales) {
  df2 <- df_err %>%
    mutate(err_scaled = err_orig * s)
  
  rmse <- sqrt(mean(df2$err_scaled^2))
  mae  <- mean(abs(df2$err_scaled))
  bias <- mean(df2$err_scaled)
  
  res <- bind_rows(res, tibble(scale = s, RMSE = rmse, MAE = mae, Bias = bias))
}

print(res)

# Gràfic dels resultats
res_long <- res %>%
  pivot_longer(cols = -scale, names_to = "metric", values_to = "value")

ggplot(res_long, aes(x = factor(scale), y = value, group = metric, color = metric)) +
  geom_point(size = 3) +
  geom_line(aes(group = metric), linewidth = 1) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Sensibilitat de la lectura de pressió al nivell de soroll",
    x = "Factor d'escala del soroll",
    y = "Valor de la mètrica"
  ) +
  theme_minimal()


# Carregar dades
df <- read_csv("dataset_fugas_gas_realista_v4.csv") %>%
  mutate(fuga = factor(fuga, levels=c(0,1)))

# Variables essencials
vars_esenciales <- c(
  "delta_flujo_nodo",
  "flujo_entrada_tuberia",
  "zscore_presion_local",
  "flujo_salida_tuberia"
)

# Separació Train/Test
set.seed(42)
idx <- createDataPartition(df$fuga, p=0.7, list=FALSE)
train_raw <- df[idx, ]
test_raw  <- df[-idx, ]

# Aplicar SMOTE al conjunt d'entrenament
Xtr_df <- train_raw %>% dplyr::select(all_of(vars_esenciales))
ytr    <- as.numeric(as.character(train_raw$fuga))
sm     <- SMOTE(X=Xtr_df, target=ytr, K=5, dup_size=2)
train_bal <- sm$data %>%
  rename(fuga = class) %>%
  mutate(fuga = factor(fuga, levels=c("0","1")))

# Preparar matrius per a XGBoost
X_tr <- as.matrix(train_bal %>% dplyr::select(-fuga))
y_tr <- as.numeric(as.character(train_bal$fuga))
dtrain <- xgb.DMatrix(data=X_tr, label=y_tr)

X_te <- as.matrix(test_raw %>% dplyr::select(all_of(vars_esenciales)))
y_te <- as.numeric(as.character(test_raw$fuga))
dtest <- xgb.DMatrix(data=X_te, label=y_te)

# Validació creuada
set.seed(42)
cv <- xgb.cv(
  params = list(objective="binary:logistic", eval_metric="auc"),
  data   = dtrain,
  nfold  = 5,
  nrounds = 100,
  early_stopping_rounds = 10,
  verbose = FALSE
)

# Entrenar el model final
best_n <- cv$best_iteration
mod_xgb_small <- xgb.train(
  params   = list(objective="binary:logistic", eval_metric="auc"),
  data     = dtrain,
  nrounds  = best_n
)

# Prediccions i mètriques
preds <- predict(mod_xgb_small, dtest)
roc_obj <- roc(y_te, preds)

metrics <- tibble(
  AUC     = as.numeric(roc_obj$auc),
  Brier   = mean((preds - y_te)^2),
  LogLoss = -mean(y_te*log(preds) + (1-y_te)*log(1-preds)),
  Accuracy    = mean((preds>=0.5) == y_te),
  Sensitivity = sum((preds>=0.5) & (y_te==1))/sum(y_te==1),
  Specificity = sum((preds< 0.5) & (y_te==0))/sum(y_te==0),
  Precision   = sum((preds>=0.5) & (y_te==1))/sum(preds>=0.5),
  F1          = {
    prec <- sum((preds>=0.5)&(y_te==1))/sum(preds>=0.5)
    rec  <- sum((preds>=0.5)&(y_te==1))/sum(y_te==1)
    2*prec*rec/(prec+rec)
  }
)

print(metrics)


# Corba ROC
plot.roc(roc_obj, main="ROC - XGBoost (4 variables essencials)")

# Exemple: buscar t a [0,1] que minimitzi el cost
cost_fun <- function(t, truth, prob){
  preds <- ifelse(prob>=t,1,0)
  fp <- sum(preds==1 & truth==0)
  fn <- sum(preds==0 & truth==1)
  100*fn + 5*fp
}

ts <- seq(0,1,by=0.01)
costs <- sapply(ts, cost_fun, truth=y_te, prob=preds)
best_t <- ts[which.min(costs)]
plot(ts, costs, type="l", main="Cost vs llindar", xlab="t", ylab="Cost")
abline(v=best_t, col="red")
best_t


# Llindar òptim
best_t <- 0.07
pred_class <- ifelse(preds >= best_t, 1, 0)

# Matriu de confusió simple
conf_mat <- table(Predit = pred_class, Veritable = y_te)
print(conf_mat)

# Mètriques amb caret::confusionMatrix
# Convertim a factors (positiu = "1")
cm <- confusionMatrix(
  factor(pred_class, levels = c(0,1)),
  factor(y_te,       levels = c(0,1)),
  positive = "1"
)
print(cm)

# Mètriques addicionals
roc_obj <- roc(y_te, preds)
metrics <- tibble(
  AUC         = as.numeric(roc_obj$auc),
  Precisió    = cm$overall["Accuracy"],
  Sensibilitat = cm$byClass["Sensitivity"],
  Especificitat = cm$byClass["Specificity"],
  Precisió   = cm$byClass["Precision"],
  F1          = cm$byClass["F1"]
)
print(metrics)


# Prediccions de probabilitat sobre el conjunt de test
prob_xgb_small <- predict(mod_xgb_small, dtest)

# Convertir a classes utilitzant el llindar òptim
best_t <- 0.07
pred_small_class <- factor(
  if_else(prob_xgb_small >= best_t, "Fuga", "NoFuga"),
  levels = c("NoFuga","Fuga")
)

# Valors reals com a factor
truth_small <- factor(
  y_te,
  levels = c(0,1),
  labels = c("NoFuga","Fuga")
)

# Matriu de confusió i mètriques
cm_small <- confusionMatrix(pred_small_class, truth_small, positive = "Fuga")
print(cm_small)

precision_small <- cm_small$byClass["Precision"]
recall_small    <- cm_small$byClass["Sensitivity"]
f1_small        <- F1_Score(y_true = truth_small,
                            y_pred = pred_small_class,
                            positive = "Fuga")

cat("\nPrecisió =", precision_small,
    "\nSensibilitat =", recall_small,
    "\nF1-score  =", f1_small, "\n")

# Corbes ROC i Precisió–Recall
roc_small <- roc(truth_small, prob_xgb_small)
auc_small <- as.numeric(roc_small$auc)
cat("\nAUC-ROC =", round(auc_small, 4), "\n")

pr_small <- pr.curve(
  scores.class0 = prob_xgb_small[truth_small=="Fuga"],
  scores.class1 = prob_xgb_small[truth_small=="NoFuga"],
  curve = TRUE
)

# Gràfica ROC
plot(roc_small, main = "ROC - XGB essencial")

# Gràfica PR
plot(pr_small, main = "PR - XGB essencial")

# Anàlisi d'errors
results_small <- tibble(
  prob  = prob_xgb_small,
  pred  = pred_small_class,
  truth = truth_small
)

fn_small <- results_small %>% filter(truth=="Fuga",   pred=="NoFuga")
fp_small <- results_small %>% filter(truth=="NoFuga", pred=="Fuga")

cat("\nFalsos negatius:", nrow(fn_small),
    "\nFalsos positius:", nrow(fp_small), "\n")

# Veure els primers casos de cada tipus d'error
head(fn_small)
head(fp_small)

X_te_small <- test_raw %>%
  dplyr::select(delta_flujo_nodo,
                flujo_entrada_tuberia,
                zscore_presion_local,
                flujo_salida_tuberia) %>%
  as.matrix()

dtest_small <- xgb.DMatrix(data = X_te_small)

prob_node <- predict(mod_xgb_small, dtest_small)

test_scored <- test_raw %>%
  mutate(
    prob_fuga  = prob_node,
    pred_fuga  = ifelse(prob_fuga >= best_t, 1L, 0L)
  )

top_events <- test_scored[order(-test_scored$prob_fuga), 
                          c("nodo","prob_fuga","pred_fuga", 
                            names(test_scored)[-(1:3)])][1:20, ]

head(top_events, 50)

# Resum per node
node_summary <- test_scored %>%
  dplyr::as_tibble() %>%
  group_by(nodo) %>%
  summarise(
    n_instants = n(),
    n_alertes   = sum(pred_fuga),
    pct_alertes = mean(pred_fuga),
    prob_mitjana  = mean(prob_fuga),
    prob_max    = max(prob_fuga)
  ) %>%
  arrange(desc(prob_mitjana)) %>%
  slice_head(n = 10)

print(node_summary)


font_add(family = "Arial", regular = "Arial.ttf") 
showtext_auto()
theme_set(theme_minimal(base_family = "Arial"))

# Dades de canonades
pipes_data <- tribble(
  ~pipe, ~from, ~to,  ~length_m, ~diam_in,
  ...
)

vars_esenciales <- c(
  "delta_flujo_nodo",
  "flujo_entrada_tuberia",
  "zscore_presion_local",
  "flujo_salida_tuberia"
)
best_t <- 0.07

# Generar prediccions sobre el conjunt de test 
df_test_pred <- test_raw %>%
  select(nodo, all_of(vars_esenciales)) %>%
  mutate(
    prob_fuga = predict(mod_xgb_small,
                        xgb.DMatrix(as.matrix(select(., all_of(vars_esenciales))))),
    pred_fuga = as.integer(prob_fuga >= best_t)
  )

# Resum per node
node_summary <- df_test_pred %>%
  group_by(nodo) %>%
  summarise(
    n_instants = n(),
    n_alertes   = sum(pred_fuga),
    pct_alertes = n_alertes / n_instants,
    prob_mitjana  = mean(prob_fuga),
    prob_max    = max(prob_fuga),
    .groups = "drop"
  )

all_nodes <- tibble(name = sort(unique(c(pipes_data$from, pipes_data$to))))
nodes     <- all_nodes %>%
  left_join(node_summary %>% rename(name = nodo), by = "name") %>%
  replace_na(list(
    n_instants = 0,
    n_alertes   = 0,
    pct_alertes = 0,
    prob_mitjana = 0,
    prob_max    = 0
  ))

g_igraph <- graph_from_data_frame(
  d = pipes_data %>% select(from, to),
  vertices = nodes,
  directed = FALSE
)
graph <- as_tbl_graph(g_igraph)

set.seed(42)
lay <- create_layout(graph, layout = "fr")

ggraph(lay) +
  geom_edge_link(color = "grey80", alpha = 0.6) +
  geom_node_point(aes(size = pct_alertes, color = prob_mitjana)) +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 size = 3, 
                 color = "black") +
  scale_size_continuous(range = c(2, 8), name = "% alertes") +
  scale_color_viridis_c(name = "Prob. mitjana", option = "C") +
  theme_graph() +
  labs(
    title    = "Mapa de risc de fuita en la xarxa"
  )
