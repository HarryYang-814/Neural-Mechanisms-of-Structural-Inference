library(lmerTest)
library(lme4)
library(dplyr)
library(readxl)
library(effectsize)
library(emmeans)
library(car)
library(ggplot2)
library(ggsignif)
library(ggpubr)     
library(patchwork)  
###### loading reflexives data, 31 ID
data_raw <- read_excel("G:/Research Data/1_Projects/2_Merge and composition/Data and analysis/data_analysis_revised/data/0_behavioural/behavioural data.xlsx",sheet = 1) 

##response accuracy for all items 
## (NP-correct/ungrammatical/unpragmatic,
## VP-correct/ungrammatical/unpragmatic,
## control-correct/ungrammatical/unpragmatic)

rawq <- data_raw
rawq$Accuracy < as.numeric(rawq$Accuracy)
accuracy_means <- aggregate(Accuracy~Subject, data = rawq, FUN = mean)
range(accuracy_means$Accuracy)
## mean accuracy by participants: 80.7%-98.2%

##mean accuracy for each condition

raw.NP <- subset(data_raw, data_raw$Condition == "NP") 
raw.VP <- subset(data_raw, data_raw$Condition == "VP") 
raw.CC <- subset(data_raw, data_raw$Condition == "CC") 

raw.NP$Accuracy <- as.numeric(raw.NP$Accuracy)
raw.VP$Accuracy <- as.numeric(raw.VP$Accuracy)
raw.CC$Accuracy <- as.numeric(raw.CC$Accuracy)

raw.NP_acc <- aggregate(Accuracy~Subject, data = raw.NP, FUN = mean)
range(raw.NP_acc$Accuracy) ## 73.7-100% 
raw.VP_acc <- aggregate(Accuracy~Subject, data = raw.VP, FUN = mean)
range(raw.VP_acc$Accuracy) ## 71.1-100% 
raw.CC_acc <- aggregate(Accuracy~Subject, data = raw.CC, FUN = mean)
range(raw.CC_acc$Accuracy) ## 92.1-100% 

nrow(rawq) ## 3762 data points

# get RTs,analyses only on correct trials
data_RT <- subset(data_raw, Accuracy == 1)

# trimming outliers and 3 sd by conditions
data_RT.norm <- data_RT[data_RT$RT > 100 & data_RT$RT < 20000 & !is.na(data_RT$RT), ]
(nrow(data_RT)-nrow(data_RT.norm))/nrow(data_RT)*100 ##0.25% of the data are outliers

data_RT.trim <- data_RT.norm %>% 
  group_by(Condition) %>% 
  filter(
    RT <  mean(RT) + 3*sd(RT),
    RT >  mean(RT) - 3*sd(RT)
  )
(nrow(data_RT)-nrow(data_RT.trim))/nrow(data_RT)*100 ##2.06% of the data trimmed




boxcox <- boxCox(data_RT$RT~ data_RT$Condition * data_RT$Subject)
##lamda peaks close to 0, thus use logRT to normalize the distribution of RT

qqnorm(data_RT.trim$logRT)
qqline(data_RT.trim$logRT)
## log transformation of RTs
data_RT.trim$logRT <- log(data_RT.trim$RT)
##variable naming and contrast coding
data_RT.trim$Subject   <- as.factor(data_RT.trim$Subject)
data_RT.trim$Item      <- as.factor(data_RT.trim$Item)
data_RT.trim$Condition <- as.factor(data_RT.trim$Condition)
data_RT.trim$Type      <- as.factor(data_RT.trim$Type)

## define sum contrast coding
contrasts(data_RT.trim$Condition) <- contr.sum(3)
contrasts(data_RT.trim$Condition)


contrasts(data_RT.trim$Type) <- contr.sum(3)
contrasts(data_RT.trim$Type)

----


----{"r LMM on RTs"}
#### RTs modeling using LMMs ----

## model comparison:
rt_model_null <- lmer(
  logRT ~ Condition + (1 | Subject) + (1 | Item),
  data = data_RT.trim)

rt_model_full <- lmer(
  logRT ~ Condition*Type + (1 | Subject) + (1 | Item),
  data = data_RT.trim)

anova(rt_model_null,rt_model_full)

## final model
rt_model <- lmer(logRT ~ Condition*Type + (1 + Type|Subject)
                  + (1 |Item),
                  data_RT.trim,
                  REML = FALSE)
summary(rt_model)
anova(rt_model)

## planned pairwise comparisons
emm_rt_cond <- emmeans(rt_model,~ Condition,pbkrtest.limit = 4000)
emm_rt_cond
pairs(emm_rt_cond,adjust = "tukey")

emm_rt_type <- emmeans(rt_model,~ Type,pbkrtest.limit = 4000)
emm_rt_type
pairs(emm_rt_type,adjust = "tukey")

emm_rt_cond_by_type <- emmeans(rt_model,~ Condition|Type,pbkrtest.limit = 4000)
emm_rt_cond_by_type
pairs(emm_rt_cond_by_type,adjust = "tukey")

----
  
----{"r GLMM on Accuracy"}
#### Accuracy modeling using GLMM ----
data_acc <- data_raw
##variable naming and contrast coding
data_acc$Subject   <- as.factor(data_acc$Subject)
data_acc$Item      <- as.factor(data_acc$Item)
data_acc$Condition <- as.factor(data_acc$Condition)
data_acc$Type      <- as.factor(data_acc$Type)

## define sum contrast coding
contrasts(data_acc$Condition) <- contr.sum(3)
contrasts(data_acc$Condition)

contrasts(data_acc$Type) <- contr.sum(3)
contrasts(data_acc$Type)



## final model
acc_model <- glmer(Accuracy~ Condition*Type + (1  |Subject),
                   family = binomial,
                   data = data_acc)
summary(acc_model)
Anova(acc_model)


## planned pairwise comparisons
emm_acc_cond <- emmeans(acc_model, ~Condition, pbkrtest.limit = 4000)
emm_acc_cond
pairs(emm_acc_cond,adjust = "tukey")

emm_acc_type <- emmeans(acc_model, ~Type, pbkrtest.limit = 4000)
emm_acc_type
pairs(emm_acc_type,adjust = "tukey")

emm_acc_cond_by_type <- emmeans(acc_model,~ Condition|Type,pbkrtest.limit = 4000)
emm_acc_cond_by_type
pairs(emm_acc_cond_by_type,adjust = 'tukey')


##``````#visualization
##``````  

## Helper: standard error -------------------------------------------------
se <- function(x) {
  x <- x[is.finite(x)]
  sd(x) / sqrt(length(x))
}
## ================================================================
## Recode factors for RT (data_RT.trim) and Accuracy (data_acc)
## ================================================================
data_RT.trim <- data_RT.trim %>%
  mutate(
    Condition = as.character(Condition),
    Type      = as.character(Type)
  ) %>%
  mutate(
    Condition = dplyr::recode(Condition, "CC" = "control"),
    Type      = dplyr::recode(Type, "unpragmatic" = "implausible"),
    Condition = factor(Condition, levels = c("control", "NP", "VP")),
    Type      = factor(Type,      levels = c("correct", "implausible", "ungrammatical"))
  )

data_acc <- data_acc %>%
  mutate(
    Condition = as.character(Condition),
    Type      = as.character(Type)
  ) %>%
  mutate(
    Condition = dplyr::recode(Condition, "CC" = "control"),
    Type      = dplyr::recode(Type, "unpragmatic" = "implausible"),
    Condition = factor(Condition, levels = c("control", "NP", "VP")),
    Type      = factor(Type,      levels = c("correct", "implausible", "ungrammatical"))
  )

## ================================================================
## Subject-level means
## ================================================================
subj_rt <- data_RT.trim %>%
  group_by(Subject, Condition, Type) %>%
  summarise(RT = mean(RT, na.rm = TRUE), .groups = "drop")

subj_acc <- data_acc %>%
  group_by(Subject, Condition, Type) %>%
  summarise(Accuracy = mean(Accuracy, na.rm = TRUE), .groups = "drop") %>%
  mutate(AccPerc = Accuracy * 100)

## ================================================================
## Shared colour palette
## ================================================================
cond_cols <- c(
  "control" = "#ff5ec7",
  "NP"      = "#c45eff",
  "VP"      = "#6d3fbf"
)

## ================================================================
## (a) RT panel
## ================================================================
rt_min <- min(subj_rt$RT, na.rm = TRUE)

max_rt_type <- subj_rt %>%
  group_by(Type) %>%
  summarise(y_max = max(RT, na.rm = TRUE), .groups = "drop")

sig_rt <- data.frame(
  Type   = c("correct", "correct",
             "ungrammatical", "ungrammatical"),
  group1 = c("control", "NP", "control", "control"),
  group2 = c("NP",      "VP", "NP",      "VP"),
  label  = c("***",     "***", "***",    "***")
) %>%
  left_join(max_rt_type, by = "Type") %>%
  group_by(Type) %>%
  mutate(
    y.position = y_max * (1.05 + 0.07 * (row_number() - 1))
  ) %>%
  ungroup()

rt_upper <- max(c(subj_rt$RT, sig_rt$y.position), na.rm = TRUE) * 1.02

p_rt <- ggplot(subj_rt,
               aes(x = Condition, y = RT,
                   colour = Condition, fill = Condition)) +
  geom_boxplot(width = 0.45, alpha = 0.25,
               outlier.shape = NA, colour = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.9) +
  stat_summary(fun = mean,
               fun.min = function(z) mean(z) - se(z),
               fun.max = function(z) mean(z) + se(z),
               geom = "errorbar",
               width = 0.3, size = 0.7,
               colour = "black") +
  stat_summary(fun = mean,
               geom = "point",
               shape = 18, size = 1.5,
               colour = "red") +
  facet_wrap(~ Type, nrow = 1) +
  scale_colour_manual(values = cond_cols, name = "Condition") +
  scale_fill_manual(values   = cond_cols, guide = "none") +
  coord_cartesian(ylim = c(rt_min, rt_upper)) +
  labs(x = NULL, y = "RT (ms)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    plot.margin      = margin(t = 5, r = 5, b = 5, l = 5)
  ) +
  ggpubr::stat_pvalue_manual(
    sig_rt,
    label        = "label",
    y.position   = "y.position",
    tip.length   = 0.01,
    bracket.size = 0.3,
    size         = 3.5
  )

## ================================================================
## (b) Accuracy panel — y based on data + small margin
## ================================================================
acc_min <- min(subj_acc$AccPerc, na.rm = TRUE)

max_acc_type <- subj_acc %>%
  group_by(Type) %>%
  summarise(y_max = max(AccPerc, na.rm = TRUE), .groups = "drop")

sig_acc <- data.frame(
  Type   = c("correct", "correct",
             "implausible", "implausible",
             "ungrammatical", "ungrammatical", "ungrammatical"),
  group1 = c("control", "NP",
             "control", "control",
             "control", "control", "NP"),
  group2 = c("NP",      "VP",
             "NP",      "VP",
             "NP",      "VP",     "VP"),
  label  = c("***", "**",
             "*",   "*",
             "***", "***", "*")
) %>%
  left_join(max_acc_type, by = "Type") %>%
  group_by(Type) %>%
  mutate(
    y.position = y_max * (1.05 + 0.07 * (row_number() - 1))
  ) %>%
  ungroup()

acc_upper <- max(c(subj_acc$AccPerc, sig_acc$y.position), na.rm = TRUE) * 1.02

p_acc <- ggplot(subj_acc,
                aes(x = Condition, y = AccPerc,
                    colour = Condition, fill = Condition)) +
  geom_boxplot(width = 0.45, alpha = 0.25,
               outlier.shape = NA, colour = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.9) +
  stat_summary(fun = mean,
               fun.min = function(z) mean(z) - se(z),
               fun.max = function(z) mean(z) + se(z),
               geom = "errorbar",
               width = 0.3, size = 0.7,
               colour = "black") +
  stat_summary(fun = mean,
               geom = "point",
               shape = 18, size = 1.5,
               colour = "red") +
  facet_wrap(~ Type, nrow = 1) +
  scale_colour_manual(values = cond_cols, name = "Condition") +
  scale_fill_manual(values   = cond_cols, guide = "none") +
  coord_cartesian(ylim = c(acc_min, acc_upper)) +
  labs(x = NULL, y = "Accuracy (%)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    plot.margin      = margin(t = 5, r = 5, b = 5, l = 5)
  ) +
  ggpubr::stat_pvalue_manual(
    sig_acc,
    label        = "label",
    y.position   = "y.position",
    tip.length   = 0.01,
    bracket.size = 0.3,
    size         = 3.5
  )

## ================================================================
## Combine (a) RT and (b) Accuracy
##   Legend on the right; add (a)/(b) tags
## ================================================================
behavior_plot <- (p_rt / p_acc) +
  plot_layout(heights = c(1.6, 1.6), guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

# Apply global theme tweaks (works with older patchwork)
behavior_plot <- behavior_plot +
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "bold"),
    plot.tag        = element_text(face = "bold", size = 12)
  )

behavior_plot
# ggsave("Behavior_RT_Accuracy_combined_rightLegend_dynAcc.svg",
#        behavior_plot, width = 7, height = 7)






  