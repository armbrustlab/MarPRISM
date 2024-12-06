library(tidyverse)
library(readxl)
library(plotrix)
library(ggpubr)

f1_xgModel_newCont <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_xgModel_noOutliers_newContaminationMetric_xg.csv")
f1_xgModel_newCont_xgRfFeatures <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_xgModel_noOutliers_newContaminationMetric_xgRf.csv")
f1_rfModel_newCont_rfFeatures <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_rfModel_noOutliers_rfFeatures_newContaminationMetric.csv")

f1_xgModel_xgFeatures_withOutliers <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_xgModel_withOutliers.csv")
f1_rfModel_rfFeatures_withOutliers <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_rfModel_withOutliers_rfFeatures.csv")

oldModel <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_oldModel.csv")


f1 <- bind_rows(f1_xgModel_newCont %>% mutate(model = "XG model, XG features, contamination removed"), 
                f1_xgModel_newCont_xgRfFeatures %>% mutate(model = "XG model, XG+RF features, contamination removed"), 
                f1_rfModel_newCont_rfFeatures %>% mutate(model = "RF model, RF features, contamination removed"), 
                
                f1_xgModel_xgFeatures_withOutliers %>% mutate(model = "XG model, XG features, with contamination"),
                f1_rfModel_rfFeatures_withOutliers %>% mutate(model = "RF model, RF features, with contamination"), 
                oldModel %>% mutate(model = "Lambert et al. model"))

f1 %>% 
  group_by(model) %>% summarize(mean = mean(F1_Scores), se = std.error(F1_Scores)) %>% 
  arrange(desc(mean))

dat <- f1 %>% 
  mutate(model = factor(model, levels = c("XG model, XG features, contamination removed", 
                                          "RF model, RF features, contamination removed", 
                                          "XG model, XG+RF features, contamination removed", 
                                          "XG model, XG features, with contamination",
                                          "RF model, RF features, with contamination", 
                                          "Lambert et al. model")))

dat %>%
  group_by(model) %>% summarize(mean = mean(F1_Scores), se = std.error(F1_Scores)) %>%
  ggplot(aes(x = model, y = mean)) + geom_bar(stat = 'identity', fill = 'white', color = 'black') + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) + 
  labs(x = "Model", y = "Mean F1 score") +
  theme(text = element_text(size=20, color = 'black')) + 
  theme_classic() + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 18, color = 'black')) +
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(axis.title.x = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 22, color = 'black')) +
  theme(legend.title = element_text(size = 26, color = 'black')) + 
  coord_cartesian(ylim=c(.6,1)) + 
  stat_compare_means(data = dat, aes(x = model, y = F1_Scores),
                     method = "anova", label.x = .5, label.y = .99) +
  stat_compare_means(data = dat, aes(x = model, y = F1_Scores),
                     label = "p.signif", method = "t.test",
                     ref.group = ".all.")   

ggsave("~/Dropbox/grad/research/crossValidation/F1_byModel_ofInterest.svg", dpi = 600, height = 13, width = 10)




f1_xgModel_newCont_class <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_xgModel_noOutliers_newContaminationMetric_xg_byClass.csv")
f1_xgModel_newCont_xgRfFeatures_class <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_xgModel_noOutliers_newContaminationMetric_xgRf_byClass.csv")
f1_rfModel_newCont_rfFeatures_class <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_rfModel_noOutliers_rfFeatures_byClass_newContaminationMetric.csv")

f1_xgModel_xgFeatures_withOutliers_class <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_xgModel_withOutliers_byClass.csv")
f1_xgModel_rfFeatures_withOutliers_class <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_rfModel_withOutliers_rfFeatures_byClass.csv")

oldModel_class <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/f1_scores_results_oldModel_byClass.csv")


f1 <- bind_rows(f1_xgModel_newCont_class %>% mutate(model = "XG model, XG features, contamination removed"), 
                f1_xgModel_newCont_xgRfFeatures_class %>% mutate(model = "XG model, XG+RF features, contamination removed"), 
                f1_rfModel_newCont_rfFeatures_class %>% mutate(model = "RF model, RF features, contamination removed"), 
                
                f1_xgModel_xgFeatures_withOutliers_class %>% mutate(model = "XG model, XG features, with contamination"),
                f1_xgModel_rfFeatures_withOutliers_class %>% mutate(model = "RF model, RF features, with contamination"), 
                oldModel_class %>% mutate(model = "Lambert et al. model"))

f1 <- f1 %>% gather(Class_0:Class_2, key = "Trophic mode", value = "F1_Scores")
f1 <- f1 %>% mutate(`Trophic mode` = parse_number(`Trophic mode`))

f1 <- f1 %>% mutate(`Trophic mode` = ifelse(`Trophic mode` == 1, "Mixotrophic", ifelse(`Trophic mode` == 0, "Heterotrophic", "Phototrophic")))

f1 %>% distinct(model) 
f1 %>% distinct(`Trophic mode`) 

dat <- f1 %>% 
  mutate(model = factor(model, levels = c("XG model, XG features, contamination removed", 
                                          "RF model, RF features, contamination removed", 
                                          "XG model, XG+RF features, contamination removed", 
                                          "XG model, XG features, with contamination",
                                          "RF model, RF features, with contamination", 
                                          "Lambert et al. model")))

dat %>%
  group_by(model, `Trophic mode`) %>% summarize(mean = mean(F1_Scores), se = std.error(F1_Scores)) %>% 
  arrange(desc(mean)) %>% filter(`Trophic mode` == "Phototrophic")

dat %>%
  group_by(model, `Trophic mode`) %>% summarize(mean = mean(F1_Scores), se = std.error(F1_Scores)) %>% 
  arrange(desc(mean)) %>% filter(`Trophic mode` == "Mixotrophic")

dat %>%
  group_by(model, `Trophic mode`) %>% summarize(mean = mean(F1_Scores), se = std.error(F1_Scores)) %>% 
  arrange(desc(mean)) %>% filter(`Trophic mode` == "Heterotrophic")


dat %>%
  group_by(model, `Trophic mode`) %>% summarize(mean = mean(F1_Scores), se = std.error(F1_Scores)) %>%
  ggplot(aes(x = model, y = mean, fill = `Trophic mode`)) + geom_bar(stat = 'identity', position = 'dodge', alpha = .65) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) + 
  labs(x = "Model", y = "Mean F1 score") +
  theme(text = element_text(size=20, color = 'black')) + 
  scale_fill_manual(values = c("Phototrophic" = "deepskyblue2", "Heterotrophic" = "red", "Mixotrophic" = "black"))+
  theme_classic() + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 18, color = 'black')) +
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(axis.title.x = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 22, color = 'black')) +
  theme(legend.title = element_text(size = 26, color = 'black')) + 
  coord_cartesian(ylim=c(.6,1)) + 
  stat_compare_means(data = dat %>%
                       filter(`Trophic mode` == "Mixotrophic"),
                     aes(x = model, y = F1_Scores),
                     method = "anova", label.x = .5, label.y = .99, color = 'black') +
  stat_compare_means(data = dat %>%
                       filter(`Trophic mode` == "Mixotrophic"),
                     aes(x = model, y = F1_Scores),
                     label = "p.signif", method = "t.test",
                     ref.group = ".all.", color = 'black')

ggsave("~/Dropbox/grad/research/crossValidation/F1_byModel_byClass_ofInterest_mixotrophy.svg", dpi = 600, height = 13, width = 15)



dat %>%
  group_by(model, `Trophic mode`) %>% summarize(mean = mean(F1_Scores), se = std.error(F1_Scores)) %>%
  ggplot(aes(x = model, y = mean, fill = `Trophic mode`)) + geom_bar(stat = 'identity', position = 'dodge', alpha = .65) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) + 
  labs(x = "Model", y = "Mean F1 score") +
  theme(text = element_text(size=20, color = 'black')) + 
  scale_fill_manual(values = c("Phototrophic" = "deepskyblue2", "Heterotrophic" = "red", "Mixotrophic" = "black"))+
  theme_classic() + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 18, color = 'black')) +
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(axis.title.x = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 22, color = 'black')) +
  theme(legend.title = element_text(size = 26, color = 'black')) + 
  coord_cartesian(ylim=c(.6,1)) + 
  stat_compare_means(data = dat %>%
                       filter(`Trophic mode` == "Phototrophic"),
                     aes(x = model, y = F1_Scores),
                     method = "anova", label.x = .5, label.y = .99, color = 'lightskyblue') +
  stat_compare_means(data = dat %>%
                       filter(`Trophic mode` == "Phototrophic"),
                     aes(x = model, y = F1_Scores),
                     label = "p.signif", method = "t.test",
                     ref.group = ".all.", color = 'lightskyblue')

ggsave("~/Dropbox/grad/research/crossValidation/F1_byModel_byClass_ofInterest_phototrophy.svg", dpi = 600, height = 13, width = 15)

dat %>%
  group_by(model, `Trophic mode`) %>% summarize(mean = mean(F1_Scores), se = std.error(F1_Scores)) %>%
  ggplot(aes(x = model, y = mean, fill = `Trophic mode`)) + geom_bar(stat = 'identity', position = 'dodge', alpha = .65) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) + 
  labs(x = "Model", y = "Mean F1 score") +
  theme(text = element_text(size=20, color = 'black')) + 
  scale_fill_manual(values = c("Phototrophic" = "deepskyblue2", "Heterotrophic" = "red", "Mixotrophic" = "black"))+
  theme_classic() + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 18, color = 'black')) +
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 26, color = 'black')) + 
  theme(axis.title.x = element_text(size = 26, color = 'black')) + 
  theme(legend.text = element_text(size = 22, color = 'black')) +
  theme(legend.title = element_text(size = 26, color = 'black')) + 
  coord_cartesian(ylim=c(.6,1)) + 
  stat_compare_means(data = dat %>%
                       filter(`Trophic mode` == "Heterotrophic"),
                     aes(x = model, y = F1_Scores),
                     method = "anova", label.x = .5, label.y = .99, color = 'red') +
  stat_compare_means(data = dat %>%
                       filter(`Trophic mode` == "Heterotrophic"),
                     aes(x = model, y = F1_Scores),
                     label = "p.signif", method = "t.test",
                     ref.group = ".all.", color = 'red')

ggsave("~/Dropbox/grad/research/crossValidation/F1_byModel_byClass_ofInterest_heterotrophy.svg", dpi = 600, height = 13, width = 15)


percData <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/k_train_size_vs_f1_score_mean_se_newContaminationMetric_xg.csv")

colnames(percData)[1:2] <- c("k", "propOfData")

head(percData)

percData %>% filter(k == 6) %>% ggplot(aes(x = propOfData, y = mean)) + geom_point() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.01,position=position_dodge(.9)) + 
  geom_smooth() + labs(x = "Percent of data used for training", y = "Mean F1 score") +
  theme(text = element_text(size=20, color = 'black'))

ggsave("~/Dropbox/grad/research/crossValidation/F1_ByPercOfData_noOutliers.png", dpi = 600, height = 6, width = 6)


percData_class <- read_csv("~/Dropbox/grad/research/trophicModePrediction/trophic-mode-working-main/k_train_size_vs_f1_score_mean_se_newContaminationMetric_xg_by_class.csv")

percData_class <- percData_class %>% gather(class_0_f1:class_2_f1, key = "var", value = "value")

percData_class <- percData_class %>% mutate(mean = parse_number(value))

percData_class <- percData_class %>% mutate(se = str_extract(value, "se.*"))
percData_class <- percData_class %>% mutate(se = parse_number(se))

percData_class <- percData_class %>% mutate(var = ifelse(var == "class_1_f1", "Mixotrophic", ifelse(var == "class_0_f1", "Heterotrophic", "Phototrophic")))

colnames(percData_class)[1:2] <- c("k", "propOfData")

percData_class %>% filter(k == 6) %>% ggplot(aes(x = propOfData, y = mean, color = var, group = var)) + geom_point() + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.018) + 
  geom_smooth(se = FALSE) + labs(x = "Percent of data used for training", y = "Mean F1 score") + 
  scale_color_manual(values = c("Phototrophic" = "deepskyblue2", "Heterotrophic" = "red", "Mixotrophic" = "black"))+ 
  labs(color = "Trophic mode prediction") + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(axis.text.x = element_text(size = 22, color = "black"))  + 
  theme(axis.text.y = element_text(size = 22, color = 'black')) + 
  theme(axis.title.y = element_text(size = 22, color = 'black')) + 
  theme(legend.text = element_text(size = 18, color = 'black')) +
  theme(legend.title = element_text(size = 18, color = 'black')) +
  theme(axis.title.x = element_text(size = 22, color = 'black')) + 
  scale_y_continuous(limits = c(.3,1.05), breaks = c(.4,.6,.8,1))

ggsave("~/Dropbox/grad/research/crossValidation/meanF1_ByK_byClass_noOutliers.png", dpi = 600, height = 6, width = 9)
ggsave("~/Dropbox/grad/research/crossValidation/meanF1_ByK_byClass_noOutliers.svg", dpi = 600, height = 6, width = 9)

