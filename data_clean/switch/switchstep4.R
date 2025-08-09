library(ggplot2)
library(dplyr)
library(readxl)

# switch<-read_excel("/Users/tanlirou/Documents/yunfu/dataclean0616/switchclean/switch_afterage_afterPrano0.xlsx")
switch<-read_excel("D:/datasets/yunfu/raw_data/switch_results/switch_dependent_var.xlsx")


switch <- switch %>%
  mutate(Age_year = round(Age_year, 0))
switch_ana<-switch %>%
  group_by(Age_year) %>%
  summarise(
    switch_mean_rt = mean(Mean_RT,na.rm = T),
    sd_switch_mean_rt = sd(Mean_RT,na.rm = T),
    se_switch_mean_rt = sd(Mean_RT,na.rm = T)/sqrt(n()),
    switch_mean_acc = mean(Switch_acc,na.rm = T),
    sd_switch_mean_acc = sd(Switch_acc,na.rm = T),
    se_switch_mean_acc = sd(Switch_acc,na.rm = T)/sqrt(n()),
    switch_M_rt = mean(M_trials_RT,na.rm = T),
    sd_switch_M_rt = sd(M_trials_RT,na.rm = T),
    se_switch_M_rt = sd(M_trials_RT,na.rm = T)/sqrt(n()),
    switch_M_acc = mean(M_acc,na.rm = T),
    sd_switch_M_acc = sd(M_acc,na.rm = T),
    se_switch_M_acc = sd(M_acc,na.rm = T)/sqrt(n()),
    switch_P_rt = mean(P_trials_RT,na.rm = T),
    sd_switch_P_rt = sd(P_trials_RT,na.rm = T),
    se_switch_P_rt = sd(P_trials_RT,na.rm = T)/sqrt(n()),
    switch_P_acc = mean(P_acc,na.rm = T),
    sd_switch_P_acc = sd(P_acc,na.rm = T),
    se_switch_P_acc = sd(P_acc,na.rm = T)/sqrt(n()),
    switch_B_rt = mean(B_rt,na.rm = T),
    sd_switch_B_rt = sd(B_rt,na.rm = T),
    se_switch_B_rt = sd(B_rt,na.rm = T)/sqrt(n()),
    switch_B_acc = mean(B_acc,na.rm = T),
    sd_switch_B_acc = sd(B_acc,na.rm = T),
    se_switch_B_acc = sd(B_acc,na.rm = T)/sqrt(n()),
    switch_pure_acc = mean(Pure_trials_acc,na.rm = T),
    sd_pure_acc = sd(Pure_trials_acc,na.rm = T),
    se_pure_acc = sd_pure_acc/sqrt(n()),
    switch_pure_rt = mean(Pure_RT,na.rm = T),
    sd_pure_rt = sd(Pure_RT,na.rm = T),
    se_pure_rt = sd_pure_rt/sqrt(n()),
    switchcost_SR_acc = mean(SC_SR_acc,na.rm = T),
    sd_switchcost_SR_acc = sd(SC_SR_acc,na.rm = T),
    se_switchcost_SR_acc = sd_switchcost_SR_acc/sqrt(n()),
    switchcost_SR_rt = mean(SC_SR_rt,na.rm = T),
    sd_switchcost_SR_rt = sd(SC_SR_rt,na.rm = T),
    se_switchcost_SR_rt = sd_switchcost_SR_rt/sqrt(n()),
    switchcost_SP_acc = mean(SC_SP_acc,na.rm = T),
    sd_switchcost_SP_acc = sd(SC_SP_acc,na.rm = T),
    se_switchcost_SP_acc = sd_switchcost_SP_acc/sqrt(n()),
    switchcost_SP_rt = mean(SC_SP_rt,na.rm = T),
    sd_switchcost_SP_rt = sd(SC_SP_rt,na.rm = T),
    se_switchcost_SP_rt = sd_switchcost_SP_rt/sqrt(n()),    
    switchcost_BP_acc = mean(SC_BP_acc,na.rm = T),
    sd_switchcost_BP_acc = sd(SC_BP_acc,na.rm = T),
    se_switchcost_BP_acc = sd_switchcost_BP_acc/sqrt(n()),
    switchcost_BP_rt = mean(SC_BP_rt,na.rm = T),
    sd_switchcost_BP_rt = sd(SC_BP_rt,na.rm = T),
    se_switchcost_BP_rt = sd_switchcost_BP_rt/sqrt(n())
  )





purert <- ggplot(switch_ana, aes(x = Age_year, y = switch_pure_rt)) +
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +

  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size = 1),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    plot.title = element_text(size = 18)
  ) +
  scale_y_continuous(breaks = seq(from = 400, to = 1100, by = 100),
                     limits = c(400, 1100),
                     labels = function(y) y) +

  scale_x_continuous(breaks = seq(from = 11, to = 18, by = 1),
                     limits = c(11, 18), 
                     labels = function(x) x) +
  labs(title = "switch_pure_rt", x = "Age_year", y = "Mean reaction time (ms)")


# 显示图表  
print(purert)

pureacc<-ggplot(switch_ana, aes(x = Age_year, y = switch_pure_acc)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.05),
                     labels = function(y) y)+
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switch_pure_acc", x = "Age_year", y = "Mean accuracy") # 注意labs中的y应该大写  


# 显示图表  
print(pureacc)







pswitchrt<-ggplot(switch_ana, aes(x = Age_year, y = switch_mean_rt)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 780, by = 10),
                     labels = function(y) y)+
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switch_mean_rt", x = "Age_year", y = "Mean reaction time (ms)") # 注意labs中的y应该大写  

# 显示图表  
print(pswitchrt)


pswitchacc<-ggplot(switch_ana, aes(x = Age_year, y = switch_mean_acc)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.02),
                     labels = function(y) y)+
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switch_mean_acc", x = "Age_year", y = "Mean accuracy") # 注意labs中的y应该大写  

# 显示图表  
print(pswitchacc)

pmrt<-ggplot(switch_ana, aes(x = Age_year, y = switch_M_rt)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 800, by = 10),
                     labels = function(y) y)+
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switch_M_rt", x = "Age_year", y = "Mean reaction time (ms)") # 注意labs中的y应该大写  

# 显示图表  
print(pmrt)

pmacc<-ggplot(switch_ana, aes(x = Age_year, y = switch_M_acc)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.02),
                     labels = function(y) y)+
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switch_M_acc", x = "Age_year", y = "Mean accuracy") # 注意labs中的y应该大写  

# 显示图表  
print(pmacc)

pprt<-ggplot(switch_ana, aes(x = Age_year, y = switch_P_rt)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 780, by = 10),
                     labels = function(y) y)+
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switch_P_rt", x = "Age_year", y = "Mean reaction time (ms)") # 注意labs中的y应该大写  

# 显示图表  
print(pprt)

ppacc<-ggplot(switch_ana, aes(x = Age_year, y = switch_P_acc)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.02),
                     labels = function(y) y)+
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switch_P_acc", x = "Age_year", y = "Mean accuracy") # 注意labs中的y应该大写  

# 显示图表  
print(ppacc)

pbrt<-ggplot(switch_ana, aes(x = Age_year, y = switch_B_rt)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 800, by = 10),
                     labels = function(y) y)+
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switch_B_rt", x = "Age_year", y = "Mean reaction time (ms)") # 注意labs中的y应该大写  

# 显示图表  
print(pbrt)

pbacc<-ggplot(switch_ana, aes(x = Age_year, y = switch_B_acc)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.02),
                     labels = function(y) y)+
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switch_B_acc", x = "Age_year", y = "Mean accuracy") # 注意labs中的y应该大写  

# 显示图表  
print(pbacc)

pscrt<-ggplot(switch_ana, aes(x = Age_year, y = switchcost_SR_rt)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 780, by = 5),
                     labels = function(y) y)+ 
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switchcost_SR_rt", x = "Age_year", y = "time (ms)") # 注意labs中的y应该大写  

# 显示图表  
print(pscrt)

pscacc<-ggplot(switch_ana, aes(x = Age_year, y = switchcost_SR_acc)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  #scale_y_continuous(breaks = seq(from = -1, to = 0, by = 0.02),
  #                   labels = function(y) y)+ 
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switchcost_SR_acc", x = "Age_year", y = " ") # 注意labs中的y应该大写  

# 显示图
print(pscacc)


pscrt1<-ggplot(switch_ana, aes(x = Age_year, y = switchcost_SP_rt)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 780, by = 5),
                     labels = function(y) y)+ 
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switchcost_SP_rt", x = "Age_year", y = "time (ms)") # 注意labs中的y应该大写  

# 显示图表  
print(pscrt1)

pscacc1<-ggplot(switch_ana, aes(x = Age_year, y = switchcost_SP_acc)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  #scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.01),
   #                  labels = function(y) y)+ 
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switchcost_SP_acc", x = "Age_year", y = " ") # 注意labs中的y应该大写  

# 显示图表  
print(pscacc1)

pscrt2<-ggplot(switch_ana, aes(x = Age_year, y = switchcost_BP_rt)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = 0, to = 780, by = 5),
                     labels = function(y) y)+ 
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switchcost_BP_rt", x = "Age_year", y = "time (ms)") # 注意labs中的y应该大写  

# 显示图表  
print(pscrt2)

pscacc2<-ggplot(switch_ana, aes(x = Age_year, y = switchcost_BP_acc)) +  
  geom_point(color = "#FFB6C1", size = 2, alpha = 0.6) +
  
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE, color = "black", fill = "#FFB6C1", alpha = 0.2) +
  
  theme(panel.background = element_blank(),    # 去掉背景  
        axis.line = element_line(colour = "black", size = 1)) + # 设置轴线颜色为黑色  
  theme(axis.text.y = element_text(size = 18))+ 
  theme(axis.text.x = element_text(size = 18))+
  theme(plot.title = element_text(size = 18))+
  scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.01),
                     labels = function(y) y)+   
  scale_x_continuous(breaks = seq(min(switch_ana$Age_year), max(switch_ana$Age_year), by = 1),  
                     labels = function(x) x) +  # 设置x轴刻度  
  labs(title = "switchcost_BP_acc", x = "Age_year", y = " ") # 注意labs中的y应该大写  

# 显示图表  
print(pscacc2)
