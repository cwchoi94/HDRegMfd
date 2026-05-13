### simulation results

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggtext)
library(patchwork)


path0 = './QM_code/simulation/sim_result/'

# basic parameters
c.h.list = 0.1*seq(1,19,length.out=7)

tau.list = c(0.1,0.5,0.9)
beta.norm.list = c(0.5,1.0)
error.type.list = c('normal-homoscedastic','normal-heteroscedastic','t-homoscedastic','t-heteroscedastic')
penalty.list = c('LASSO','SCAD','MCP','ORACLE')

kernel = 'Gaussian'
cv.type = 'BIC'

model = 'QM'
n = 200
p = 400
s = 3


path = paste0(path0,'sensitivity_',kernel,'_',cv.type,'_',n,'_',p,'_',s,'/result_each_QM/')


# define score
score = 'RMSE'
score = 'MAD'
score.list = c('TP','FP','MSE','MAD','runtime.opt.second')
result.cols = c('tau','penalty','c.h','iter',score.list)

# save results as csv files
loop1 = expand.grid(error.type=error.type.list)
loop2 = expand.grid(tau=tau.list,c.h=c.h.list,penalty=penalty.list)

for (idx in 1:nrow(loop1)){
  error.type = loop1[idx,'error.type']
  
  print(error.type)
  
  # score results
  result.list.all = list()
  for (j in 1:nrow(loop2)){
    tau = loop2[j,'tau']
    penalty = loop2[j,'penalty']
    c.h = loop2[j,'c.h']
    
    file.name = paste0(path,'result_',error.type,'_',penalty,'_',tau,'_',c.h,'.csv')
    result.each = read.csv(file.name)[-1,-1]
    result.each['tau'] = tau
    result.each['penalty'] = penalty
    result.each['c.h'] = c.h
    rownames(result.each) = NULL
    
    result.list.all[[j]] = result.each
  }
  
  result.list.all = do.call(rbind,result.list.all)
  result.list.all = result.list.all[result.cols]
  rownames(result.list.all) = NULL
  
  
  # make long data.frame
  ch_levels <- sprintf('%.1f',sort((unique(result.list.all$c.h))))
  penalty_levels = penalty.list
  
  df_long <- result.list.all %>%
    pivot_longer(
      cols = c(TP, FP, MSE, MAD, runtime.opt.second),
      names_to = "measure",
      values_to = "value"
    ) %>%
    mutate(measure = case_when(measure=='runtime.opt.second' ~ 'Time (sec)',TRUE ~ as.character(measure))) %>%  
    mutate(tau_label = paste0("bold(tau == ", tau,')')) %>% 
    mutate(measure = factor(measure, levels = c("TP", "FP", "MSE", "MAD", "Time (sec)"))) %>% 
    mutate(c.h = sprintf("%.1f", as.numeric(as.character(c.h)))) %>%
    mutate(c.h = factor(c.h, levels = ch_levels)) %>% 
    mutate(penalty = factor(penalty, levels = penalty.list))
  
  # make dummy data (to set y-range)
  dummy_data <- data.frame(
    measure = factor(c("TP", "FP", "MSE", "MAD", "Time (sec)"), 
                     levels = c("TP", "FP", "MSE", "MAD", "Time (sec)")),
    value = c(0, 3, 
              0, 0, 
              0, 0,
              0, 0,
              0, 0), # TP range = c(0,3), lower bdds of the remainders = 0
    tau_label = unique(df_long$tau_label)[1],
    c.h = factor(ch_levels[1], levels = ch_levels),
    penalty = factor(penalty_levels[1], levels=penalty.list)
  )
  
  # draw ggplots
  plot = ggplot(df_long, aes(x = factor(c.h), y = value, color = penalty, group = penalty)) +
    geom_blank(data = dummy_data, aes(y = value)) +
    stat_summary(fun = mean, geom = "line", size = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 2) +
    # IQR bar
    stat_summary(fun.data = function(x) {
      data.frame(y = median(x), 
                 ymin = quantile(x, 0.25), 
                 ymax = quantile(x, 0.75))
    }, geom = "errorbar", width = 0.2, alpha = 0.7) +
    # 5(measure) x 3(tau) subplots
    facet_grid(measure ~ tau_label, scales = "free_y", labeller = label_parsed) +
    scale_x_discrete(breaks = unique(df_long$c.h)) +
    theme_minimal(base_size = 16) +
    theme(
      panel.spacing.x = unit(1.5, "lines"), 
      panel.spacing.y = unit(2.0, "lines"), 
      strip.text = element_text(size = 13, face='bold'),
      plot.title = element_markdown(hjust = 0.5, size = 18),
      legend.position = "bottom",
      legend.title = element_blank(), 
      legend.text = element_text(size = 14),
      axis.title.x = element_text(size = 15,face='bold',margin=margin(t=15)),
      axis.title.y = element_text(size = 15,face='bold'),
      axis.text = element_text(size = 12),
      panel.border = element_rect(color = "grey80", fill = NA, size = 1)
    ) +
    labs(x = expression(c[h]),y = '') +
    guides(color = guide_legend(byrow = TRUE))
  
  
  ggsave(paste0(path0,'sensitivity_',error.type,'.pdf'),plot,width=10,height=8)
}





