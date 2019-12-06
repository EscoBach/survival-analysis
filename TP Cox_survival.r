#http://www.sthda.com/english/wiki/cox-proportional-hazards-model

library("survival")
#library("survminer")

data("lung")
head(lung)

res.cox <- coxph(Surv(time, status) ~ sex, data = lung)
summary(res.cox)

res.cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data =  lung)
summary(res.cox)

#ggsurvplot(survfit(res.cox), color = "#2E9FDF",ggtheme = theme_minimal())

sex_df <- with(lung, data.frame(sex = c(1, 2), age = rep(mean(age, na.rm = TRUE), 2),ph.ecog = c(1, 1)))
sex_df

fit <- survfit(res.cox, newdata = sex_df)
#ggsurvplot(fit, conf.int = TRUE, legend.labs=c("Sex=1", "Sex=2"),ggtheme = theme_minimal())

res.c <- cox.zph(res.cox)

par(mfrow=c(1,3))

res.c
plot(res.c)