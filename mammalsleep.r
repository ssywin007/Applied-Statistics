library(faraway)
data(mammalsleep)
require(MASS)
require(leaps)

summary(mammalsleep)
#4 nas in lifespan, 4 nas in gestation
rowSums(is.na(mammalsleep))

pairs(mammalsleep)
#predator ~ sleep
#danger ~ sleep
#exposure ~ sleep
#gestation ~ sleep
corrm <- cor(na.omit(mammalsleep))
#nondream ~ all x's, dream ~ x's except weight
#many x's are correlated with each other
#we can use regression single imputation to deal with the missing data
#we can use pcr to fit the model

#deal with the missing data, by regression single imputation
#if we use mammalsleep data directly, the obs whose sleep are na will not be included in this
#imputation model, so we use x here
x <- mammalsleep[ , c('body', 'brain', 'lifespan', 'gestation', 'predation', 'exposure', 'danger')]
modellifespan <- lm(lifespan ~ . , x)
modelgestation <- lm(gestation ~ . , x)
x[is.na(x[ , 'lifespan']), 'lifespan'] <- predict(modellifespan, x[is.na(x[ , 'lifespan']), ])
x[is.na(x[ , 'gestation']), 'gestation'] <- predict(modelgestation, x[is.na(x[ , 'gestation']), ])
mammalsleep[ , c('lifespan', 'gestation')] <- x[ , c('lifespan', 'gestation')]
#use sleep = dream + nondream to compute missing data
mammalsleep[is.na(mammalsleep[ , 'nondream']), "nondream"] <- mammalsleep[is.na(mammalsleep[ , 'nondream']), "sleep"] - mammalsleep[is.na(mammalsleep[ , 'nondream']), "dream"]
mammalsleep[is.na(mammalsleep[ , 'dream']), "dream"] <- mammalsleep[is.na(mammalsleep[ , 'dream']), "sleep"] - mammalsleep[is.na(mammalsleep[ , 'dream']), "nondream"]

pcmammal <- prcomp(na.omit(x), scale = TRUE)
summary(pcmammal)
round(pcmammal$rot[ , 1], 2)
round(pcmammal$rot[ , 2], 2)
round(pcmammal$rot[ , 3], 2)
#however, hard to intercept as the highly collinearity, try multiple regression
#the principal component can be treated as

#criterion bases procedures for model selection is better if we know the purpose of our model
modeln <- lm(nondream ~ ., mammalsleep[ , !names(mammalsleep) %in% c('dream', 'sleep')])
step(modeln)
#AIC suggests including gestation and danger in our model
modeln <- lm(nondream ~ gestation + danger, mammalsleep)
summary(modeln)
#diagnosis
plot(fitted(modeln),residuals(modeln),xlab="Fitted",ylab="Residuals")
abline(h=0)
#non-constant errors, do transformation for y
boxcox(modeln, plotit = T)
#we see 0 in the 95% CI, and 1 not in the CI, so do log transformation for y
logy <- log(mammalsleep[, 'nondream'])
#fitting
modeln <- lm(logy ~ ., mammalsleep[ , !names(mammalsleep) %in% c('nondream', 'dream', 'sleep')])
step(modeln)
#AIC suggests including gestation and danger in our model
modeln <- lm(logy ~ gestation + danger, mammalsleep)
summary(modeln)
#diagnosis
plot(fitted(modeln),residuals(modeln),xlab="Fitted",ylab="Residuals")
abline(h=0)
plot(fitted(modeln),sqrt(abs(residuals(modeln))), xlab="Fitted",ylab=
       expression(sqrt(hat(epsilon))))
#almost constant symmetrical variation, check residuals versus predictors plots
plot(mammalsleep[!is.na(mammalsleep[, 'nondream']) & !is.na(mammalsleep[, 'lifespan']), 'body'], residuals(modeln),xlab="body",ylab="Residuals")
abline(h=0)
#outlier
stud <- rstudent(modeln)
stud[which.max(abs(stud))]
cook <- cooks.distance(modeln)
halfnorm(cook,3,ylab="Cook’s distances")
#as we have no reason to delete the obs of asin.elephant, but it is really an outlier, 
#so try to do log transformation for body, do so to brain
logbody <- log(mammalsleep$body)
logbrain <- log(mammalsleep$brain)
modeln <- lm(logy ~ . + logbody + logbrain, mammalsleep[ , !names(mammalsleep) %in% c('body', 'brain', 'nondream', 'dream', 'sleep')])
step(modeln)
#AIC suggests including logbody, gestation and danger in our model
modeln <- lm(logy ~ logbody + gestation + danger, mammalsleep)
summary(modeln)
#diagnosis
plot(fitted(modeln),residuals(modeln),xlab="Fitted",ylab="Residuals")
abline(h=0)
plot(fitted(modeln),sqrt(abs(residuals(modeln))), xlab="Fitted",ylab=
       expression(sqrt(hat(epsilon))))
plot(log(mammalsleep[!is.na(mammalsleep[, 'nondream']) & !is.na(mammalsleep[, 'lifespan']), 'body']), residuals(modeln),xlab="logbody",ylab="Residuals")
abline(h=0)
plot(log(mammalsleep[!is.na(mammalsleep[, 'nondream']) & !is.na(mammalsleep[, 'lifespan']), 'brain']), residuals(modeln),xlab="logbrain",ylab="Residuals")
abline(h=0)
plot(mammalsleep[!is.na(mammalsleep[, 'nondream']) & !is.na(mammalsleep[, 'lifespan']), 'lifespan'], residuals(modeln),xlab="lifespan",ylab="Residuals")
abline(h=0)
plot(mammalsleep[!is.na(mammalsleep[, 'nondream']) & !is.na(mammalsleep[, 'lifespan']), 'gestation'], residuals(modeln),xlab="gestation",ylab="Residuals")
abline(h=0)
plot(mammalsleep[!is.na(mammalsleep[, 'nondream']) & !is.na(mammalsleep[, 'lifespan']), 'predation'], residuals(modeln),xlab="predation",ylab="Residuals")
abline(h=0)
plot(mammalsleep[!is.na(mammalsleep[, 'nondream']) & !is.na(mammalsleep[, 'lifespan']), 'exposure'], residuals(modeln),xlab="exposure",ylab="Residuals")
abline(h=0)
plot(mammalsleep[!is.na(mammalsleep[, 'nondream']) & !is.na(mammalsleep[, 'lifespan']), 'danger'], residuals(modeln),xlab="danger",ylab="Residuals")
abline(h=0)
#almost constant symmetrical variation, residuals versus responsor and predictors
#check normality, good
qqnorm(residuals(modeln), ylab = 'residuals')
qqline(residuals(modeln))
#check outliers and influential points, good
stud <- rstudent(modeln)
stud[which.max(abs(stud))]
cook <- cooks.distance(modeln)
halfnorm(cook,3,ylab="Cook’s distances")
#treat with categorical predictor
modelcate <- lm(logy ~ logbody * as.factor(danger) + gestation * as.factor(danger), mammalsleep)
anova(modelcate)
#the ANOVA test suggests no need to include the interaction term in our model

#fit multiple linear regression model for dream, we have done log transforamtion for body and brain
modeld <- lm(dream ~ . + logbody + logbrain, mammalsleep[ , !names(mammalsleep) %in% c('body', 'brain', 'nondream', 'sleep')])
step(modeld)
#AIC suggests including logbody, logbrain, gestation, predation and danger in our model
modeld <- lm(dream ~ logbody + logbrain + gestation + predation + danger, mammalsleep)
summary(modeld)
#diagnosis
plot(fitted(modeld),residuals(modeld),xlab="Fitted",ylab="Residuals")
abline(h=0)
plot(fitted(modeld),sqrt(abs(residuals(modeld))), xlab="Fitted",ylab=
       expression(sqrt(hat(epsilon))))
#there seems an outlier in our fit, try to find it
stud <- rstudent(modeld)
stud[which.max(abs(stud))]
cook <- cooks.distance(modeld)
halfnorm(cook,3,ylab="Cook’s distances")
qqnorm(residuals(modeld), ylab = 'residuals')
qqline(residuals(modeld))
#Asia.Eleohant is an outlier and an influential points, it is natrual that this obs 
#should be far from others, but I don't have reason to delete it.
#the residual of Echidna is the largest, whose dream value is 0. It is weird, try to kick it out of our model
logbody <- log(mammalsleep$body[mammalsleep[, 'dream'] != 0])
logbrain <- log(mammalsleep$brain[mammalsleep[, 'dream'] != 0])
modeld <- lm(dream ~ . + logbody + logbrain, mammalsleep[mammalsleep[, 'dream'] != 0 , !names(mammalsleep) %in% c('body', 'brain', 'nondream', 'sleep')])
step(modeld)
#AIC suggests including logbody, logbrain, gestation, predation and danger in our model
modeld <- lm(dream ~ logbody + logbrain + gestation + predation + danger, mammalsleep[mammalsleep[, 'dream'] != 0, ])
summary(modeld)
#diagnosis
plot(fitted(modeld),residuals(modeld),xlab="Fitted",ylab="Residuals")
abline(h=0)
plot(fitted(modeld),sqrt(abs(residuals(modeld))), xlab="Fitted",ylab=
       expression(sqrt(hat(epsilon))))
#also non-constant errors, do transformation for y
boxcox(modeld, plotit = T)
boxcox(modeld, plotit = T, lambda = seq(-0.2, 0.8, by = 0.1))
#we see 0 in the 95% CI, and 1 not in the CI, so do log transformation for y
#we should include the dropped observation in our new fit, as the value of the response of the
#observation is 0, so do log(y + 0.01) transformation for y
logy <- log(mammalsleep$dream + 0.01)
logbody <- log(mammalsleep$body)
logbrain <- log(mammalsleep$brain)
#fitting
modeld <- lm(logy ~ . + logbody + logbrain, mammalsleep[ , !names(mammalsleep) %in% c('body', 'brain', 'nondream', 'dream', 'sleep')])
step(modeld)
#AIC suggests including lifespan and danger in our model
modeld <- lm(logy ~ lifespan + danger, mammalsleep)
summary(modeld)
#diagnosis
plot(fitted(modeld),residuals(modeld),xlab="Fitted",ylab="Residuals")
abline(h=0)
plot(fitted(modeld),sqrt(abs(residuals(modeld))), xlab="Fitted",ylab=
       expression(sqrt(hat(epsilon))))
stud <- rstudent(modeld)
stud[which.max(abs(stud))]
cook <- cooks.distance(modeld)
halfnorm(cook,3,ylab="Cook’s distances")
qqnorm(residuals(modeld), ylab = 'residuals')
qqline(residuals(modeld))
#outlier Echidna
logbody <- log(mammalsleep$body[mammalsleep[, 'dream'] != 0])
logbrain <- log(mammalsleep$brain[mammalsleep[, 'dream'] != 0])
logy <- log(mammalsleep$dream[mammalsleep[, 'dream'] != 0])
#fitting
modeld <- lm(logy ~ . + logbody + logbrain, mammalsleep[mammalsleep[, 'dream'] != 0 , !names(mammalsleep) %in% c('body', 'brain', 'nondream', 'dream', 'sleep')])
step(modeld)
#AIC suggests including logbody, logbrain, gestation, predation and danger in our model
modeld <- lm(logy ~ logbody + logbrain + gestation + predation + danger, mammalsleep[mammalsleep[, 'dream'] != 0, ])
summary(modeld)
#suffer from outliers
plot(fitted(modeld),residuals(modeld),xlab="Fitted",ylab="Residuals")
abline(h=0)
plot(fitted(modeld),sqrt(abs(residuals(modeld))), xlab="Fitted",ylab=
       expression(sqrt(hat(epsilon))))
#Asian.Elephant is an outlier and influential obs, try to delete it
stud <- rstudent(modeld)
stud[which.max(abs(stud))]
cook <- cooks.distance(modeld)
halfnorm(cook,3,ylab="Cook’s distances")
logbody <- log(mammalsleep$body[mammalsleep[, 'dream'] != 0 & mammalsleep[, 'body'] < 2000])
logbrain <- log(mammalsleep$brain[mammalsleep[, 'dream'] != 0 &  mammalsleep[, 'body'] < 2000])
logy <- log(mammalsleep$dream[mammalsleep[, 'dream'] != 0 &  mammalsleep[, 'body'] < 2000])
#fitting
modeld <- lm(logy ~ . + logbody + logbrain, mammalsleep[mammalsleep[, 'dream'] != 0 & mammalsleep[, 'body'] < 2000 , !names(mammalsleep) %in% c('body', 'brain', 'nondream', 'dream', 'sleep')])
step(modeld)
#AIC suggests including logbody, logbrain, gestation, predation and danger in our model
modeld <- lm(logy ~ logbody + logbrain + gestation + predation + danger, mammalsleep[mammalsleep[, 'dream'] != 0 & mammalsleep[, 'body'] < 2000, ])
summary(modeld)
#coefs of logbody and logbrain almost equal but with opposed signals, which make nonsense
#diagnosis
plot(fitted(modeld),residuals(modeld),xlab="Fitted",ylab="Residuals")
abline(h=0)
plot(fitted(modeld),sqrt(abs(residuals(modeld))), xlab="Fitted",ylab=
       expression(sqrt(hat(epsilon))))
#almost constant symmetrical variation, check residuals versus predictors plots
plot(log(mammalsleep[!is.na(mammalsleep[, 'dream']) & !is.na(mammalsleep[, 'lifespan']) & mammalsleep[, 'dream'] != 0 & mammalsleep[, 'body'] < 2000, 'body']), residuals(modeld),xlab="logbody",ylab="Residuals")
abline(h=0)
plot(log(mammalsleep[!is.na(mammalsleep[, 'dream']) & !is.na(mammalsleep[, 'lifespan']) & mammalsleep[, 'dream'] != 0 & mammalsleep[, 'body'] < 2000, 'brain']), residuals(modeld),xlab="logbrain",ylab="Residuals")
abline(h=0)
plot(mammalsleep[!is.na(mammalsleep[, 'dream']) & !is.na(mammalsleep[, 'lifespan']) & mammalsleep[, 'dream'] != 0 & mammalsleep[, 'body'] < 2000, 'lifespan'], residuals(modeld),xlab="lifespan",ylab="Residuals")
abline(h=0)
plot(mammalsleep[!is.na(mammalsleep[, 'dream']) & !is.na(mammalsleep[, 'lifespan']) & mammalsleep[, 'dream'] != 0 & mammalsleep[, 'body'] < 2000, 'gestation'], residuals(modeld),xlab="gestation",ylab="Residuals")
abline(h=0)
plot(mammalsleep[!is.na(mammalsleep[, 'dream']) & !is.na(mammalsleep[, 'lifespan']) & mammalsleep[, 'dream'] != 0 & mammalsleep[, 'body'] < 2000, 'predation'], residuals(modeld),xlab="predation",ylab="Residuals")
abline(h=0)
plot(mammalsleep[!is.na(mammalsleep[, 'dream']) & !is.na(mammalsleep[, 'lifespan']) & mammalsleep[, 'dream'] != 0 & mammalsleep[, 'body'] < 2000, 'exposure'], residuals(modeld),xlab="exposure",ylab="Residuals")
abline(h=0)
plot(mammalsleep[!is.na(mammalsleep[, 'dream']) & !is.na(mammalsleep[, 'lifespan']) &  mammalsleep[, 'dream'] != 0 & mammalsleep[, 'body'] < 2000, 'danger'], residuals(modeld),xlab="danger",ylab="Residuals")
abline(h=0)
#almost constant symmetrical variation, residuals versus responsor and predictors
#check normality, good
qqnorm(residuals(modeld), ylab = 'residuals')
qqline(residuals(modeld))
#check outliers and influential points, good
stud <- rstudent(modeld)
stud[which.max(abs(stud))]
cook <- cooks.distance(modeld)
halfnorm(cook,3,ylab="Cook’s distances")
#predict
data <- mammalsleep["Kangaroo", ]
data$body <- log(data$body)
data$brain <- log(data$brain)
colnames(data)[1] <- 'logbody'
colnames(data)[2] <- 'logbrain'
lognk <- predict(modeln, data, interval = 'prediction')
nk <- exp(lognk)
logdk <- predict(modeld, data, interval = 'prediction')
dk <- exp(logdk)


require(leaps)
b <- regsubsets(logy ~ . + logbody + logbrain, mammalsleep[mammalsleep[, 'dream'] != 0 , !names(mammalsleep) %in% c('body', 'brain', 'nondream', 'dream', 'sleep')])
rs <- summary(b)
rs$which
AIC <- 48*log(rs$rss/48) + (2:8) *2
AIC
plot(AIC ~ I(1:7), ylab="AIC", xlab="Number of Predictors")
