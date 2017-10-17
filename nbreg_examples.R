dispersion = function(model){
  #Calculate the dispersion statistic
  return(sum(residuals(model, type = "pearson")^2)/model$df.residual)
}

###1 Affairs
affairs = read.csv('affairs.csv')
str(affairs)
affairs$relig = as.factor(affairs$relig)
affairs$ratemarr = as.factor(affairs$ratemarr)


#Poisson Model Selection
m1 = glm(naffairs ~ male + age + educ + occup + affair + kids + avgmarr + hapavg + vryhap + notrel + slghtrel + smerel + vryrel + yrsmarr, data=affairs, family=poisson)
summary(m1)
m2 = glm(naffairs ~ kids + avgmarr + hapavg + vryhap + notrel + slghtrel + smerel + vryrel + yrsmarr, data=affairs, family=poisson)
summary(m2)

dispersion(m2)

#Legrange Multiplier Test for overdisp
mu = predict(m2, type='response')
obs = nrow(affairs)
mmu = mean(mu)
n.ybar = obs*mmu
musq = mu^2
mu2 = mean(musq)*obs
chi.val = (mu2 - n.ybar)^2 / (2*mu2)
chi.val
pchisq(chi.val,1,lower.tail=F)

#NB2
m1nb = glm.nb(naffairs ~ kids + avgmarr + hapavg + vryhap + notrel + slghtrel + smerel + vryrel + yrsmarr, data=affairs)
summary(m1nb)

#boundary LRT
LR = -2*logLik(m2) + m1nb$twologlik
pchisq(LR,1,lower.tail=F)/2

#reduced NB2 model
m2nb = glm.nb(naffairs ~  avgmarr + hapavg + vryhap + smerel + vryrel + yrsmarr, data=affairs)
summary(m2nb)



###2 Heart Procedures
heart = read.csv('azpro.csv')

#graphical visualization of los by procedure
par(mfrow=c(1,2))
hist(heart$los[heart$procedure == 0], main = 'PTCA',
     ylim = c(0,1250),
     xlim = c(0,50), xlab = 'LOS')
hist(heart$los[heart$procedure == 1], main = 'CABG',
     ylim = c(0,1250),
     xlim = c(0,50), xlab = 'LOS')

#poisson model
m1 = glm(los ~ procedure + sex + admit + age75, data=heart, family=poisson)
summary(m1)

dispersion(m1)

#score test for overdisp
mu = as.numeric(predict(m1,type='response'))
z = ((heart$los - mu)^2 - heart$los) / (mu * sqrt(2))
zscore = lm(z ~ 1)
summary(zscore)

#NB2
m1nb = glm.nb(los ~ procedure + sex + admit + age75, data=heart)
summary(m1nb)

#boundary LRT
LR = -2*logLik(m1) + m1nb$twologlik
pchisq(LR,1,lower.tail=F)/2



###3 Titanic
titanic = read.csv('titanic.csv')

#convert to grouped counts
class1=NULL;class2=NULL;class3=NULL
for(i in 1:dim(titanic)[1]){
  
  if(titanic$age[i] == 1){
    titanic$age[i] = 'Adult'
  }else{
    titanic$age[i] = 'Child'
  }
  
  if(titanic$sex[i] == 1){
    titanic$sex[i] = 'Male'
  }else{
    titanic$sex[i] = 'Female'
  }
  
  if(titanic$class[i] == 1){
    class1[i] = 1
  }else if(titanic$class[i] != 1){
    class1[i] = 0
  }
  
  if(titanic$class[i] == 2){
    class2[i] = 1
  }else if(titanic$class[i] != 2){
    class2[i] = 0
  }
  
  if(titanic$class[i] == 3){
    class3[i] = 1
  }else if(titanic$class[i] != 3){
    class3[i] = 0
  }
}

titanic = cbind(titanic,class1,class2,class3)

survive=NULL
cases=NULL

for(i in c('Adult','Child')){
  for(j in c('Female','Male')){
    for(k in 1:3){
      
      survive = c(survive,sum(titanic$survived[titanic$age == i & titanic$sex == j & titanic$class == k]))
      cases = c(cases,length(titanic$survived[titanic$age == i & titanic$sex == j & titanic$class == k]))
      
    }
  }
}

#create df Titanic grouped data set
age = c(rep('adults',6),rep('child',6))
sex = c(rep('F',3),rep('M',3),rep('F',3),rep('M',3))
class = c(rep(seq(1,3,1),4))

titanic.count = data.frame(survive,cases,age,sex,class)

#poisson model
m1 = glm(survive ~ age + sex + factor(class) + offset(log(cases)), data=titanic.count, family=poisson)
summary(m1)

dispersion(m1)

#NB2
m1nb = glm.nb(survive ~ age + sex + factor(class) + offset(log(cases)), data=titanic.count)
summary(m1nb)

#boundary LRT
LR = -2*logLik(m1) + m1nb$twologlik
pchisq(LR,1,lower.tail=F)/2
