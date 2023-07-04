### Script for the computation of the hysteresis index by G. Zuecco
### last update: September 2015

source("Hysteresis_Functions.r")

### discarding the smallest value

data <- read.csv("hysteresis_examples.csv")### Input data should be imported by a text or Excel file with two columns (i.e., Q and y) ###

Q <- data[,2]### INPUT DATA ###
y <- data[,3]### INPUT DATA ###
Q_fixed <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00) ### INPUT DATA - The user should be able to change/add/delete these values ###

### Data normalization
Qnorm <- (Q-min(Q))/(max(Q)-min(Q))
ynorm <- (y-min(y))/(max(y)-min(y))

### The function indice_Qnorm is applied to find discharge indices for observations
Qrise_h <- rep(0,length(Q_fixed))
Qrise_l <- rep(0,length(Q_fixed))
Qfall_h <- rep(0,length(Q_fixed))
Qfall_l <- rep(0,length(Q_fixed))

#last Qfall h should be 122

for (k in 1:length(Q_fixed)){
    if(isTRUE(which(Qnorm==Q_fixed[k]))){
      a <- which(Qnorm==Q_fixed[k])
    }else{
      a <- NA
    }  
#    a=ifelse(!is.null(which(Qnorm==Q_fixed[k])),which(Qnorm==Q_fixed[k]),NaN)
    emp <- is.na(a)
    #emp <- double(emp)
    if (emp==TRUE){
        Qrise_h[k]<- indice_Qnorm(Q_fixed[k],Qnorm)[[1]]
        Qrise_l[k]<- indice_Qnorm(Q_fixed[k],Qnorm)[[2]]
        Qfall_h[k]<- indice_Qnorm(Q_fixed[k],Qnorm)[[3]]
        Qfall_l[k] <- indice_Qnorm(Q_fixed[k],Qnorm)[[4]]
    } else {
        Qrise_h[k] <- a[1]
        Qrise_l[k] <- Qrise_h[k]
        Qfall_h[k] <- a[length(a)]
        Qfall_l[k] <- Qfall_h[k]
    }
}

if (any(is.na(Qfall_l)) | any(is.na(Qfall_l))){
   print('ERROR: the independent variable, x, does not reach the minimum selected value in the rising and/or the falling curve')
   check_Qfixed <- 1
} else {
   check_Qfixed <- 0
}

### Computation of slope and intercept to find the y values corresponding to Q_fixed
m_rise <- (ynorm[Qrise_h]-ynorm[Qrise_l])/(Qnorm[Qrise_h]-Qnorm[Qrise_l])# slope
m_fall <- (ynorm[Qfall_l]-ynorm[Qfall_h])/(Qnorm[Qfall_l]-Qnorm[Qfall_h])# slope
q_rise <- ((Qnorm[Qrise_h]*ynorm[Qrise_l])-(Qnorm[Qrise_l]*ynorm[Qrise_h]))/(Qnorm[Qrise_h]-Qnorm[Qrise_l])# intercept
q_fall <- ((Qnorm[Qfall_l]*ynorm[Qfall_h])-(Qnorm[Qfall_h]*ynorm[Qfall_l]))/(Qnorm[Qfall_l]-Qnorm[Qfall_h])# intercept

y_fixed_rise <- rep(0,length(Q_fixed))
y_fixed_fall <- rep(0,length(Q_fixed))



for (k in 1:length(Q_fixed)){
    if (is.na(m_rise[k])){
        y_fixed_rise[k] <- ynorm[Qrise_h[k]]
    } else {
        y_fixed_rise[k] <- m_rise[k]*Q_fixed[k]+q_rise[k]
    }
    if (is.na(m_fall[k])){
        y_fixed_fall[k] <- ynorm[Qfall_h[k]]
    } else {
        y_fixed_fall[k] <- m_fall[k]*Q_fixed[k]+q_fall[k]
    }
}

# 
# Q_fixed <- Q_fixed
# y_fixed_rise <- y_fixed_rise
# y_fixed_fall <- y_fixed_fall'


####################################################################################################################
####################################################################################################################
### Index h: computation based on integrals (or trapezoid areas)



rise_area <- area_hysteresis_index(Q_fixed,y_fixed_rise,y_fixed_fall)[[1]]
fall_area <- area_hysteresis_index(Q_fixed,y_fixed_rise,y_fixed_fall)[[2]]
diff_area <- area_hysteresis_index(Q_fixed,y_fixed_rise,y_fixed_fall)[[3]]
h <- area_hysteresis_index(Q_fixed,y_fixed_rise,y_fixed_fall)[[4]]

min_dA <- min(diff_area)
max_dA <- max(diff_area)

for (k in 1:(length(Q_fixed)-1)){
    if (is.na(rise_area[k])){  # forcing linearity (no loop) = 0
        rise_area[k] <- 0
    }
}

for (k in 1:(length(Q_fixed)-1)){
    if (is.na(fall_area[k])){# forcing linearity (no loop) = 0
        fall_area[k] <- 0
    }
}

for (k in 1:(length(Q_fixed)-1)){
    if (is.na(diff_area[k])){# forcing linearity (no loop) = 0
        diff_area[k] <- 0
    }
}

if (is.na(h)){
    h <- 0# forcing linearity (no loop) = 0
}
####################################################################################################################
####################################################################################################################
### Hysteresis class
maxQ <- max(Q)# streamflow peak
maxQ_ind=which(Q==maxQ)
min_y_rise <- min(y[1:maxQ_ind[1]])
max_y_rise <- max(y[1:maxQ_ind[1]])
change_max_y_rise <- abs(max_y_rise-y[1])
change_min_y_rise <- abs(min_y_rise-y[1])

if (change_max_y_rise>change_min_y_rise){
    if (min_dA>0 && max_dA>0){
        hyst_class <- 1
    } else {
        if (min_dA<0 && max_dA<0){
            hyst_class <- 4
        } else {
            if (min_dA<=0 && max_dA>0 && h>=0){
                hyst_class <- 2
            } else {
                if (min_dA<0 && max_dA>=0 && h<0){
                    hyst_class <- 3
                } else {
                    hyst_class <- 0# linearity
                }
            }
        }
    }
}

if (change_max_y_rise<change_min_y_rise){
    if (min_dA>0 && max_dA>0){
        hyst_class <- 5
    } else {
        if (min_dA<0 && max_dA<0){
            hyst_class <- 8
        } else {
            if (min_dA<=0 && max_dA>0 && h>=0){
                hyst_class <- 6
            } else {
                if (min_dA<0 && max_dA>=0 && h<0){
                    hyst_class <- 7
                } else {
                    hyst_class <- 0# linearity
                }
            }
        }
    }
}


if (change_max_y_rise==change_min_y_rise){
    min_y_fall <- min(y[maxQ_ind:length(y)])
    max_y_fall <- max(y[maxQ_ind:length(y)])
    change_max_y_fall <- abs(max_y_fall-y[1])
    change_min_y_fall <- abs(min_y_fall-y[1])
    if (change_max_y_fall>change_min_y_fall){
        if (min_dA>0 && max_dA>0){
            hyst_class <- 1
        } else {
            if (min_dA<0 && max_dA<0){
                hyst_class <- 4
            } else {
                if (min_dA<=0 && max_dA>0 && h>=0){
                    hyst_class <- 2
                } else {
                    if (min_dA<0 && max_dA>=0 && h<0){
                        hyst_class <- 3
                    } else {
                        hyst_class <- 0# linearity
                    }
                }
            }
        }
    } else {
        if (change_max_y_fall<change_min_y_fall){
            if (min_dA>0 && max_dA>0){
                hyst_class <- 5
            } else {
                if (min_dA<0 && max_dA<0){
                    hyst_class <- 8
                } else {
                    if (min_dA<=0 && max_dA>0 && h>=0){
                        hyst_class <- 6
                    } else {
                        if (min_dA<0 && max_dA>=0 && h<0){
                            hyst_class <- 7
                        } else {
                            hyst_class <- 0# linearity
                        }
                    }
                }
            }
        }
    }
    if (change_max_y_fall==change_min_y_fall){
        hyst_class <- 0
    }
}


### Summary
summary <- list("diff_area"=diff_area,"h"=h, "hyst_class"=hyst_class)### OUTPUT DATA

####################################################################################################################

### OUTPUT - Figure 1 ###
library(ggplot2)

ggplot()+geom_path(aes(x=Q,y=y),colour="cornflowerblue",lwd=1)+theme_light()+
  xlab("Streamflow")+ylab("y")+
  ggtitle("Hysteretic plot (input data)")

### OUTPUT - Figure 2 ###

diff_area_reshape <- c(diff_area,diff_area[length(diff_area)])

ggplot()+geom_step(aes(x=Q_fixed,y=diff_area_reshape),lwd=1)+theme_light()+
  xlim(0,1)+
  geom_hline(yintercept = 0,colour="red",lwd=1)+
  xlab("Streamflow (-)")+
  ylab(expression(paste(Delta,"A (-)")))+
  ggtitle("Difference between the integrals")


