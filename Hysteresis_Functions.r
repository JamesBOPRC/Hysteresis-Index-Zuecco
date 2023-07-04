area_hysteresis_index <- function(Q_fixed,y_rising,y_falling){
# Computation of the difference between integrals on the rising and the
# falling curve. The areas of trapezoid are computed because the functions
# are discontinuous
du <- length(Q_fixed)-1
step <- rep(0,du)
areas_rising <- rep(0,du)
areas_falling <- rep(0,du)
diff_area <- rep(0,du)


for (j in 1:length(step)){
    step[j] <- Q_fixed[j+1]-Q_fixed[j]
    areas_rising[j] <- (y_rising[j+1]+y_rising[j])*step[j]/2
    areas_falling[j] <- (y_falling[j+1]+y_falling[j])*step[j]/2
    diff_area[j] <- areas_rising[j] - areas_falling[j]
}

h <- sum(diff_area)

return(list(areas_rising,areas_falling,diff_area,h))

}

indice_Qnorm <- function(Qselected,Qnorm){
  # Zuecco G. - This function finds the indices for different values of the normalized Q (=discharge)
  # The script considers the selected discharge point (for example: 0.15, 0.30, 0.45,
  # 0.60, 0.75, 0.90). The output of the function are the observed limits for the selected discharge both
  # in the rising and in the falling limb of the hydrograph. Qnorm represents the normalized discharge.
  
  # Find the streamflow peak
  maxQ <- max(Qnorm)# Streamflow peak
  
  maxQnorm=which(Qnorm==maxQ)# Index
  
  
  #########################################################################
  # Rising limb of the hydrograph
  delta_rise <- rep(0,length(Qnorm))*NaN
  for (i in 1:maxQnorm){
    delta_rise[i] <- Qnorm[i]-Qselected# First part to find the closest values to Qselected
    abs_delta_r <- abs(delta_rise)# Absolute values (vector)
  }
  
  
  # Lower than "0" = Qrise_minor
  negative <- delta_rise<0
  
  #not sure about this part!!!
  if (all(negative==FALSE)){
    Qrise_major <- NaN
    Qrise_minor <- NaN
    Qfall_major <- NaN
    Qfall_minor <- NaN
    
    return()
  }
  
  value_negative=max(delta_rise[which(negative==1)])
  Qrise_minor=max(which(delta_rise==value_negative))
  
  
  # Higher than "0" = Qrise_major
  if (delta_rise[Qrise_minor+1]>delta_rise[Qrise_minor]){
    Qrise_major <- Qrise_minor+1
  } else {
    positive <- delta_rise>0
    value_positive=min(delta_rise[which(positive==1)])# Find closest value to Qselected in the upper part of the vector
    Qrise_major=min(which(delta_rise==value_positive))# Select the closest value in time to the lower part of the vector
  }
  
  
  #########################################################################
  # Falling limb of the hydrograph
  delta_fall <- rep(0,length(Qnorm))*NaN
  
  
  
  for (j in maxQnorm:length(Qnorm)){
    delta_fall[j] <- Qnorm[j]-Qselected
    abs_delta_f <- abs(delta_fall)
  }
  
  
  # Higher than "0" = Qfall_major
  positive_f=delta_fall>=0
  
  value_positive_f=min(delta_fall[which(positive_f==TRUE)])
  Qfall_major=max(which(delta_fall==value_positive_f))
  
  if (Qfall_major==length(Qnorm)){
    Qrise_major <- NaN
    Qrise_minor <- NaN
    Qfall_major <- NaN
    Qfall_minor <- NaN
    return()
  }
  
  # Lower than "0" = Qfall_minor
  negative_f <- delta_fall<0
  
  if (all(negative_f==FALSE)){
    Qrise_major <- NaN
    Qrise_minor <- NaN
    Qfall_major <- NaN
    Qfall_minor <- NaN
    return()
  }
  
  if (delta_fall[Qfall_major+1]<delta_fall[Qfall_major]){
    Qfall_minor <- Qfall_major+1
  } else {
    negative_f <- delta_fall<0
    value_negative_f=max(delta_fall[which(negative_f==1)])
    Qfall_minor=min(which(delta_fall==value_negative_f))
  }
  
  
  ########################################################################
  return(list(Qrise_major, Qrise_minor, Qfall_major, Qfall_minor))
  
}

