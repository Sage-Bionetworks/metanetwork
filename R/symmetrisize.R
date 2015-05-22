symmetrisize <- function(x){
  require(dplyr)
  x/2 + t(x)/2 %>% return
}