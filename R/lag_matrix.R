# convolve the DL function with exposure to give outcome
lag_matrix	<-	function(rain, p, start.at.zero = T){
  slider	<-	function(n)		rev(rain[n:(n+p-!start.at.zero)])
  nums		<-	as.list(1:(length(rain) - p))
  matrix(unlist(lapply(nums, slider)), nrow = (length(rain) - p), ncol = p  + start.at.zero, byrow = TRUE)
}