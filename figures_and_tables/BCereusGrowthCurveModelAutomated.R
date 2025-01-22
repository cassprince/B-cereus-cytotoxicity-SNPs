#BEFORE USING THIS CODE:
##set working directory and change file name
##change output .csv and .pdf file names and location (location is optional)


library(growthcurver)

##loading the data into R, data for three isolates are contained in this set
setwd("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab")
BCTest <- read.csv("32 dC Growth Curve Avgs 7_8_21.csv", header= TRUE, stringsAsFactors = FALSE)

##counters for column number and array index, respectively
count<-2
counter<-1
TimeCol<-as.numeric(na.omit(BCTest[,1]))
TimeColLength<-length(TimeCol)
##created empty array for early stationary phase time values
EStatTime<-0
EarlyStatTime<-array(c(EStatTime))


##while loop to run the data of each isolate and add the time values to the array "EarlyStatTime"
pdf(file = "C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\32 dC Growth Models 7_8_21.pdf")
while (count<(ncol(BCTest)+1))
{
ColMax<-max(BCTest[,count], na.rm = TRUE)
ColMaxInd<-match(ColMax, BCTest[,count])
BCTestCol<-as.numeric(na.omit(BCTest[,count]))
ColShort<-as.numeric(na.omit(BCTestCol))
Threshhold<-ColMax-0.2

##omit any rows without data based on omitted data from death phase
NAOmitted<-na.omit(BCTestCol)
ShortLength<-length(NAOmitted)
TimeLength<-length(TimeCol)
ColVal<-0

StopCon<- FALSE

##this inner loop crops the column so that die off does not negatively affect the growth curve fit
##it will run when the maximum absorbance value is not the last data point in the set, when the stopping condition is false, and when the index being compared to the threshhold does not exceed the length of the column
##upon entering the inner while loop, "ColMaxInd" acts as a counter
##the code will scan down the column, index by index, comparing each absorbance value after the maximum value to a threshhold
##if the index being compared to the threshhold DOES exceed the length of the column, the while loop will end
##the loop will also end when the stopping condition (the Boolean variable "StopCon" TRUE), therefore no longer meeting the requirements for the while statement

while(ColMaxInd < TimeLength && StopCon==FALSE && ColMaxInd<ShortLength+1){
  
ColVal<-BCTestCol[as.integer(ColMaxInd)]

if(ColVal < Threshhold){
ColShort<-BCTestCol[1:ColMaxInd-1]
ShortLength<-length(ColShort)

StopCon<-TRUE

}

ColMaxInd<-ColMaxInd+1
}

##determine number of rows in column being read
ShortTime <-TimeCol[1:nrow(as.matrix(ColShort))]

##fitting the metrics and curve to isolate 
gc_fit <- SummarizeGrowth(ShortTime, ColShort)

##defined "gc_function" as the logistic curve with metrics determined by growthcurver
##"t" is used as the variable for time
gc_function<-function(t){
gc_var<-(gc_fit$vals$k/(1 +((gc_fit$vals$k-gc_fit$vals$n0)/gc_fit$vals$n0)*exp(-(gc_fit$vals$r)*t)))
return(gc_var)
}

##solving for the slope of the tangent line at the midpoint of the logistic curve
slope<-(((gc_fit$vals$k*0.5)-(gc_function((gc_fit$vals$t_mid)+0.01)))/((gc_fit$vals$t_mid)-((gc_fit$vals$t_mid)+0.01)))
int_b<-(gc_fit$vals$k*0.5)-(slope*gc_fit$vals$t_mid)

##"int_lines" is the intersection of the tangent lines of stationary phase and log phase midpoint 
##int_lines+2 is the value used by Fermanian et. al. to define the entry into early stationary phase
int_lines<-((gc_fit$vals$k-int_b)/slope)
EStatTime<-int_lines+2

gc_fit <- SummarizeGrowth(ShortTime, ColShort)
plot(gc_fit, xlab="Time (in hours)", ylab="OD", main=(colnames(BCTest[count])) )

EarlyStatTime[counter]<-EStatTime

counter<-counter+1
count<-count+1

}

dev.off()

##Make a data frame including names of isolates and corresponding times to early stationary phase.
EarlyStatTimeDF<-data.frame(EarlyStatTime, row.names=(colnames(BCTest[,2:(ncol(BCTest))])))

##I like to print the values just to make sure the code worked and the values are reasonable 
print(EarlyStatTimeDF)
write.csv(EarlyStatTimeDF,'C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\32 dC EStatTimes 7_8_21.csv')

#other useful values
##str(gc_fit$vals)
#plot an isolate to see the fit
##gc_fit <- SummarizeGrowth(TimeCol, BCTest[,16])
##plot(gc_fit, xlab="Time (in hours)", ylab="OD")