.affyILM_scale<-function(data.all,method)
{
	methods<-c("linear","linear.stat")
	if(!sum(methods==method))
		{
			message("Requested method not known, using method linear")
			method<-"linear"
		}
	if(method=="linear.stat")
		{
		medians<-apply(data.all,2,median,na.rm=T)
		mads<-apply(data.all,2,mad,na.rm=T)
		medians.med<-median(medians)
		#we use the value which is closest to the median value
		diffs<-abs(medians-medians.med)
		med.ref<-which(diffs==min(diffs))
		median.all<-medians[med.ref]
		mad.all<-mads[med.ref]
		#precaution if two values are both as close to the median
		if(length(median.all)>1)
			{
			median.all<-median(median.all);
			mad.all<-median(mad.all);
			}
		data.scaled<-data.all;
		for(i in 1:ncol(data.all))
			{
			data.scaled[,i]<-(((data.all[,i]-medians[i])/(1.4826*mads[i]))*1.4826*mad.all)+median.all;
			}
		}
	if(method=="linear")
		{
		data.medians<-apply(data.all,2,median,na.rm=T);
		#computing pairwise slopes
		slopes<-c();
		for(i in 1:ncol(data.all))
			{
			index<-c(1:ncol(data.all));
			slopes.temp<-data.medians[index]/data.medians[i];
			slopes<-cbind(slopes,slopes.temp);
			}
		slopes.medians<-apply(slopes,2,median,na.rm=T);
		data.scaled<-data.all;
		for(i in 1:ncol(data.all))
			{
			data.scaled[,i]<-data.all[,i]*slopes.medians[i];
			}
		#check: compute slopes after modification
#		data2.medians<-apply(data.scaled,2,median,na.rm=T);
#		slopes2<-c();
#		for(i in 1:ncol(data.all))
#			{
#			index<-c(1:ncol(data.all));
#			slopes2.temp<-data2.medians[index]/data2.medians[i];
#			slopes2<-cbind(slopes2,slopes2.temp);
#			}
#		slopes2.medians<-apply(slopes2,2,median,na.rm=T);
#		message("posterior pairwise slopes");
#		print(slopes2)
#		message("posterior medians of slopes");
#		print(slopes2.medians)
		}
	return(data.scaled);
}



#check ok, all slopes = 1.
