package edu.deakin.timo;

/*
	Helper class to store code shared by the different entropy estimation methods
*/

public abstract class BaseClass{
	/**Standard deviation
		This method of calculating sd explained in e.g. http://www.johndcook.com/blog/2008/09/26/comparing-three-methods-of-computing-standard-deviation/
		
		@param	arr	array for which the sd is to be calculated
		@return standard deviation
	*/
	protected double std(double[] arr){
		double sum = 0,sum2 = 0;
		for (int i = 0; i<arr.length;++i){
			sum+=arr[i];
			sum2+=arr[i]*arr[i];
		}
		return Math.sqrt((sum2-sum*sum/((double) arr.length))/(((double) arr.length)-1d));		
	}
	
	/**Create coarse-grained time series
		@param arr			original time series in
		@param scaleFactor	the length of mean for the coarse-grained time series
		@return				Coarse-grained time series
	*/
	protected double[] coarseGraining(double[] arr, int scaleFactor){
		double[] y = new double[arr.length/scaleFactor];
		for (int i = 0; i<arr.length/scaleFactor;++i){
			y[i] = 0;
			for (int k=0;k<scaleFactor;++k){
				y[i] += arr[i*scaleFactor+k];
			}
			y[i] /= (double) scaleFactor;
		}
		return y;
	}
}