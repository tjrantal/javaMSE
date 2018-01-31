package edu.deakin.timo;
import java.util.List;
import java.util.ArrayList;

/**
	Class to calculate refined composite multiscale entropy (RCME). Ported from Ihlen et al J Biomechanics 49 (2016) 1420-1428 supplementary materials Matlab code to java by Timo Rantalainen 2018 tjrantal at gmail dot com
	Licensed with the GPL.
	
	I left out parameter t as Ihlen et al used t = 1 -> has no effect. Downsample data prior to calling RCME if using t is desired

	@author Timo Rantalainen	
*/

public class RefinedCompositeMultiscaleEntropy extends BaseClass{
	private double[] rcmeResult;	/*<Array for MSE results, accessed with @link #getMSE() after construction*/
	private double[] dataIn;		/*<Array to contain the signal in*/
	private int m;					/*<Number of time lags*/
	private double r;				/*<Relative tolerance*/
	
	private int scale;				/*<Number of scales*/
	
	/**Constructor
		@param	dataIn	array for which the RCME is to be calculated
		@param	m		Number of time lags in the template vector (m = 4 used by Ihlen 2016)
		@param	r		tolerance (maximum distance, will be multiplied with dataIn SD)
		@param	scale	Number of scales (equivalent to tau in MSE)
	*/
	public RefinedCompositeMultiscaleEntropy(double[] dataIn, int m, double r, int scale){
		this.dataIn = dataIn;
		this.m = m;
		this.r = r;
		this.scale = scale;
		runAnalysis();
	}

	/**Start RCME analysis at different scales in threads*/
	private void runAnalysis(){
		double sd = std(dataIn);
		rcmeResult = new double[scale];
		List<Thread> threads = new ArrayList<Thread>();
		List<CMERunnable> cmeRunnables = new ArrayList<CMERunnable>();
		//The CME threads for various coarseness levels created here
		for (int i = 1; i<=scale;++i){
			cmeRunnables.add(new CMERunnable(dataIn,m,r*sd,i,scale));	//Thread for scale = i
			threads.add(new Thread(cmeRunnables.get(i-1)));
			threads.get(i-1).start();

		}
		
		//Join (=wait for completion) the MSE threads of various coarseness levels
		for (int i =0; i<threads.size();++i){
				try{
					((Thread) threads.get(i)).join();
					rcmeResult[i] = cmeRunnables.get(i).getCE();
				}catch(Exception er){
					rcmeResult[i] = Double.NaN;	
				}

		}
		
	}
	
	/**return RCME results array
		@return The RCME results array
	*/
	public double[] getRCME(){
		double[] ret = new double[rcmeResult.length];
		for (int i = 0;i<rcmeResult.length;++i){
				ret[i] = rcmeResult[i];
		}
		return ret;
	}

	

	
	/**A subclass to enable multi threading different coarseness level SE calculations*/
	public class CMERunnable implements Runnable{
		private double[] dataIn;
		private int m;
		private double tolerance;
		private double[] coarseGrain;
		private int tau;
		private int scale;
		private double ce;	//composite entropy
		

		/**Constructor
			@param	dataIn 	time series in (coarse-grained or otherwise)
			@param	m		the timelags to calculate entropy for
			@param	tolerance		tolerance (maximum distance)
			@param	tau		the length of the mean for the coarse-grained time series		
		*/
		public CMERunnable(double[] dataIn, int m, double tolerance, int tau,int scale){
			this.dataIn = dataIn;
			this.m = m;
			this.tolerance = tolerance;
			this.tau = tau;
			this.scale = scale;
		}	
		
		/**Implement the runnable interface*/
		public void run(){
			
			ce = compositeEntropy(dataIn,m,tolerance,tau,scale);
		}
		
		/**	Return the calculated SE
			@return SE
		*/
		public double getCE(){
			return ce;
		}
		
		
		/**Calculate composite entropy
		@param	arr 	time series in (coarse-grained or otherwise)
		@param	m		the timelags to calculate entropy for
		@param	tolerance		tolerance (maximum distance)
		@param	tau		the length of the mean for the coarse-grained time series		
		*/
		private double compositeEntropy(double[] arr, int m, double tolerance, int tau,int scale){
			double[][] n = new double[2][tau];
			//Iterate with all possible starting points of the coarse grainings
			for (int k = 0;k<tau;++k){
				//Get data with the current starting point
				double[] tempData = new double[arr.length-scale];
				for (int i = 1+k; i<arr.length-scale+k; ++i){
					tempData[i-k-1] = arr[i];
				}
				//Calculate the coarse-graining means
				double[] coarseGrain = coarseGraining(tempData, tau);
				int N = coarseGrain.length;
				double[][] count = new double[2][N-m];
				//Step 2 Iteration of template vectors. m and m+1 handled in the same loop
				for (int i = 0; i<N-m;++i){
					//Check all remaining bits of data against the current bit of data (=template)
					for (int j = i+1; j<N-m;++j){
						int kk= 0;
						//Check if data points match the template, stop incrementing at first failure
						int tempCount = 0;
						while (kk < m && Math.abs(coarseGrain[i+kk] - coarseGrain[j+kk]) <= tolerance){
							++kk;
						}
						//if kk == m all data points matched the template -> increment count, and check m+1
						if (kk == m){
							count[0][i]+=1d;	//increment m, all were under tolerance
							if (j+m < coarseGrain.length && Math.abs(coarseGrain[i+m] - coarseGrain[j+m]) <= tolerance){
								count[1][i]+=1d;	//increment m+1, all were under tolerance
							}
							
						}
						
					}
					//Normalise counts
					count[0][i]/=(double) (N-m);
					count[1][i]/=(double) (N-m);
				}
				//Calculate mean count of template matches across all template vectors j
				for (int i = 0; i<N-m;++i){
					n[0][k] +=count[0][i];
					n[1][k] +=count[1][i];
				}
				//Normalise the counts
				n[0][k]/=(double) (N-m);
				n[1][k]/=(double) (N-m);
				
			}
			//Calculate refined composite entropy
			return Math.log(mean(n[0])/mean(n[1]));
		}		
	}
	
	private double mean(double[] a){
		double b = 0;
		for (int i = 0;i<a.length;++i){
			b+=a[i];
		}
		return b/((double) a.length);		
	}
}
