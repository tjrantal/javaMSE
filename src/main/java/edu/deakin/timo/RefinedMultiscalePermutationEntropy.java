package edu.deakin.timo;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

/**
	Class to calculate refined multiscale permutation entropy (RMPE). Ported from Ihlen et al J Biomechanics 49 (2016) 1420-1428 supplementary materials Matlab code to java by Timo Rantalainen 2018 tjrantal at gmail dot com
	Licensed with the GPL.
	
	I left out parameter t as Ihlen et al used t = 1 -> has no effect. Downsample data prior to calling RCME if using t is desired

	@author Timo Rantalainen	
*/

public class RefinedMultiscalePermutationEntropy extends BaseClass{
	private double[] rmpeResult;	/*<Array for MSE results, accessed with @link #getMSE() after construction*/
	private double[] dataIn;		/*<Array to contain the signal in*/
	private int m;					/*<Number of time lags*/
	private int scale;				/*<Number of scales*/
	
	/**Constructor
		@param	dataIn	array for which the RCME is to be calculated
		@param	m		Number of time lags in the template vector (m = 4 used by Ihlen 2016)
		@param	r		tolerance (maximum distance, will be multiplied with dataIn SD)
		@param	scale	Number of scales (equivalent to tau in MSE)
	*/
	public RefinedMultiscalePermutationEntropy(double[] dataIn, int m, int scale){
		this.dataIn = dataIn;
		this.m = m;
		this.scale = scale;
		runAnalysis();
	}

	/**Start RMPE analysis at different scales in threads*/
	private void runAnalysis(){
		int[][] permlist = permutations(new int[]{0,1,2,3});
		rmpeResult = new double[scale];
		List<Thread> threads = new ArrayList<Thread>();
		List<PERunnable> peRunnables = new ArrayList<PERunnable>();
		//The PE threads for various coarseness levels created here
		for (int i = 1; i<=scale;++i){
			peRunnables.add(new PERunnable(dataIn,m,i,scale,permlist));	//Thread for scale = i
			threads.add(new Thread(peRunnables.get(i-1)));
			threads.get(i-1).start();
		}
		
		//Join (=wait for completion) the PE threads of various coarseness levels
		for (int i =0; i<threads.size();++i){
				try{
					((Thread) threads.get(i)).join();
					rmpeResult[i] = peRunnables.get(i).getPE();
				}catch(Exception er){
					rmpeResult[i] = Double.NaN;	
				}

		}
		
	}
	
	/**return RCME results array
		@return The RCME results array
	*/
	public double[] getRMPE(){
		double[] ret = new double[rmpeResult.length];
		for (int i = 0;i<rmpeResult.length;++i){
				ret[i] = rmpeResult[i];
		}
		return ret;
	}

	

	
	/**A subclass to enable multi threading different coarseness level SE calculations*/
	public class PERunnable implements Runnable{
		private double[] dataIn;
		private int m;
		private double[] coarseGrain;
		private int[][] permlist;
		private int tau;
		private int scale;
		private double pe;	//permutation entropy
		

		/**Constructor
			@param	dataIn 	time series in (coarse-grained or otherwise)
			@param	m		the timelags to calculate entropy for
			@param	tau		the length of the mean for the coarse-grained time series		
			@param	scale	the number of different coarse grains considered altogether
		*/
		public PERunnable(double[] dataIn, int m, int tau,int scale, int[][] permlist){
			this.dataIn = dataIn;
			this.m = m;
			this.tau = tau;
			this.scale = scale;
			this.permlist = permlist;
		}	
		
		/**Implement the runnable interface*/
		public void run(){
			
			pe = permutationEntropy();
		}
		
		/**	Return the calculated SE
			@return SE
		*/
		public double getPE(){
			return pe;
		}
		
		
		/**
			Calculate permutation entropy
		*/
		private double permutationEntropy(){
			double[][] n = new double[permlist.length][tau];
			//Iterate with all possible starting points of the coarse grainings
			for (int k = 0;k<tau;++k){
				//Get data with the current starting point
				double[] tempData = new double[dataIn.length-scale];
				for (int i = k; i<dataIn.length-scale+k; ++i){  	//This differs from Ihlen implementation (Ihlen has i = 1+k)
					tempData[i-k] = dataIn[i];	//This differs from Ihlen implementation (Ihlen has tempData[i-k-1])
				}
				//Calculate the coarse-graining means
				double[] coarseGrain = coarseGraining(tempData, tau);
				int N = coarseGrain.length;
				//Up to this point all but identical with RCME
				//Ordering of values is compared to the permutations (template orderings)

				Sortable[] currentEpoch = new Sortable[m];
				int[] currentOrder = new int[m];
				for (int i = 0; i<N-m+1;++i){
					//Get the current epoch of data to consider
					for (int j = 0;j<m;++j){
						currentEpoch[j] = new Sortable(coarseGrain[i+j],j);
					}
					//Get the ordering of the elements in the current epoch
					Arrays.sort(currentEpoch);
					for (int j = 0;j<m;++j){
						currentOrder[j] = currentEpoch[j].b;	//Get the ordering of the epoch
					}
					
					//Check whether this matches any template
					
					for (int p = 0; p<permlist.length;++p){
						int ma = 0;
						//Check whether the template matches the current epoch order
						while (ma<m && permlist[p][ma] == currentOrder[ma]){
							++ma;
						}
						//If ma == m then all indices have matched
						if (ma == m){
							n[p][k]++;	//Matching template, increment count
							break;	//Found match, no other permutation can match, escape the loop
						}
					}
					
				}
			}
			
			//Calculate mean counts across all k iterations, and leave out templates with 0 matches
			ArrayList<Double> meanCounts = new ArrayList<Double>();
			double sum = 0;
			for (int p = 0;p<permlist.length;++p){
				double meanVal = mean(n[p]);
				if (meanVal > 0d){
					meanCounts.add(meanVal);
					sum+=meanVal;	//Used in finding the probabilities
				}
			}
			//Calculate refined permutation entropy
			double pe = 0;
			for (int p = 0; p<meanCounts.size(); ++p){
				meanCounts.set(p,meanCounts.get(p)/sum);	//calculate probabilities
				pe+=(-meanCounts.get(p)*Math.log(meanCounts.get(p)));	//Accumulate refined permutation entropy
			}
			return pe;
		}

		//helper class for sorting
		public class Sortable implements Comparable<Sortable>{
			public double a;	//Used to store the value
			public int b;	//Used store the original index of the value
			public Sortable(double a, int b){
				this.a = a;
				this.b = b;
			}
			@Override
			public int compareTo(Sortable b){
				return  this.a <= b.a ? -1 : 1;
			}
		}
	}
	
	private double mean(double[] a){
		double b = 0;
		for (int i = 0;i<a.length;++i){
			b+=a[i];
		}
		return b/((double) a.length);		
	}
	
	/**Helper functions for permutation*/
	/**
		Return all permutations of an array
		Implementation of Heap's algorithm https://en.wikipedia.org/wiki/Heap%27s_algorithm
		@params a the array to get permutations for. Note that the number of permutations is the factorial of the length of the array, i.e. n!
		@returns all permutations of array @a
	*/
	public static int[][] permutations(int[] a){
		int[][] b = new int[factorial(a.length)][a.length];
		//Heap's algorithm in non-recursive format
		int n = a.length;
		int[] c = new int[n];
		int currentPerm = 0;	//Keeps track of the current permutation number
		b = setRow(b,a,currentPerm);	//First permutation is the original
		int i = 0;
		while (i<(n)){
			if (c[i] < i){
				if (i % 2 == 0){
					a = swap(a,0,i);
				}else{
					a = swap(a,c[i],i);
				}
				++currentPerm;
				b = setRow(b,a,currentPerm);	//Set the permutation
				c[i]++;
				i = 0;
			}else{
				c[i] = 0;
				++i;
			}
		}
		return b;
	}

	/**
		Helper function to swap to items in an array
	*/
	public static int[] swap(int[] a, int b, int c){
		int tempB = a[b];
		a[b] = a[c];
		a[c] = tempB;
		return a;
	}
	
	/**
		Helper function to store the current permutation in the result array
	*/
	public static int[][] setRow(int[][] a,int[] b, int c){
		for (int i = 0;i<b.length;++i){
			a[c][i] = b[i];
		}
		return a;
	}
	
	/**
		Calculate the factorial (n! = 1*2* ... *n)
		@params a the value to calculate factorial for
		@returns factorial of @a
	*/
	public static int factorial(int a){
		int b;		
		b = 1;
		for (int i = 2; i<=a; ++i){
			b = b*i;
		}
		return b;
	}
	
	
}
