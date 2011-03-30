package fbi.genome.sequencing.rnaseq.reconstruction;

public class Kernel {

	public static final byte KERNEL_EPANECHNIKOV= 1;
	

	public static void main(String[] args) {
		getEpanechnikovKernel(10);
	}
	
	public static double[] getEpanechnikovKernel(int w) {
		double[] f= new double[w];
		// http://en.wikipedia.org/wiki/Kernel_(statistics)
		// 3/(2* windowsize) * (1- u^2)
		// u= -1+ (idx/ (w-1))* 2
		double csum= 0;
		for (int i = 0; i < f.length; ++i) {
			double u= ((2* i)/ (double) (w- 1d)) -1d;
			f[i]= (1- (u*u))* 3d/ (2d* w);
			csum+= f[i];
		}
		if (csum!= 1) {
			double csum2= 0;
			for (int i = 0; i < f.length; ++i) {
				f[i]/= csum;
				csum2+= f[i];
			}
			csum= csum2;
			//System.err.println(csum2);
		}
		// must be, sometimes its 0.9999999999...
		//assert(csum== 1);
		
		return f;
	}
	
	public static int smoothen(byte kernel, int w, int[] b) {
		
		double[] f= getKernel(kernel, w);
		if (f== null)
			return -1;
		
		double[] a= new double[b.length];
		for (int i = 0; i < a.length; i++) 
			a[i]= 0;
		
		for (int i = -w/ 2; i < a.length- (w/ 2); ++i) {
			int med= b[i+ (w/ 2)];
			for (int j = 0; j < w; ++j) {
				if (i+ j< 0|| i+ j>= a.length)
					continue;
				a[i+ j]+= f[j]* med;
			}
		}
		
		int sum= 0;
		for (int i = 0; i < a.length; i++) {
			b[i]= (int) Math.round(a[i]);
			sum+= b[i];
		}
		
		return sum;
	}
	
	public static double[] getKernel(byte kernel, int w) {
		if (kernel== KERNEL_EPANECHNIKOV)
			return getEpanechnikovKernel(w);
		System.err.println("Not implemented kernel "+ kernel);
		return null;
	}
	
}
