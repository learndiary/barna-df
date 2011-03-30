import fbi.genome.model.splicegraph.Graph;


public class Test {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		float f= 0.00001f;
//		String s= String.format("%1$f", f);
//
//		System.out.println(s+"*");
//		
//		System.err.println("3 mod 0"+(3%0));
		
		long[] a= new long[2];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < 64; j++) {
				a[i]|= (int) Math.pow(2, j);
			}
		}
		
		for(int idx= -1; (idx= Graph.getNextTxIdx(a, idx))!= -1; ){
			System.out.println(idx);
		}
		
		long x= 1;
		for (int i = 0; i < 64; i++, x*= 2) {
			System.err.println(i+"\t"+x);
		}
	}

	void kstest() {
		//StatisticalComparison.compare(null, null);
	}
}
