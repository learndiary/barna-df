package barna.astalavista.statistics;

import java.util.Arrays;

public class PreCluster {

	static boolean normal= true;	// only for debug purposes
	
	int[] a;	// array with indices of rowNr that are in cluster
	
	public PreCluster(int[] a) {
		this.a= a;
	}
	
	@Override
	public int hashCode() {
//		int hash= a.hashCode(); 	// not toString().hashCode(), suffix "..@memAddr"
		int hash= Arrays.hashCode(a);
		return hash;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof PreCluster))
			return false;
		
		int[] b= ((PreCluster) obj).a;
//		if (a.length!= b.length)
//			return false;
//		for (int i = 0; i < b.length; i++) {
//			if (a[i]!= b[i])
//				return false;
//		}
//		
//		return true;
		
		if (normal)
			return Arrays.equals(a, b);
		else 
			return false;
	}
}
