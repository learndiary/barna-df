package fbi.commons;

public class MyArrays {
	int[] ia; long[] la; double[] da; float[] fa; byte[] ba;
	int len;
	public MyArrays(Object array) {
		if (array instanceof int[]) {
			ia= (int[]) array;
			len= ia.length;
		} else if (array instanceof long[]) {
			la= (long[]) array;
			len= la.length;
		} else if (array instanceof double[]) {
			da= (double[]) array;
			len= da.length;
		} else if (array instanceof float[]) {
			fa= (float[]) array;
			len= fa.length;
		} else if (array instanceof byte[]) {
			ba= (byte[]) array;
			len= ba.length;
		}
	}
	public int length() {
		return len;
	}
	public Number elementAt(int pos) {
		if (ia!= null)
			return ia[pos];
		if (la!= null)
			return la[pos];
		if (da!= null)
			return da[pos];
		if (ba!= null)
			return ba[pos];
		//if (fa!= null)
		return fa[pos];		
	}
}
