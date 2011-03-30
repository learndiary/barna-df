package fbi.genome.sequencing.rnaseq.reconstruction;

public class BinVector {

	public int[] ax, ay;
	public int size, mass;
	
	public BinVector() {
		ax= new int[100];
		ay= new int[100];
		size= 0;
		mass= 0;
	}
	
	public void incrTuple(int len) {
		
		++mass;
		int p= 0;
		if (size> 0)
			p= binarySearch(ax, 0, size, len);
		else {
			ax[0]= len;
			ay[0]= 1;
			++size;
			return;
		}
		if (p>= 0) 
			++ay[p];
		else {
			p= -p- 1;	
			if (size+ 1>= ax.length)
				extend();
			System.arraycopy(ax, p, ax, p+1, size- p+ 1);
			System.arraycopy(ay, p, ay, p+1, size- p+ 1);
			++size;
			ax[p]= len;
			ay[p]= 1;
		}
		
	}
	
	public int getQuartile(float x) {
		float nr= x* mass, sum= 0;
		int i;
		for (i = 0; sum< nr&& i < size; sum+= ay[i++]); 
		return ax[i];
	}
	
	private void extend() {
		int[] newx= new int[ax.length* 2];
		System.arraycopy(ax, 0, newx, 0, size);
		int[] newy= new int[ay.length* 2];
		System.arraycopy(ay, 0, newy, 0, size);
		ax= newx;
		ay= newy;
		System.gc();
	}
	
	public static int binarySearch(int[] a, int fromIndex, int toIndex,
			   long key) {
		rangeCheck(a.length, fromIndex, toIndex);
		return binarySearch0(a, fromIndex, toIndex, key);
	}
	private static void rangeCheck(int arrayLen, int fromIndex, int toIndex) {
        if (fromIndex > toIndex)
            throw new IllegalArgumentException("fromIndex(" + fromIndex +
                       ") > toIndex(" + toIndex+")");
        if (fromIndex < 0)
            throw new ArrayIndexOutOfBoundsException(fromIndex);
        if (toIndex > arrayLen)
            throw new ArrayIndexOutOfBoundsException(toIndex);
    }
	
	// Like public version, but without range checks.
	private static int binarySearch0(int[] a, int fromIndex, int toIndex,
				     long key) {
		int low = fromIndex;
		int high = toIndex - 1;
		
		while (low <= high) {
			 int mid = (low + high) >>> 1;
			 long midVal = a[mid];
			
			 if (midVal < key)
				low = mid + 1;
			 else if (midVal > key)
				high = mid - 1;
			 else
				return mid; // key found
		}
		return -(low + 1);  // key not found.
	}

	
	public String toString() {
		
		StringBuilder sb= new StringBuilder();
		for (int i = 0; i < size; i++) {
			sb.append(ax[i]);
			sb.append("\t");
			sb.append(ay[i]);
			sb.append("\n");
		}
		
		return sb.toString();
	}
}
