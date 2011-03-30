package fbi.genome.sequencing.rnaseq.reconstruction;

import java.util.Comparator;

public class Tuple {

	public static class TupleByXComparator implements Comparator<Tuple> {
		public int compare(Tuple o1, Tuple o2) {
			return (o1.x- o2.x);
		}
	}
	
	public TupleByXComparator defaultTupleByXComparator= new TupleByXComparator();
	
	public int x, y;
	public Tuple(int x, int y) {
		this.x= x;
		this.y= y;
	}
}
