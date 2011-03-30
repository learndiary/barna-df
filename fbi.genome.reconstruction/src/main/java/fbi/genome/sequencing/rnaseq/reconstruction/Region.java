package fbi.genome.sequencing.rnaseq.reconstruction;

import java.util.Comparator;
import java.util.Vector;

public class Region {
	
	public static class PositionComparator implements Comparator<Region> {
		//@Override
		public int compare(Region o1, Region o2) {
			if (o1.exStart< o2.exStart)
				return -1;
			if (o1.exStart> o2.exStart)
				return 1;
			if (o1.exEnd< o2.exEnd)
				return -1;
			if (o1.exEnd> o2.exEnd)
				return 1;
			return 0;
		}
	}
	
	public static PositionComparator defaultPositionComparator= new PositionComparator(); 
	
	public int exStart, exEnd;
	Vector<Variation> inVarVec= new Vector<Variation>(), outVarVec= new Vector<Variation>();
}
