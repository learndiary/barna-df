package fbi.genome.io.bed;

import java.util.Iterator;

import fbi.genome.model.bed.BEDobject2;

public class BEDiteratorArray implements BufferedBEDiterator{
	
	BEDobject2[] a;
	int p, s;
	
	public BEDiteratorArray(BEDobject2[] a) {
		this.a= a;
		p= 0;
	}

	@Override
	public Iterator<BEDobject2> iterator() {
		return new BEDiteratorArray(a);
	}

	@Override
	public boolean hasNext() {		
		return (a!= null&& p< a.length);
	}

	
	@Override
	public BEDobject2 next() {
		
		if (a== null|| p>= a.length)
			return null;
		
		return a[p++];
	}
	
	/**
	 * @deprecated remove, not in Iterable/Iterator interface
	 * @return
	 */
	public int size() {
		if (a== null)
			return 0;
		return a.length;
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
	}
	
	@Override
	public void mark() {
		s= p;
	}
	
	@Override
	public boolean reset() {
		
		if(s< 0|| s>= a.length)
			return false;
		p= s;
		return true;
	}
}
