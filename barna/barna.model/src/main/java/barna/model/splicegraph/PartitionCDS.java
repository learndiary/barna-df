/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.model.splicegraph;

import barna.model.Translation;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

public class PartitionCDS extends Partition {

	boolean cdsValid53;
	int currFrame5= Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI];
	int frame3= Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI];
	HashMap<Integer, Vector<PartitionCDS>> mapConfirmedCDSflanks;
	
	public PartitionCDS() {
		super();
		currFrame5= frame3= Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI];
		cdsValid53= false;
	}
	public PartitionCDS(int frame3) {
		this();
		this.frame3= frame3;
	}
	public PartitionCDS(long[] tx, int frame3) {
		this();
		this.transcripts= tx;
		this.frame3= frame3;
	}

	@Override
	public boolean equals(Object obj) {
		
		if (!(obj instanceof PartitionCDS))
			return false;
		
		if (!super.equals(obj))
			return false;
		
		PartitionCDS p= (PartitionCDS) obj;
		if (p.currFrame5!= currFrame5
				|| p.frame3!= frame3
				|| p.cdsValid53!= p.cdsValid53)
			return false;
		return true;
	}
	
	@Override
	public Partition clonePartitionWithoutTx() {
		PartitionCDS p= new PartitionCDS();
		p.parents= (HashMap<PartitionSet, PartitionSet>) parents.clone();
		Iterator<PartitionSet> iter= parents.keySet().iterator();
		while (iter.hasNext())
			iter.next().partitions.put(p, p);
		p.frame3= frame3;
		p.currFrame5= currFrame5;
		p.cdsValid53= cdsValid53;
		return p;
	}
	public int getCurrFrame5() {
		return currFrame5;
	}
	public void setCurrFrame5(int currFrame5) {
		this.currFrame5 = currFrame5;
		int combi= Translation.getCombinedFrame(getCurrFrame5(), getFrame3());
		if (combi== -65528)
			System.currentTimeMillis();
	}
	public int getFrame3() {
		return frame3;
	}
	public void setFrame3(int frame3) {
		this.frame3 = frame3;
		int combi= Translation.getCombinedFrame(getCurrFrame5(), getFrame3());
		if (combi== -65528)
			System.currentTimeMillis();
	}
	
	public void addConfirmedCDSflanks(PartitionCDS anotherPartition, int combined) {
		if (mapConfirmedCDSflanks== null) {
			mapConfirmedCDSflanks= new HashMap<Integer, Vector<PartitionCDS>>(2);
		}
		Vector<PartitionCDS> v= mapConfirmedCDSflanks.get(combined);
		if (v== null) {
			v= new Vector<PartitionCDS>(1,1);
			mapConfirmedCDSflanks.put(combined, v);
		}
		v.add(anotherPartition);
	}
	
	public void removeConfirmedCDSflanks() {
		if (mapConfirmedCDSflanks== null|| mapConfirmedCDSflanks.size()== 0)
			return;
		mapConfirmedCDSflanks.clear();		
		assert(mapConfirmedCDSflanks.size()== 0);
	}
	public HashMap<Integer, Vector<PartitionCDS>> getMapConfirmedCDSflanks() {
		return mapConfirmedCDSflanks;
	}
	
	
}
