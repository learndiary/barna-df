/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.model.splicegraph;

import barna.model.Translation;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

public class PartitionCDS extends Partition {

	public boolean cdsValid53;
	public int currFrame5= Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI];
	public int frame3= Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI];
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
