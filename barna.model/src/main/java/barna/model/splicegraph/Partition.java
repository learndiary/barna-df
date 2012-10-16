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

import java.util.HashMap;
import java.util.Iterator;


public class Partition {
	public HashMap<PartitionSet, PartitionSet> parents;
	public long[] transcripts;
	public Partition() {
		transcripts= null;
		parents= new HashMap<PartitionSet,PartitionSet>(4,1f);
	}
	public Partition(long[] tx) {
		this();
		this.transcripts= tx;
	}
	
	
	public void addParent(PartitionSet newParent) { 
		newParent.partitions.put(this,this);
		parents.put(newParent,newParent);
	}
	
	@Override
	public String toString() {
		StringBuilder sb= new StringBuilder("[");
		if (transcripts== null)
			sb.append("null]");
		else {
			for (int i = 0; i < transcripts.length; i++) 
				sb.append(Long.toString(transcripts[i])+",");
			sb.deleteCharAt(sb.length()- 1);
			sb.append("]");
		}
		return sb.toString();
	}
	
	@Override
	public int hashCode() {
		int sum= 0;
		for (int i = 0; transcripts!= null&& i < transcripts.length; i++) {
			sum+= transcripts[i];
		}
		return sum;
	}
	
	@Override
	public boolean equals(Object obj) {
		
		if (!(obj instanceof Partition))
			return false;
		
		Partition p= (Partition) obj;
		if (transcripts== null) {
			if (p.transcripts== null)
				return true;	// ??
			else 
				return false;
		}
		if (transcripts.length!= p.transcripts.length)
			return false;
		
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i]!= p.transcripts[i])
				return false;
		}
				
		return true;
	}
	
	public Partition clonePartitionWithoutTx() {
		
		Partition p= new Partition();
		p.parents= (HashMap<PartitionSet, PartitionSet>) parents.clone();
		Iterator<PartitionSet> iter= parents.keySet().iterator();
		while (iter.hasNext())
			iter.next().partitions.put(p, p);
		return p;
	}
	public HashMap<PartitionSet, PartitionSet> getParents() {
		return parents;
	}
	public void setParents(HashMap<PartitionSet, PartitionSet> parents) {
		this.parents = parents;
	}
	public long[] getTranscripts() {
		return transcripts;
	}
	public void setTranscripts(long[] transcripts) {
		this.transcripts = transcripts;
	}
}


