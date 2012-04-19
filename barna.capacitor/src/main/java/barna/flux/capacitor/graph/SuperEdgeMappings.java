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

package barna.flux.capacitor.graph;

import barna.model.Transcript;
import barna.model.constants.Constants;
import barna.model.splicegraph.AbstractEdge;
import barna.model.splicegraph.SuperEdge;

public class SuperEdgeMappings extends SuperEdge implements MappingsInterface {

	Mappings mappings= null;
	
	public SuperEdgeMappings(AbstractEdge[]  edges, long[] support, boolean pairedEnd) {
		super(edges, support, pairedEnd);
		mappings= new Mappings();
	}
	
	
	public float[] getCoverage(boolean sense, int minMapLen, int maxMapLen) {
//		if (isPend()) {
//			if (sense)
//				return edges[0].getCoverage(sense, minMapLen, maxMapLen);
//			else
//				return edges[edges.length- 1].getCoverage(sense, minMapLen, maxMapLen);
//		} else
		//return super.getCoverage(sense, minMapLen, maxMapLen);
		return null;
	}

	public float getCoverage(int p, boolean sense, int minMapLen, int maxMapLen) {
		
		assert(isExonic());
		
		int start= getGpos(sense, true, minMapLen, maxMapLen);
		int end= getGpos(sense, false, minMapLen, maxMapLen);
		if(p< start|| p> end)
			return -1;
		
		float[] a= getCoverage(sense, minMapLen, maxMapLen);
		if (a== null)
			return -1;
		
		int sStart= getGpos(sense, true, minMapLen, maxMapLen), sEnd= getGpos(sense, false, minMapLen, maxMapLen);
		if (p< sStart|| p> sEnd)
			return -1;
		if ((sense&& mappings.coverage== null)|| ((!sense)&& mappings.coverageRev== null))
			return 0;
		return sense? mappings.coverage[p- sStart]:mappings.coverageRev[p- sStart]; 
		
	}

	@Override
	public int getMapLength(Transcript tx, int maxMapLength) {
		if (isPend()) {
			System.err.println("what about pends");
			return -1;
			
		} else {
			int start= tx.getExonicPosition(getEdges()[0].getTail().getSite().getPos()),
				first= tx.getExonicPosition(getEdges()[0].getHead().getSite().getPos()),
				last= tx.getExonicPosition(getEdges()[getEdges().length- 1].getTail().getSite().getPos()),
				end= tx.getExonicPosition(getEdges()[getEdges().length- 1].getHead().getSite().getPos());
			return Math.min(Math.min(last- start, end- first)+ 1, maxMapLength);
		}
	}

	public double getNtCoverage(Transcript tx, byte dir, int mapLenMax) {
		
		if (isPend()) {
			System.err.println("check pend");
			double cov= mappings.getReadNr();
			if (dir== Constants.DIR_BOTH) {
				cov*= 2;
				cov/= getEffLength(tx, Constants.DIR_FORWARD, mapLenMax)+ 
						getEffLength(tx, Constants.DIR_BACKWARD, mapLenMax); 
			} else if (dir== Constants.DIR_FORWARD)
				cov/= getEffLength(tx, Constants.DIR_FORWARD, mapLenMax);
			else if (dir== Constants.DIR_BACKWARD)
				cov/= getEffLength(tx, Constants.DIR_FORWARD, mapLenMax);
			else 
				assert(false);
			return cov;
		} else
			//return super.getNtCoverage(tx, dir, mapLenMax);
			// TODO
			return -1d;
		
	}


	@Override
	public Mappings getMappings() {
		return mappings;
	}

	
}
