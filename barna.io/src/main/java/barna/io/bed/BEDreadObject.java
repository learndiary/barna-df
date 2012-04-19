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

package barna.io.bed;

import barna.model.bed.BEDobject;

public class BEDreadObject extends BEDobject {
	public static char SEP_PE= '_', SEP= ':';
	byte pEnd= 0, mm= -1;
	int tot= -1;
	
	public BEDreadObject(String chromName, byte newStrand, int newStart, int newEnd) {
		super(chromName, newStrand, newStart, newEnd);
	}
	
	@Override
	public void setName(String name) {
		
		int lp= name.length(), p= lp- 1;
		while (p>= 0&& mm< 0) {
			while(p>= 0&& name.charAt(p)!=SEP&& name.charAt(p)!=SEP_PE)
				--p;
			if (p< 0)
				break;
			if (name.charAt(p)== SEP_PE)
				pEnd= Byte.parseByte(name.substring(p+1,lp));
			else if (tot< 0)
				tot= Integer.parseInt(name.substring(p+1,lp));
			else 
				mm= Byte.parseByte(name.substring(p+1,lp));
			lp= p;
		}		
		
		super.setName(name.substring(0,p));
	}

	public byte getPEnd() {
		return pEnd;
	}

	public byte getMm() {
		return mm;
	}

	public int getTot() {
		return tot;
	}
}
