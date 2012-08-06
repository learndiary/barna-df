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

package barna.flux.capacitor.closure;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
class MultiFragExtension extends MultiFrag {
	
	protected String[] sequences= null;
	protected int matches= -1;
	protected byte[] matchTable= null;
	protected int applicables= -1;
	protected byte[] applicTable= null;

	public MultiFragExtension(MultiFrag base) {
		
		this.b= base.b;
		this.cs= base.cs;
		this.ext= base.ext;
		this.it= base.it;
		this.next= base.next;
		this.number= base.number;
		this.ow= base.ow;
		this.pred= base.pred;
		this.s= base.s;
		this.sel= base.sel;
		this.trans= base.trans;
		this.weight= base.weight;
	}

	public String[] getSequences() {
		return sequences;
	}
	
	public void setSequences(String seq1, String seq2) {
		
		this.sequences= new String[2];
		sequences[0]= seq1;
		sequences[1]= seq2;
		
		countMatches();
	}	
	
	public void setSequence(boolean first, String seq) {
		
		if (sequences== null)
			sequences= new String[2];
			
		if (first)
			sequences[0]= seq;
		else
			sequences[1]= seq;

		if (sequences[0]!= null&& sequences[1]!= null)
			countMatches();			
	}
	
	public int countMatches() {
		
		if (sequences== null|| sequences[0]== null|| sequences[1]== null
			|| sequences[0].length()!= sequences[1].length())
			return matches;
		
		matches= 0;
		matchTable= new byte[sequences[0].length()];
		for (int i= 0; i< sequences[0].length(); ++i) {
			if (Character.toUpperCase(sequences[0].charAt(i))==
				Character.toUpperCase(sequences[1].charAt(i))) {
				++matches;
				matchTable[i]= 1;
			} else
				matchTable[i]= 0;
		}
		
		return matches;
	}
	
	public byte[] getApplicTable() {
		return applicTable;
	}
	
	public void setApplicables(byte[] newApplicTable) {
		this.applicTable= newApplicTable;
		this.applicables= countApplicables(newApplicTable);
	}
	
	public int countApplicables(byte[] newApplicTable) {
		
		if (newApplicTable== null)
			return (-1);
		
		int counter= 0;
		for (int i= 0; i< newApplicTable.length; ++i)
			if (newApplicTable[i]== 1)
				++counter;
				
		return counter;
	}
	
	public String toString() {
		
		String result= super.toString()+ barna.commons.system.OSChecker.NEW_LINE;
		
		if (matchTable== null)
			return result;
		result+= "matches: "+ matches;
		result+= ",\tmatch ratio: "+ ((float) matches/ (float) ext);
		if (applicTable!= null) {
			result+= ";\t\tapplicables: "+ applicables;
			result+= ",\tapplicable ratio: "+ ((float) applicables/ (float) ext);
		}
		result+= barna.commons.system.OSChecker.NEW_LINE;
		
		result+= sequences[0]+ "   "+ sequences[0]+ barna.commons.system.OSChecker.NEW_LINE;
		for (int i= 0; i< matchTable.length; ++i)
			result+= (matchTable[i]== 0)? "X":"|";
		if (applicTable== null)
			return result;
		result+= "   ";
		for (int i= 0; i< applicTable.length; ++i)
			result+= (applicTable[i]== 0)? "*":((applicTable[i]== (-1))?"X":"|");
		result+= barna.commons.system.OSChecker.NEW_LINE+ sequences[1]+ "   "+ sequences[1]+ barna.commons.system.OSChecker.NEW_LINE;

		return result;
	}
}
