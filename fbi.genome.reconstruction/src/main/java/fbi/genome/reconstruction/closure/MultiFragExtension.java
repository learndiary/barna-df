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

package fbi.genome.reconstruction.closure;

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
		
		String result= super.toString()+ "\n";
		
		if (matchTable== null)
			return result;
		result+= "matches: "+ matches;
		result+= ",\tmatch ratio: "+ ((float) matches/ (float) ext);
		if (applicTable!= null) {
			result+= ";\t\tapplicables: "+ applicables;
			result+= ",\tapplicable ratio: "+ ((float) applicables/ (float) ext);
		}
		result+= "\n";
		
		result+= sequences[0]+ "   "+ sequences[0]+ "\n";
		for (int i= 0; i< matchTable.length; ++i)
			result+= (matchTable[i]== 0)? "X":"|";
		if (applicTable== null)
			return result;
		result+= "   ";
		for (int i= 0; i< applicTable.length; ++i)
			result+= (applicTable[i]== 0)? "*":((applicTable[i]== (-1))?"X":"|");
		result+= "\n"+ sequences[1]+ "   "+ sequences[1]+ "\n";

		return result;
	}
}
