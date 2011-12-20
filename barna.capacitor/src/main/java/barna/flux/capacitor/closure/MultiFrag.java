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

package barna.flux.capacitor.closure;

/**
 * fragments outside function `pairalign'.
 * 
 * @author micha
 */
public class MultiFrag {

		// b[0], b[1]:  begin of the diagonal
	protected int[] b= new int[2];
		// s[0], s[1]:  sequences, to which diagonal belongs
	protected int[] s= new int[2];
		// ext:         length of the diagonal
	protected int ext;
		// it:          iteration step 
	protected int it;
		// weight:      individual weight of the diagonal
	protected float weight;
		// ow:          overlap weight of the diagonal
	protected float ow; 
		// sel:         1, if accepted in filter proces, 0 else
	protected short sel;
		// trans:       translation
	short trans;
		// cs:          crick strand 
	short cs; 
		// *next:       next diagonal 
	protected MultiFrag next;
		// predeciding diagonal 
	protected MultiFrag pred;


		// ADDITIONAL ATTRIBUTES
	protected int number= 0;

	public boolean isAccepted() {
		
		if (sel== 1)
			return true;
			
		return false;
	}
	
	public boolean isConsistent() {
		
		return isAccepted();
	}
	
	public void setAccepted(boolean newAcception) {
		
		if (newAcception)
			sel= 1;
		else
			sel= 0;
	}
	
	public void setConsistent(boolean newConsistent) {
		
		setAccepted(newConsistent);
	}
		/**
		 * Returns the ow.
		 * @return float
		 */
		public float getOverlapWeight() {
			return ow;
		}

		/**
		 * Returns the weight.
		 * @return float
		 */
		public float getWeight() {
			return weight;
		}

		/**
		 * Sets the ow.
		 * @param ow The ow to set
		 */
		public void setOverlapWeight(float ow) {
			this.ow= ow;
		}

		/**
		 * Sets the weight.
		 * @param weight The weight to set
		 */
		public void setWeight(float weight) {
			this.weight= weight;
		}
		
		public String toString() {
			
			String result= getNumber()+ ") seq: ";
			if (getSequenceNo(true)< 10)
				result+= " "+ getSequenceNo(true);
			else
				result+= getSequenceNo(true);
			result+= " ";
			if (getSequenceNo(false)< 10)
				result+= " "+ getSequenceNo(false);
			else
				result+= getSequenceNo(false);
			result+= "\t";

			result+= "beg: ";
			if (getSequenceStart(true)< 10)
				result+= "  "+ getSequenceStart(true);
			else if (getSequenceStart(true)< 100)
				result+= " "+ getSequenceStart(true);
			else 
				result+= getSequenceStart(true);
			result+= " ";
			if (getSequenceStart(false)< 10)
				result+= "  "+ getSequenceStart(false);
			else if (getSequenceStart(false)< 100)
				result+= " "+ getSequenceStart(false);
			else 
				result+= getSequenceStart(false);
			result+= "\t";
			
			result+= "len: ";
			if (getLength()< 10)
				result+= "  "+ getLength();
			else if (getLength()< 100)
				result+= " "+ getLength();
			else 
				result+= getLength();
			result+= "\t";

			result+= "wgt: "+ getWeight()+ "\t";
			result+= "olw: "+ getOverlapWeight()+ "\t";
			result+= "it: "+ getIteration()+ "\t";
			
			if (isConsistent())
				result+= "cons  \t";
			else
				result+= "incons\t";
				
				// ??! check with Burkhard
			if (trans== 0)
				result+= "N-frg\t";
			else
				result+= "?-frg\t";
				
			return result;
		}
		
		public void setSequenceNos(int seqNo1, int seqNo2) {
			
			s[0]= seqNo1;
			s[1]= seqNo2;
		}
		
		public void setSequenceNo(boolean firstSeq, int seqNo) {
			
			if (firstSeq)
				s[0]= seqNo;
			else 
				s[1]= seqNo;
		}
		
		public void setSequenceNos(int[] newSeqNos) {
			
				// error, don't set
			if (newSeqNos.length!= 2)
				return;
				
			s= newSeqNos;
		}
		
		public int[] getSequenceNos() {
			
			return s;
		}
		
		public int getSequenceNo(boolean firstSeq) {
			
			if (firstSeq)
				return s[0];
			return s[1];
		}
		
		public void setSequenceStarts(int start1, int start2) {
		
				// create new, if necessary
			if (b== null)
				b= new int[2];
				
			b[0]= start1;
			b[1]= start2;
		}
		
		public void setSequenceStart(boolean newSeq, int newStart) {
			
				// create new, if necessary
			if (b== null)
				b= new int[2];
				
			if (newSeq)
				b[0]= newStart;
			else
				b[1]= newStart;
		}
		
		public int[] getSequenceStarts() {
			
			return b;
		}
		
		public int getSequenceStart(boolean firstSeq) {
			
			if (firstSeq)
				return b[0];
			return b[1];
		}		

		/**
		 * Returns the ext.
		 * @return int
		 */
		public int getLength() {
			return ext;
		}

		/**
		 * Sets the ext.
		 * @param ext The ext to set
		 */
		public void setLength(int ext) {
			this.ext= ext;
		}

		/**
		 * Returns the it.
		 * @return int
		 */
		public int getIteration() {
			return it;
		}

		/**
		 * Sets the it.
		 * @param it The it to set
		 */
		public void setIteration(int it) {
			this.it= it;
		}

		/**
		 * Returns the trans.
		 * @return short
		 */
		public short getTranslation() {
			return trans;
		}

		/**
		 * Sets the trans.
		 * @param trans The trans to set
		 */
		public void setTranslation(short trans) {
			this.trans = trans;
		}

		public void setTranslation(int trans) {
			
			setTranslation((short) trans);
		}

		/**
		 * Returns the number.
		 * @return int
		 */
		public int getNumber() {
			return number;
		}

		/**
		 * Sets the number.
		 * @param number The number to set
		 */
		public void setNumber(int number) {
			this.number = number;
		}

}
