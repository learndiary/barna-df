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

package barna.model.bed;

import barna.commons.ByteArrayCharSequence;
import barna.model.Graph;
import barna.model.Mapping;

import java.util.Comparator;

public class BEDobject2 extends ByteArrayCharSequence implements Mapping{

	public static final byte BYTE_PLUS= 43, BYTE_COMMA= 44, BYTE_MINUS= 45, BYTE_DOT= 46, BYTE_ZERO= 48;
	
	// chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRGB, bCount, bSizes, bStarts
	public static final byte[] DEFAULT_VALUE= new byte[] {BYTE_DOT, BYTE_ZERO, BYTE_ZERO,
		BYTE_DOT, BYTE_ZERO, BYTE_DOT, 		// name, score, strand
		BYTE_ZERO, BYTE_ZERO, BYTE_ZERO,	// thick start/end, itemRGB 
		BYTE_ZERO, BYTE_ZERO, BYTE_ZERO};	// bCount, sizes, starts
	
	// 0-based
	public static final int FN_NAME= 3, FN_BLOCK_SIZES= 10, FN_BLOCK_STARTS= 11;	
	public static final int DEFAULT_CAPACITY= 150;

	public static final BedIDComparator DEFAULT_ID_COMPARATOR = new BedIDComparator();
	
	int bedStart= -1, bedEnd= -1, score= -1;
	byte strand= Byte.MIN_VALUE, blockCount= -1;

    int chrP2= -1,
		nameP1= -1, nameP2= -1; 
	int blockSizeP1= -1, blockSizeP2= -1,
		lastBsize= -1, lastBstart= -1;	// can exceed byte
	BEDobject2 next= null;
	
	public static class BedIDComparator implements Comparator<BEDobject2> {
		public int compare(BEDobject2 o1, BEDobject2 o2) {			
			return o1.getName().compareTo(o2.getName());
		}
	}
	
	
	public BEDobject2() {
		super(DEFAULT_CAPACITY);
	}
	public BEDobject2(int capacity) {
		super(capacity);
	}
	
	public BEDobject2(ByteArrayCharSequence cs) {
		super(cs);
	}
	
	public void clear() {
		resetFind();		// p1= p2= cnt= 
		start= end= 0;		// ByteArrayCharSequence		
		bedStart= bedEnd= score= -1;	// BED object
		strand= Byte.MIN_VALUE;
		chrP2= nameP1= nameP2= (byte) (blockSizeP1= blockSizeP2= -1);
		resetBlocks();		
	}
	
	public void resetBlocks() {
		lastBsize= lastBstart= -1;
	}
	
	public void init() {
		if (end!= 0) {
			// 3 mandatory fields
			resetFind();
			find(0);
			if (cnt== 0)
				chrP2= p2;
			else 
				return;
			bedStart= getTokenInt(1);
			if (bedStart< 0) {
				bedStart= -1;
				return;
			}
			bedEnd= getTokenInt(2);
			if (bedEnd< 0) {
				bedEnd= -1;
				return;
			}
			
			// 9 opt fields
			find(3);	// name
			if (cnt== 3) {
				nameP1= p1;
				nameP2= p2;
			} else
				return;
			
			score= getScore();
			if (score< 0) {
				score= -1;
			}
			
			// field 5
			if (p2+ 1>= end)
				return;
			if (chars[p2+ 1]== BYTE_PLUS)
				strand= 1;
			else if (chars[p2+ 1]== BYTE_MINUS)
				strand= -1;
			else 
				strand= 0;
			
			// 6: thickStart, 7: thickEnd, 8: itemRGB
			if (p2+ 2>= end)
				return;
			blockCount= 0;
			try {
				blockCount= (byte) getTokenInt(9);
			} catch (IllegalArgumentException e) {
				blockCount= 0;
			}
			if (blockCount< 0) {
				blockCount= 0;	// no blocks
				return;
			}
			
			find(10);
			if (cnt== 10) {
				blockSizeP1= p1;
				if (p2< end)
					blockSizeP2= p2;
			}
		}
	}
	
	public boolean isInited() {
		return (bedStart>= 0&& bedEnd>= 0&& chrP2>= 0&& nameP1>= 0&& nameP2>= 0);
	}
	
	public void setChromosome(CharSequence chr) {
		if (end== start) {
			assert(end== start);
			append(chr);
			ensureLength(start, 1);
			cnt= 0;
			p1= start;
			p2= end;
			chrP2= end;
			chars[end++]= BYTE_TAB;	// mandatory, add fs after
		} else {
			p1= start;
			p2= chrP2;
			cnt= 0;
			replaceCurrField(chr);
		}
	}
	
	/**
	 * Sets the &quot;start&quot;field of the bed object.<br>
	 * @param x new start position
	 */
	public void setStart(int x) {
		bedStart= x;
		if (countTokens()< 2) {
			p1= end+ 1;
			append(x);
			ensureLength(end, 1);
			chars[end++]= BYTE_TAB;	// mandatory, add fs after
			p2= end;
            cnt= 1;
		} else {
			p1= start;
			p2= chrP2;
			cnt= 0;
			replace(1, start);
		}
	}
	
	/**
	 * Sets the &quot;end&quot;field of the bed object.<br>
	 * <b>Note</b>: counts the number of fields that are initialized 
	 * so far, consider calling append when successively building 
	 * a BED string.
	 * @param x new end position
	 */
	public void setEnd(int x) {
		bedEnd= x;
		if (countTokens()< 3) {
			assert(cnt== 1&& p2== end);
			p1= end+ 1;
			append(x);
			ensureLength(end, 1);
			p2= end;
			cnt= 2;
			// last mandatory, no fs after
		} else {
			p1= start;
			p2= chrP2;
			cnt= 0;
			replace(2, x);
		}
	}
	
	// optional fields
	public void setName(CharSequence name) {
		if (cnt== 2&& p2== end) {
			assert(nameP1< 0|| nameP2< 0);
			ensureLength(end, 1);
			chars[end++]= BYTE_TAB;	// optional, fs before
			nameP1= (p1= end);
			append(name);
			nameP2= (p2= end);
			cnt= 3;
		} else {
			if (!isInited())
				init();
			p1= nameP1; p2= nameP2;
			cnt= 3;
			replaceCurrField(name);
			nameP1= p1;
			nameP2= p2;
		}
	}
	
	public void setScore(int score) {
		this.score= score;
		if (cnt== 3&& p2== end) {
			ensureLength(end, 1);
			chars[end++]= BYTE_TAB;
			p1= end;
			append(score);
			p2= end;
			cnt=4;
		} else {
			if (nameP1>= 0&& nameP2>= 0) {
				p1= nameP1;
				p2= nameP2;
				cnt= 3;
			}
			find(4);
			if (cnt== 4)
				replaceCurrField(score);
		}
	}
	
	public void setStrand(byte strand) {
		this.strand= strand;
		byte strandByte= strand== 1? BYTE_PLUS: (strand== -1? BYTE_MINUS: BYTE_DOT);
		if (cnt== 4&& p2== end) {
			ensureLength(end, 2);
			chars[end++]= BYTE_TAB;
			p1= end;
			chars[end++]= strandByte;
			p2= end;
			cnt= 5;
		} else {
			if (nameP1>= 0&& nameP2>= 0) {
				p1= nameP1;
				p2= nameP2;
				cnt= 3;
			}
			find(5);
			if (cnt== 5) {
				replaceCurrField(strandByte);
			}
		}
	}
	
	public void readSequence(ByteArrayCharSequence cs) {
		int bedStart= getStart(); 
		try {
			int f= getStrand()>= 0? 1: -1;
			CharSequence chrom= getChr();
			int headerLen= cs.end;
			int bcount= getBlockCount();

            // contingent read
			if (bcount<= 1) {
				int first= f* (bedStart+ 1), 
					last= f* bedEnd,
					len= Math.abs(last- first)+ 1; 
				cs.ensureLength(cs.end, len);
				Graph.readSequence(chrom, getStrand()>= 0, 
						first, last, 
						cs, cs.end, cs.end+ len);
				cs.end+= len;

            // split read
			} else for (int i = 0; i < bcount; i++) { 
				int nextSt= getNextBlockStart(),
					nextSi= getNextBlockSize();
				assert(nextSt>= 0|| nextSi>= 0);
				int first= f* (bedStart+ nextSt+ 1), 
					last= f* (bedStart+ nextSt+ nextSi),
					len= Math.abs(last- first)+ 1; 
				if (getStrand()>= 0) {
					cs.ensureLength(cs.end, len);
					Graph.readSequence(chrom, getStrand()>= 0, 
							first,	// f* (start+ getBlockStart(i)+ 1) 
							last,	// f* (start+ getBlockStart(i)+ getBlockSize(i)) 
							cs, cs.end, cs.end+ len);
					cs.end+= len;
				} else {
					cs.ensureLength(cs.end, len);
					cs.end+= len;
					System.arraycopy(cs.chars, headerLen, cs.chars, headerLen+ len, cs.end- headerLen- len);
					Graph.readSequence(chrom, getStrand()>= 0, 
							first, last, 
							cs, headerLen, headerLen+ len);
				}
			}
		} catch (Exception e) {
			throw new RuntimeException("Problems reading BED object "+ toString(), e);
		}
	}

	
	public void setBlockCount(int count) {
		this.blockCount= (byte) count;
		if (cnt< 3&& nameP1>= 0&& nameP2>= 0) {
			p1= nameP1; p2= nameP2; cnt= 3;
		}
		if (cnt!= 9)
			find(9);
		if (cnt== 9)
			replaceCurrField(count);
		else {
			assert(p1== end);
			ensureLength(end, 1);
			chars[end++]= BYTE_TAB;
			p1= end;
			append(count);
			p2= end;
			cnt= 9;
		}
	}

	/**
	 * puts end right before the field to come
	 * @param fn
	 */
	void extend2field(int fn) {
		ensureLength(end, (fn- cnt)); 	// -1 +1
		for (; cnt < fn; ++cnt) {
			//p1= end;
			//a[end++]= DEFAULT_VALUE[cnt];
			chars[end++]= BYTE_TAB;
			p1= p2= end;	// empty fields
		}
	}
	
	public void setNextBlockSize(int x) {
		if (blockSizeP1< 0|| blockSizeP2< 0) 
			forceBlocksField();
		p1= blockSizeP1; 
		p2= blockSizeP2; 
		cnt= FN_BLOCK_SIZES;
		if (p2> p1) 
			appendCurrField(BYTE_COMMA);
		appendCurrField(x);
		blockSizeP2= p2;
	}

	private void forceBlocksField() {
		find(FN_BLOCK_STARTS);
		if (cnt< FN_BLOCK_STARTS)
			extend2field(FN_BLOCK_STARTS);
		find(FN_BLOCK_SIZES);
		assert(cnt== FN_BLOCK_SIZES);
		blockSizeP1= p1;
		blockSizeP2= p2;
	}
	
	public void setNextBlockStart(int x) {
		if (blockSizeP1< 0|| blockSizeP2< 0) 
			forceBlocksField();
		p1= blockSizeP2+ 1; 
		p2= end; 
		cnt= FN_BLOCK_STARTS;
		
		if (p2> p1)
			appendCurrField(BYTE_COMMA);
		appendCurrField(x);
		
	}
	
	public int getNextBlockSize() {
		if (blockSizeP1< 0|| blockSizeP2< 0) {
			find(FN_BLOCK_SIZES);
			if (cnt!= FN_BLOCK_SIZES)
				return -1;
		} 
		p1= blockSizeP1; 
		p2= blockSizeP2; 
		cnt= FN_BLOCK_SIZES;
		
		if (lastBsize== p2)
			lastBsize= -1;
		int x= lastBsize= Math.max(p1, lastBsize+ 1);
		while (++lastBsize < p2) {
			if (chars[lastBsize]== BYTE_COMMA)
				break;
		}
		return parseInt(x, lastBsize);
	}

    @Override
    public CharSequence getSequence() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public CharSequence getCigar() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double getCount(boolean weighted) {
        return -1;
    }

    @Override
    public byte getReadStrand(String readStrand) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getNextBlockStart() {
		if (blockSizeP1< 0|| blockSizeP2< 0) {
			find(FN_BLOCK_SIZES);
			if (cnt!= FN_BLOCK_SIZES)
				return -1;
		} 
		p1= blockSizeP2+ 1; 
		p2= end; 
		cnt= FN_BLOCK_STARTS;
		if (lastBstart== p2) 
			lastBstart= -1;
		int x= lastBstart= Math.max(p1, lastBstart+ 1);
		while (++lastBstart < p2) {
			if (chars[lastBstart]== BYTE_COMMA)
				break;
		}
		return parseInt(x, lastBstart);
	}

	public byte getStrand() {
		if (!isInited())
			init();
		return strand;
	}

    @Override
    public byte getMateFlag() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public ByteArrayCharSequence getName() {
		if (!isInited())
			init();
		if (nameP1< 0|| nameP2< 0)
			return null;
		return subSequence(nameP1, nameP2);
	}
	
	public int getStart() {
		if (!isInited())
			init();
		return bedStart;
	}
	
	public int getLength() {
		if (!isInited())
			init();
		if (blockSizeP1> 0) {
			int sum= 0;
			for (int i = 0; i < getBlockCount(); i++) 
				sum+= getNextBlockSize();	
			resetBlocks();
			return sum;
		} else 
			return getEnd()- getStart();
		
	}
	
	public int getEnd() {
		if (!isInited())
			init();
		return bedEnd;
	}

	public ByteArrayCharSequence getChr() {
		if (!isInited())
			init();
		if (chrP2< start)
			return null;
		return subSequence(start, chrP2);
	}
	
	public int getBlockCount() {
		if (!isInited())
			init();
		return blockCount;
	}

	@Override
	protected boolean find(int fieldNr) {
        boolean b = super.find(fieldNr);
        if (fieldNr== FN_NAME) {
			nameP1= p1; nameP2= p2;
		}
        return b;
	}
	
	public int getNameP2() {
		if (nameP2< 0) {
			if (cnt!= FN_NAME)
				find(FN_NAME);
		}
		return nameP2;
	}
	
	public int getNameP1() {
		if (nameP1< 0) {
			if (cnt!= FN_NAME)
				find(FN_NAME);
		}
		return nameP1;
	}
	
	@Override
	public String toString() {
		return super.toString();
	}
	public int getScore() {
		
		int x= -1;
		
		if (find(4)) {
			if (p1+ 1== p2&& chars[p1]== 46)
				return 0; // value in case of a '.'
		} else 
			return 0;	// value in case of no score field
		
		try {
			x= getTokenInt(4);
		} catch (NumberFormatException e) {
			float f= getTokenFloat(4);
			x= (int) f; 
		}
		
		return x; 
	}
	public BEDobject2 getNext() {
		return next;
	}
	public void setNext(BEDobject2 next) {
		this.next = next;
	}


    @Override
    public CharSequence getName(Boolean appendMateNumber) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public CharSequence getChromosome() {
        return getChr();
    }
}
