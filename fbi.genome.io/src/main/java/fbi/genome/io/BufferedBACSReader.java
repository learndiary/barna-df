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

package fbi.genome.io;

import fbi.commons.ByteArrayCharSequence;
import fbi.genome.model.constants.Constants;

import java.io.IOException;
import java.io.InputStream;

public class BufferedBACSReader {
	
	private static final int DEFAULT_BUFFER_SIZE = 1024;
	public static void main(String[] args) {
		
	}
	
	public static byte BYTE_NL= '\n', BYTE_CR= '\r';
	int minVol= 100;
	InputStream in;
	byte[] buf;
	int pos= 0;
	long maxBytes= 0l, totBytes= 0l;
	
	public BufferedBACSReader(InputStream in) {
		this(in, DEFAULT_BUFFER_SIZE);
	}
	
	public BufferedBACSReader(long max, InputStream in) {
		this(in);
		this.maxBytes= max;
		totBytes= 0;
	}
	
	public BufferedBACSReader(InputStream in, int capacity) {
		this.in= in;
		buf= new byte[capacity];
		pos= 0;
	}
	
	public int readLine(ByteArrayCharSequence cs) {
		byte n= BYTE_NL, r= BYTE_CR;
		byte[] b= buf, a= cs.chars; // , bb= null;
		synchronized (b) {
//			while (pos[idx]== 0) 
//				fill(idx);
			if (pos< a.length) {
				fill();
			}
			int p= pos;
			if (p== 0) {
				return - 1;
			}
			int i = 0;
			for (; i < p&& b[i]!= n; ++i);	// find lsep

			assert(i!= p);
			if (i> 0&& b[i-1]== r)
				--i;
			//bb= new byte[i];
			System.arraycopy(b, 0, cs.chars, 0, i);
			cs.start= 0;
			cs.end= i;
			if (i== 0)
				System.currentTimeMillis();
			while(++i< p&& (b[i]== r|| b[i]== n));
			System.arraycopy(b, i, b, 0, p- i);
			pos-= i;
			//System.err.println("read "+ i+ " from "+ idx);
			return i;
		}
	}
	
	
	private int fill() {
		
		byte[] b= buf;
		int len= -1;
		synchronized (b) {
			int p= pos;
			len= b.length- p;
			if (maxBytes> 0)
				len= (int) Math.min(len, maxBytes- p);
			try {
				//if (in.available()> 0)
				len= in.read(b, p, len);
			} catch (IOException e) {
				if (Constants.verboseLevel> Constants.VERBOSE_DEBUG)
					System.err.println(e);
				len= -1;
			}
			if (len< 0) {
				try {
					in.close();
				} catch (IOException e) {
					;
				}
				return -1;
				//System.err.println("closed "+ x);
			} else {
				pos+= len;
				//System.err.println("filled "+ len+ " from "+ x);
			}
		}
		//System.err.println("read "+ len);
		return len;
	}
	
	public void close() throws IOException {
		if (in!= null)
			in.close();
	}

	public static final Integer INTEGER_0= new Integer(0);
}
