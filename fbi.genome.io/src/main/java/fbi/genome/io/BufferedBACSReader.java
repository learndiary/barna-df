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

/**
 * A buffered reader class based on <code>ByteArrayCharSequence</code>.
 * @author Micha Sammeth (gmicha@gmail.com)
 * @see ByteArrayCharSequence
 */
public class BufferedBACSReader {
	
	private static final int DEFAULT_BUFFER_SIZE = 1024;
	public static void main(String[] args) {
		
	}
	
	public static byte BYTE_NL= '\n', BYTE_CR= '\r';
	int minVol= 100;
	InputStream in;
	byte[] buf;
	/**
	 * Position up to where buffer is filled
	 */
	int pos= 0;
	/**
	 * Reading position in buffer
	 */
	int lastI= 0;
	long maxBytes= 0l;
	long totBytes= 0l; 
	public long currBytes= 0l;
	long markBytes= -1l;
	
	public BufferedBACSReader(InputStream in) {
		this(in, DEFAULT_BUFFER_SIZE);
	}
	
	/**
	 * @deprecated check whether needed
	 * @param max
	 * @param in
	 */
	public BufferedBACSReader(long max, InputStream in) {
		this(in);
		this.maxBytes= max;
		totBytes= 0;
	}
	
	public BufferedBACSReader(InputStream in, int capacity) {
		this.in= in;
		buf= new byte[capacity];
		
		init();
	}
	
	private void init() {
		pos= 0;
		lastI= 0;
		markBytes= -1;
		currBytes= 0;
	}
	
	/**
	 * Returns the number of bytes that are available in the buffer
	 * before it has to be refilled. If no bytes are left, the number
	 * of bytes after a successful <code>fill()</code> operation is
	 * returned.
	 * @return the number of bytes available in the buffer
	 */
	public int available() {
		
		int diff= pos- lastI;	
		if (diff== 0) {
			if (pos< buf.length) {
				fill();
			} else if (lastI== pos&& lastI> 0) {
				System.arraycopy(buf, lastI, buf, 0, pos- lastI);
				pos-= lastI;
				fill();
				lastI= 0;
			}
			diff= pos- lastI;
		}
		return diff;
	}
	
	/**
	 * Reads from the buffer, fills if needed, and returns a 
	 * <code>ByteArrayCharSequence</code>, potentially reusing
	 * an instance provided.
	 * @param cs <code>ByteArrayCharSequence</code> for object
	 * reuse
	 * @return <code>ByteArrayCharSequence</code> with next line,
	 * or <code>null</code> if the end of the stream is reached
	 */
	public ByteArrayCharSequence readLine(ByteArrayCharSequence cs) {
		
		byte n= BYTE_NL, r= BYTE_CR;
		byte[] b= buf; 
		synchronized (b) {
			
			// find line-sep
			int p= pos;
			int i = lastI;
			for (; i < p&& b[i]!= n; ++i, ++currBytes);	

			// nothing found?
			if (i== p) {
				
				// fill
				if (pos< b.length) 
					fill();
				p= pos;
				if (p== lastI) {
					return null;
				}
				for (; i < p&& b[i]!= n; ++i, ++currBytes);	
				
				// last chance: omit marked position, if any
				if (i== p&& lastI> 0) {
					System.arraycopy(b, lastI, b, 0, p- lastI);
					pos-= lastI;
					fill();
					p= pos;
					i-= lastI;
					lastI= 0;
					for (; i < p&& b[i]!= n; ++i, ++currBytes);
				}
			}
			
			assert(i!= p);
			
			// parse line-sep
			if (i> 0&& b[i-1]== r) {
				--i;
				--currBytes;
			}
			
			// prepare result
			int lineLen= (i-lastI);
			if (cs== null) {
				byte[] b1= new byte[lineLen];
				System.arraycopy(b, lastI, b1, 0, lineLen);
				cs= new ByteArrayCharSequence(b1);
				cs.start= 0;
				cs.end= lineLen;
			} else {
				if (lineLen> cs.chars.length)
					cs.extend(lineLen- cs.chars.length);
				System.arraycopy(b, lastI, cs.chars, 0, lineLen);
				cs.start= 0;
				cs.end= lineLen;
			}
			
			
			// make sure only one lsep is read \n, \r\n, or \n\r
			// ==> enable reading empty lines
			// while(++i< p&& (b[i]== r|| b[i]== n));
			boolean noN= true, noR= true;
			while(i< p) {
				if (b[i]== n&& noN) {
					noN= false;
					++i;
					++currBytes;
				} else if (b[i]== r&& noR) {
					noR= false;
					++i;
					++currBytes;
				} else
					break;
			}
			
			// try to keep mark in buffer
			if (markBytes< currBytes- i) {
				System.arraycopy(b, i, b, 0, p- i);
				pos-= i;
				lastI= 0;
			} else {
				lastI= i;
			}
			
			return cs;
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
	
	public void mark() {
		markBytes= currBytes;
	}
	
	/**
	 * Tries to jump back to the last marked position, iff
	 * this is still within the buffer
	 * @return byte position to which could not be jumped
	 * (if value >= 0)
	 */
	public long reset() {
		
		if (currBytes- pos<= markBytes)  {
			
			long diff= currBytes- markBytes;
			lastI-= diff;
			markBytes= -1;
		} 
		
		return markBytes;
			
	}
	
}
