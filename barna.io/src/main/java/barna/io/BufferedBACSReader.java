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

package barna.io;

import barna.commons.ByteArrayCharSequence;

import java.io.BufferedReader;
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
			// note: i== lastI can hold true, 
			// eg empty lines \n\n produce an
			// empty string
			if (p== lastI) 
				return null;	// nothing more, EOS
			
			// parse line-sep
			if (i> lastI&& b[i-1]== r) {
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
				cs.resetFind();
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
		
		int len= -1;
		byte[] b= buf;
		synchronized (b) {
			int p= pos;
			len= b.length- p;
			if (maxBytes> 0)
				len= (int) Math.min(len, maxBytes- p);
			try {
				if (in.available()> 0)
					len= in.read(b, p, len);
				else
					return 0;	// nothing read, maybe EOS reached
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			if (len< 0) {
				try {
					in.close();
				} catch (IOException e) {
					;
				}
				return -1;
			} else {
				pos+= len;
			}
		}
		
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
	 * (if value &ge; 0)
	 */
	public long reset() {
		
		if (currBytes- pos<= markBytes)  {
			
			long diff= currBytes- markBytes;
			lastI-= diff;
			markBytes= -1;
			currBytes-= diff;
		} 
		
		return markBytes;
			
	}
	
}
