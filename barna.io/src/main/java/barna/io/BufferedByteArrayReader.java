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
import barna.model.constants.Constants;

import java.io.IOException;
import java.io.InputStream;

/**
 * @deprecated marked for deletion
 * @see BufferedBACSReader
 * @author micha
 *
 */
public class BufferedByteArrayReader {

	static final byte BYTE_NL= '\n', BYTE_CR= '\r';
	byte[] b;
	int pos;
	
	public BufferedByteArrayReader(int size) {
		b= new byte[size];
		pos= 0;
	}
	public BufferedByteArrayReader() {
		this(1024);
	}
	
	public ByteArrayCharSequence readLine(InputStream idx, ByteArrayCharSequence cs) {
			byte n= BYTE_NL, r= BYTE_CR;
			int i= 0;
			for (; i < pos&& b[i]!= n; ++i);	// find lsep
			while (i== pos) {
				int nu= 0;
				try {
					nu= idx.read(b, pos, b.length- pos);
				} catch (IOException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
						e.printStackTrace();
					return null;
				}
				if (nu<= 0)
					break;
				pos+= nu;
				for (; i < pos&& b[i]!= n; ++i);
				if (i== pos) {
					byte[] c= new byte[b.length* 2];
					System.arraycopy(b, 0, c, 0, b.length);
					b= c;
				}
			}
			if (i== pos)
				return null;
			if (i> 0&& b[i-1]== r)
				--i;
			System.arraycopy(b, 0, cs.chars, 0, i);
			cs.start= 0;
			cs.end= i;
			while(++i< pos&& (b[i]== r|| b[i]== n));
			System.arraycopy(b, i, b, 0, pos- i);
			pos-= i;
			return cs;
		}

}
