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

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.Hashtable;

public class ColorSpacer {

	// mapping["0"] = {"T":"T","A":"A","C":"C","G":"G"}
    // mapping["1"] = {"T":"G","A":"C","C":"A","G":"T"}
    // mapping["2"] = {"T":"C","A":"G","C":"T","G":"A"}
    // mapping["3"] = {"T":"A","A":"T","C":"G","G":"C"}
	
	public static Hashtable<Character, char[]> map2char= new Hashtable<Character, char[]>(4,1f);
	static {
		map2char.put('T', new char[]{'T', 'G', 'C', 'A'});
		map2char.put('A', new char[]{'A', 'C', 'G', 'T'});
		map2char.put('C', new char[]{'C', 'A', 'T', 'G'});
		map2char.put('G', new char[]{'G', 'T', 'A', 'C'});
	}
	
	public static void main(String[] args) {
		
	}
	
	public static void encodeFASTA(File fin, File fout) {
		try {
			ThreadedBufferedByteArrayStream tbbs= new ThreadedBufferedByteArrayStream(10000, fin, true);
			BufferedOutputStream bossy= new BufferedOutputStream(new FileOutputStream(fout));
			ByteArrayCharSequence cs= new ByteArrayCharSequence(10000);
			tbbs.readLine(cs);
			char c= ' ';
			int out= 0;
			byte[] outBuf= new byte[10000];
			while (tbbs.readLine(cs)!= null) {
				if (cs.charAt(0)== '>') {
					c= ' ';
					out= 0;
					continue;
				}
				
				int start= 0;
				if (c== ' ') 
					c= cs.charAt(start++);
				for (int i= start; i < cs.length(); ++i, ++out) {
					outBuf[out]= cs.byteAt(i);
				}
				bossy.write(outBuf, 0, out);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
