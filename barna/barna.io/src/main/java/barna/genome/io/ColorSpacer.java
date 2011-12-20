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

package barna.genome.io;

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
