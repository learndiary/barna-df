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
