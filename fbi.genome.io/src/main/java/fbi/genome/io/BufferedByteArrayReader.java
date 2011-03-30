package fbi.genome.io;

import commons.ByteArrayCharSequence;
import fbi.genome.model.constants.Constants;

import java.io.IOException;
import java.io.InputStream;

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
			System.arraycopy(b, 0, cs.a, 0, i);
			cs.start= 0;
			cs.end= i;
			while(++i< pos&& (b[i]== r|| b[i]== n));
			System.arraycopy(b, i, b, 0, pos- i);
			pos-= i;
			return cs;
		}

}
