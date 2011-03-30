package fbi.genome.io;

import commons.ByteArrayCharSequence;

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
