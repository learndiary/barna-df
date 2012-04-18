import java.io.File;
import java.io.FileInputStream;

import barna.commons.ByteArrayCharSequence;
import barna.io.BufferedBACSReader;
import barna.model.Graph;
import barna.model.bed.BEDobject2;



public class Bed2Freq {

	public static void main(String[] args) {
		File genDir= new File("");
		File bedFile= new File("");
		//long[][] freqs= getFreq(bedFile,0,10);
		//printFreq(freqs);
	}

	protected static long[][] getFreq(File bedFile, int fivePrime, int threePrime) throws Exception {
		
		int mLen= (-fivePrime)+ threePrime;
		long[][] f= new long[4][];
		for (int k = 0; k < f.length; k++) {
			f[k]= new long[mLen];
			for (int j = 0; j < f[k].length; j++) 
				f[k][j]= 0;
		}
		
		BEDobject2 bed= null;
		ByteArrayCharSequence cs= new ByteArrayCharSequence(150),
			buf= new ByteArrayCharSequence(mLen);
		BufferedBACSReader reader= new BufferedBACSReader(new FileInputStream(bedFile));
		while (reader.readLine(cs)!= null) {
			bed= new BEDobject2(cs);
			if (bed.getBlockCount()> 1)
				continue;
			Graph.readSequence(bed.getChr(), 
					bed.getStrand()>= 0, 
					bed.getStart()+ 1, 
					bed.getEnd(), 
					buf, 
					0, mLen);
			buf.toUpperCase();
			
		}
		
		return null;
	}
}
