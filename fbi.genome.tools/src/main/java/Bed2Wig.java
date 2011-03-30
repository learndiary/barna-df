import commons.ByteArrayCharSequence;
import fbi.genome.io.BufferedByteArrayReader;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.commons.IntVector;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;


public class Bed2Wig {

	public static void main(String[] args) {
		
		if (args== null|| args.length!= 2) {
			System.err.println("Usage: Bed2Wig [input.bed] [output.wig]");
			System.exit(-1);
		}
		
		
		File bedFile= new File(args[0]);
		File wigFile= new File(args[1]);
		if (!bedFile.exists()) {
			System.err.println(args[0]+ " not found.");
			System.exit(-1);
		}
		
		bed2wig(bedFile, wigFile);
	}
	static void bed2wig(File in, File out) {
		try {
			//BufferedReader buffy= new BufferedReader(new FileReader(in));
			BufferedByteArrayReader buffy= new BufferedByteArrayReader();
			FileInputStream istream= new FileInputStream(in);
			ByteArrayCharSequence cs= new ByteArrayCharSequence(300);
			IntVector v= new IntVector(10000, 10000);
			int start= -1;
			int ctr= 0;
			String chr= null;
			while(buffy.readLine(istream, cs)!= null) {
				++ctr;
				if (ctr% 1000== 0) {
					System.err.print("*");
					System.err.flush();
				}
				BEDobject2 obj= new BEDobject2(cs);
				if (start== -1) {
					start= obj.getStart();
					chr= obj.getChr().toString();
				}
				if (obj.getBlockCount()< 2) {
					for (int i = obj.getStart(); i < obj.getEnd(); i++) {
						v.set(i- start, v.get(i- start)+ 1);
					}
				} else {
					int bedstart= obj.getStart();
					for (int i = 0; i < obj.getBlockCount(); i++) {
						int bsize= obj.getNextBlockSize();
						int bstart= obj.getNextBlockStart();
						for (int j = 0; j < bsize; j++) {
							int gpos= bedstart+ bstart+ j;
							v.set(gpos- start, v.get(gpos- start)+ 1);
						}
					}
				}
			}
			istream.close();
			
			// write
			BufferedWriter writer= new BufferedWriter(new FileWriter(out));
			++start;
			writer.write("track type=wiggle_0 name="+in.getName()+"\n");
			writer.write("fixedStep chrom="+chr+" start="+start+" step=1\n");
			for (int i = 0; i < v.size(); i++) {
				writer.write(Integer.toString(v.get(i)));
				writer.write("\n");
			}
			writer.flush();
			writer.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
