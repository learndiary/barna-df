import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;

import fbi.commons.ByteArrayCharSequence;
import fbi.genome.io.BufferedByteArrayReader;
import fbi.genome.io.FileHelper;
import fbi.genome.io.bed.BEDwrapper;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.commons.IntVector;


public class Bed2Wig {

	static void myMain() {
		File dir= new File("/Users/micha/projects/demassy/download_new/regions");
		String[]  ss= dir.list();
		for (int i = 0; i < ss.length; i++) {
			if (!ss[i].endsWith("_reads.sorted.bed"))
				continue;
			File fin= new File(dir.getAbsolutePath()+ File.separator+ ss[i]);
			File fout= FileHelper.replaceSfx(fin, ".bgraph");
			
			System.err.println(fin.getName());
			bed2wig(fin, fout);
		}
	}
	
	public static void main(String[] args) {
		
//		if (1== 1)
//			myMain();
		
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
	
	/**
	 * assumes sorted by (1) chr, (2) start position
	 * @param in
	 * @param out
	 */
	static void bed2wig(File in, File out) {
		try {
			
			BEDwrapper bedreader= new BEDwrapper(in);
			if (!bedreader.isApplicable()) {
				File tmp= FileHelper.createTempFile(Bed2Wig.class.getName(), in.getName());
				tmp.deleteOnExit();
				bedreader.sort(tmp);
				in= tmp;
			}
			
			//BufferedReader buffy= new BufferedReader(new FileReader(in));
			if (out.exists())
				out.delete();
			BufferedByteArrayReader buffy= new BufferedByteArrayReader();
			FileInputStream istream= new FileInputStream(in);
			ByteArrayCharSequence cs= new ByteArrayCharSequence(300);
			IntVector v= new IntVector(10000, 10000);
			int start= -1;
			int ctr= 0;
			String chr= null, lastChr= null;			
			while(buffy.readLine(istream, cs)!= null) {
				++ctr;
				if (ctr% 1000== 0) {
					System.err.print("*");
					System.err.flush();
				}
				BEDobject2 obj= new BEDobject2(cs);
				chr= obj.getChr().toString();
				
				if (!chr.equals(lastChr)) {
					writeBedGraph(in.getName(), lastChr, start, v, out);
					v.removeAll();
					start= obj.getStart();
					lastChr= chr;
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
							if (gpos- start< 0)
								System.currentTimeMillis();
							v.set(gpos- start, v.get(gpos- start)+ 1);
						}
					}
				}
			}
			istream.close();

			writeBedGraph(in.getName(), lastChr, start, v, out);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void writeFixedStep(String name, String chr, int start, IntVector v, File out) {
		
		// write
		try {
			boolean header= !out.exists();
			BufferedWriter writer= new BufferedWriter(new FileWriter(out, true));
			++start;
			if (header) {
				writer.write("track type=wiggle_0 name="+ name+ "\n");
				writer.write("fixedStep chrom="+chr+" start="+start+" step=1\n");
			}
			for (int i = 0; i < v.size(); i++) {
				writer.write(Integer.toString(v.get(i)));
				writer.write("\n");
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	/**
	 * @deprecated not implemented
	 * @param name
	 * @param chr
	 * @param start
	 * @param v
	 * @param out
	 */
	private static void writeVariableStep(String name, String chr, int start, IntVector v, File out) {
		
		// write
		try {
			boolean header= !out.exists();
			BufferedWriter writer= new BufferedWriter(new FileWriter(out, true));
			++start;
			if (header) {
				writer.write("track type=wiggle_0 name="+ name+ "\n");
				writer.write("variableStep chrom="+chr+" start="+start+" step=1\n");
			}
			for (int i = 0; i < v.size(); i++) {
				writer.write(Integer.toString(v.get(i)));
				writer.write("\n");
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private static void writeBedGraph(String name, String chr, int start, IntVector v, File out) {
		
		// write
		try {
			boolean header= !out.exists();
			BufferedWriter writer= new BufferedWriter(new FileWriter(out, true));
			++start;
			if (header) {
				writer.write("track type=bedGraph name=\""+ name+ "\" description=\""+ name+"\" " +
						"visibility=dense color=200,100,0 altColor=0,100,200 priority=20\n");
			}
			for (int i = 0; i < v.size(); i++) {
				int val= v.get(i);
				int j = i+1;
				for (; j < v.size(); j++) 
					if (v.get(j)!= val)
						break;
				
				if (val!= 0) 
					writer.write(chr+ " "+ (start+ i)+ " "+ (start+ j)+" "+ Integer.toString(val)+ "\n");
				
				i= (j-1);
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}


}
