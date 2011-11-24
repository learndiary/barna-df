import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;


public class ReadExtractor {

	public static void main(String[] args) {
		if (args== null|| args.length< 2) {
			System.err.println("[USAGE] ReadExtractor <file> <reads.bed> [tolerance]");
			System.err.println("\nwhere\n");
			System.err.println("<file> is the name of the chromosome");
			System.err.println("and contains lines of the form \"start\\tend\"");
			System.err.println("\nand\n");
			System.err.println("<reads.bed> is the read file (BED format)");
			System.exit(0);
		}
		
		try {
			System.err.println("Hi!");
			int tol= 0;
			if (args.length> 2) {
				try {
					tol= Integer.parseInt(args[2]);
					System.err.println("tolerance= "+ tol);
				} catch (Exception e) {
					System.err.println("\tnot an integer "+ args[2]);
				}
			}
			System.err.println("\n[LOCI] getting the loci");
			HashMap<Integer, Integer> map= getMap(args[0]);
			int[] starts= new int[map.size()], ends= new int[starts.length];
			System.err.println("\tfound "+starts.length+" loci.");
			Iterator<Integer> iter= map.keySet().iterator();
			int p= 0; 
			while(iter.hasNext()) 
				starts[p++]= iter.next();
			Arrays.sort(starts);
			for (int i = 0; i < starts.length; i++) { 
				ends[i]= map.get(starts[i])+ tol;
				starts[i]-= tol;
			}
			System.err.println("\n[EXTRACT] getting the reads");
			int nr= extractReads(args[0], args[1], starts, ends);
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}

	private static HashMap<Integer, Integer> getMap(String fname) throws Exception {
		HashMap<Integer, Integer> map= new HashMap<Integer, Integer>(1000);
		BufferedReader buf= new BufferedReader(new FileReader(fname));
		String s= null;
		while ((s= buf.readLine())!= null) {
			String[] ss= s.split("\t");
			int start= Integer.parseInt(ss[0]), end= Integer.parseInt(ss[1]);
			if (start> end) {
				int h= start; start= end; end= h;
			}
			map.put(start, end);
		}
		buf.close();
		return map;
	}

	private static int extractReads(String chrName, String bedFname, int[] starts, int[] ends) throws Exception {
		
		chrName= chrName.substring(chrName.lastIndexOf(File.separator)+ 1);
		BufferedReader buf= new BufferedReader(new FileReader(bedFname));
		int p= bedFname.lastIndexOf('.');
		String outName= bedFname;
		if (p>= 0)
			outName= outName.substring(0, p);
		outName+= "_set.bed";
		BufferedWriter wri= new BufferedWriter(new FileWriter(outName));
		String line= null;
		while ((line= buf.readLine())!= null) {
			if (line.startsWith(chrName))
				break;
		}
		if (line== null)
			return 0;
		wri.write(line);
		wri.write("\n");
		int nr= 1, tot= 1;
		long t0= System.currentTimeMillis();		
		while ((line= buf.readLine())!= null&& line.startsWith(chrName)) {
			++tot;
			String[] ss= line.split("\t");
			int start= Integer.parseInt(ss[1]);
			p= Arrays.binarySearch(starts, start);
			if (p< 0)
				p= -(p+ 1)- 1;
			if (p< 0|| p>= starts.length) {
				int end= Integer.parseInt(ss[2]);
				p= Arrays.binarySearch(starts, end);
				if (p< 0)
					p= -(p+ 1)- 1;
				if (p< 0|| p>= starts.length|| end> ends[p])
					continue;
			} else if (start> ends[p])
				continue;
			
			wri.write(line);
			wri.write("\n");
			++nr;
		}
		wri.flush();
		wri.close();
		buf.close();
		System.err.println("\textracted "+ nr+ " of "+ tot+" reads ("
				+ (System.currentTimeMillis()- t0)/1000+" sec.)");
		return nr;
	}
}
