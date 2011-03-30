import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.HashMap;

public class SequenceRetriever {
	
	static final HashMap<Character, Character> mapRC= new HashMap<Character, Character>(22);
	static{
		mapRC.put('a', 't');
		mapRC.put('c', 'g');
		mapRC.put('g', 'c');
		mapRC.put('t', 'a');
		mapRC.put('A', 'T');
		mapRC.put('C', 'G');
		mapRC.put('G', 'C');
		mapRC.put('T', 'A');
		
		mapRC.put('n', 'n');
		mapRC.put('N', 'N');
		mapRC.put('r', 'y');
		mapRC.put('R', 'Y');
		mapRC.put('y', 'r');
		mapRC.put('Y', 'R');
		mapRC.put('s', 's');
		mapRC.put('S', 'S');
		mapRC.put('w', 'w');
		mapRC.put('W', 'W');
		mapRC.put('k', 'm');
		mapRC.put('K', 'M');
		mapRC.put('m', 'k');
		mapRC.put('M', 'K');
	}
	
	public static void getSequenceFromGTF(File fileGTF, File fastaDir, PrintStream p) {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fileGTF)), bufasta= null;
			StringBuffer sb= new StringBuffer();
			FileInputStream istream= null;
			String lchr= null, lstr= null, lline= null;
			int pos= -1, lineLen= -1;
			for (String s= null; (s= buffy.readLine())!= null; ) {
				String[] ss= s.split("\\s");
				int start= Integer.parseInt(ss[3]), end= Integer.parseInt(ss[4]); // incl.
				sb.ensureCapacity(end- start+ 1);
				if ((!ss[0].equals(lchr))|| (!ss[6].equals(lstr))|| start< pos) {
					if (istream!= null) {
						bufasta.close();
						istream.close();
					}
					istream= new FileInputStream(new File(fastaDir.getAbsolutePath()+ File.separator+ ss[0]+ ".fa"));
					bufasta= new BufferedReader(new InputStreamReader(istream));
					bufasta.readLine(); // >
					lline= bufasta.readLine();
					lineLen= lline.length();
					pos= 1;	// 1-based
					lchr= ss[0];
					lstr= ss[6];					
				}

				// skip to start line
				int skipLines= (start- pos)/ lineLen;
				if (skipLines> 0) {
					bufasta.skip((skipLines- 1)* (lineLen+ 1)); // FS len
					pos+= skipLines* lineLen;
					lline= bufasta.readLine();
				}
				// first line
				sb.append(lline.substring(start- pos, Math.min(end- pos+ 1, lineLen)));
				// intermediate
				while(pos+ lineLen<= end) {
					lline= bufasta.readLine();
					pos+= lineLen;
					sb.append(lline.substring(0, Math.min(end- pos+ 1, lineLen)));
				}
				
				if (ss[6].equals("-")) 
					reverseComplement(sb);
				
				p.println(">"+ ss[0]+ ":"+ ss[3]+ "-"+ ss[4]+ ""+ ss[6]);
				p.println(sb.toString());
				sb.setLength(0);
				//p.println(ss[0]+"\t"+ss[3]+"\t"+ss[4]+"\t"+seq.length()+" "+(end- start+ 1)+"\n"+seq);	
				//System.currentTimeMillis();
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	static void reverseComplement(StringBuffer sb) {
		int len= sb.length(), half= len/ 2;		
		for (int i = 0; i < half; ++i) {
			int i2= len- 1- i;
			char c1= sb.charAt(i), c2= sb.charAt(i2);
//			if (Character.isLowerCase(c1)|| Character.isLowerCase(c2)) {
//				char d1= mapRC.get(c1), d2= mapRC.get(c2);
//				System.currentTimeMillis();
//			}
			sb.setCharAt(i, mapRC.get(c2));
			sb.setCharAt(i2, mapRC.get(c1));
		}
		if (len% 2!= 0) 
			sb.setCharAt(half, mapRC.get(sb.charAt(half)));		
	}
	
	public static void main(String[] args) {
		try {
			getSequenceFromGTF(new File("/Users/micha/annotation/hg19_RefSeq_fromUCSC100615_introns_uniq.gtf"), new File("/genomes/hg19"), 
					new PrintStream(new File("/Users/micha/annotation/hg19_RefSeq_fromUCSC100615_introns_uniq.mfasta")));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
}
