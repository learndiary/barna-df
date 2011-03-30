package fbi.genome.astalavista;

import java.io.File;

import fbi.genome.io.SpliceGraphIO;
import fbi.genome.model.IntronModel;
import fbi.genome.model.splicegraph.Graph;

public class SJextractor {

	
	public static void main(String[] args) {
		
		int usFlank= 0, dsFlank= 0;
		IntronModel iModel= null;
		File inFile= null, outFile= null;
		try {
			usFlank= Integer.parseInt(args[0]);
			dsFlank= Integer.parseInt(args[1]);
			for (int i = 2; i < args.length; i++) {
				if (args[i].equalsIgnoreCase("--iModel")) {
					iModel= new IntronModel();
					iModel.read(new File(args[++i]));
				} else if (args[i].equalsIgnoreCase("--ann")) {
					inFile= new File(args[++i]);
				} else if (args[i].equalsIgnoreCase("--genome")) {
					fbi.genome.model.Graph.overrideSequenceDirPath= args[++i];
				} else if (args[i].equalsIgnoreCase("--out")) {
					outFile= new File(args[++i]);
				}
			}
		} catch (Exception e) {
			System.err.println("[NONO] Use this: asta -c extractSJ usFlankLen dsFlankLen --ann <file> --genome <path> --iModel <file> [--out <file>]");
			System.err.println("\twhere");
			System.err.println("usLength\tis the length of the upstream flank of the splice junction");
			System.err.println("dsLength\tis the length of the downstream flank of the splice junction");
			System.err.println("--ann\tthe annotation file (GTF format)");
			System.err.println("--genome\tis the path to the genomic sequence (i.e., the chromFa directory in UCSC convention)");
			System.err.println("--iModel\ta file with the intron model");
			System.err.println("--out\tspecifies an optional output file (default: stdout)");
			System.exit(-1);
		}
		
		SpliceGraphIO.extractSpliceJunctions(usFlank, dsFlank, iModel, inFile, outFile);
	}
}
