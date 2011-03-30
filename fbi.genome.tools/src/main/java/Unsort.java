
import commons.ByteArrayCharSequence;
import fbi.genome.io.ThreadedBufferedByteArrayStream;
import fbi.genome.io.UnixStreamSort;
import fbi.genome.io.UnixStreamSort.DesignatedHierarchicalFieldComparator;
import fbi.genome.io.gff.GFFSorter.Cocs;
import fbi.genome.model.constants.Constants;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;


public class Unsort {

	void unsort(File f) throws IOException {
		
		String eol= "\n";
		long size= f.length();
		ByteArrayCharSequence cs= null;
		DesignatedHierarchicalFieldComparator comp =	null; 
			// TODO new UnixStreamSort.RandomComparator(1);	
		

		ThreadedBufferedByteArrayStream buffy = new ThreadedBufferedByteArrayStream(10000, f, true);			
		PipedOutputStream out = new PipedOutputStream();
		PipedInputStream in = new PipedInputStream(out);
		BufferedWriter writer= new BufferedWriter(new OutputStreamWriter(
				out));
		UnixStreamSort sorter = new UnixStreamSort(in, comp);
		
		OutputStream outStr = System.err;
		File outFile = null;
		if (true) {
			outFile = File.createTempFile(f.getName() + "_", "_sorted");
			outStr = new FileOutputStream(outFile);
		}
		Cocs pipe = new Cocs(sorter.getOutInStream(), outStr);
		pipe.setSkipFields(new int[] {});
		pipe.setSepChar("\t");
		pipe.start();
		sorter.start();
		//fieldNrs= new int[]{fieldNrs[0], fieldNrs[1]};
		long bytesRead= 0;
		int lastPerc= 0, rowCtr= 0;
		for (ByteArrayCharSequence line = buffy.readLine(cs); line.end!= 0; line = buffy
				.readLine(cs)) {
			++rowCtr;
			bytesRead += line.length() + eol.length();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				int perc = (int) ((bytesRead * 10d) / size);
				if (perc > lastPerc) {
					++lastPerc;
					if (Constants.progress!= null)
						Constants.progress.progress();
					else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						System.err.print("*");
						System.err.flush();
					}
				}
			}
							
			writer.write(line.toString());
			writer.write(eol);
		}
		System.gc();
		Thread.currentThread().yield();
		writer.flush();
		writer.close();

		try {
			pipe.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
}
