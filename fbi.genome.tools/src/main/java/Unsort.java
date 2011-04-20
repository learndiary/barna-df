import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.tools.UnixStreamSort;
import fbi.commons.tools.UnixStreamSort.DesignatedHierarchicalFieldComparator;
import fbi.genome.io.ThreadedBufferedByteArrayStream;
import fbi.genome.io.gff.GFFSorter.Cocs;

import java.io.*;


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
            Log.progress(bytesRead, size);
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
