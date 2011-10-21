package fbi.genome.io;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

public abstract class AbstractIOWrapper implements IOWrapper {

	/**
	 * Sorts and writes the output to the file provided
	 * @param outputFile handle to which the output is
	 * sent
	 */
	public void sort(File outputFile) {
		
		FileOutputStream fos= null;
		try {
			fos= new FileOutputStream(outputFile);
			sort(fos);
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			try {
				if (fos!= null)
					fos.close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
	}

}
