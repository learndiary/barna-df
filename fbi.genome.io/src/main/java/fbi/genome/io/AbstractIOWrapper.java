package fbi.genome.io;

import java.io.File;
import java.io.FileOutputStream;

public abstract class AbstractIOWrapper implements IOWrapper {

	/**
	 * Sorts and writes the output to the file provided
	 * @param outputFile handle to which the output is
	 * sent
	 */
	public void sort(File outputFile) {
		try {
			FileOutputStream fos= new FileOutputStream(outputFile);
			sort(fos);
			fos.close();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

}
