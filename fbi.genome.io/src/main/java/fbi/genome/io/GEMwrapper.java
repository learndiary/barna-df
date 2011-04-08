package fbi.genome.io;

import java.io.File;

public class GEMwrapper {
	public static void main(String[] args) {
		GEMobject.GEMtoBED(
				new File("C:\\workspace\\Genome\\resources\\formats\\GEM_alignment_qualities.0.map"), 
				null, 
				null);
	}
}
