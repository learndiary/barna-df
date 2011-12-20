package barna.astalavista;

public class Annotator extends Thread {

	public static final String[] GENE_LOCALIZATION= new String[] {"cdsRef", "cdsMax", "cdsAll", "cdsEvent"};  
	public static final int  GENE_LOC_NO= -1, GENE_LOC_REFCDS= 0, GENE_LOC_MAXCDS= 1, GENE_LOC_ALLCDS= 2, GENE_LOC_EVENTCDS=3;
	
	int geneLoc= GENE_LOC_REFCDS;
	
	
	static void printUsage() {
		System.err.println("Howdy, I am the AStalavista Event Annotator.\n");
		System.err.println("I understand:\n\tAnnotator event.gtf annotation.gtf [options]\n");
		System.err.println("where [options] can be one or more of the following things:\n");
		
		// cds
		System.err.print("\t[");
		for (int i = 0; i < GENE_LOCALIZATION.length; i++) {
			System.err.print("-"+ GENE_LOCALIZATION[i]);
			if (i< GENE_LOCALIZATION.length- 1)
				System.err.print("|");
			else
				System.err.print("]\t");
		}
		System.err.println("annotate the cds");
	}
	
	
	
}
