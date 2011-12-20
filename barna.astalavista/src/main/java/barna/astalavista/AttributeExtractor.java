package barna.astalavista;

import org.apache.commons.cli.*;

public class AttributeExtractor {
    public static final byte VERBOSE_NORMAL = 1;
    public static final byte VERBOSE_SHUTUP = 0;

    public static byte verboseLevel= VERBOSE_NORMAL;

	private static Options options;
	private static Options getOptions() {
		if (options == null) {
			options = new Options();
			options.addOption("c", "coding", false, "extract coding transcripts");
			options.addOption("n", "noncoding", false, "extract non-coding transcripts");
		}

		return options;
	}
	
	public static void main(String[] args) {
		Parser parser= new PosixParser();
		CommandLine line= null;
		try {
			line= parser.parse(getOptions(), args);
		} catch (ParseException e) {
			if (verboseLevel> VERBOSE_SHUTUP)
				System.err.println("[UPS] "+e.getMessage());
			System.exit(-1);
		}
		
		AttributeExtractor myExtractor= new AttributeExtractor();
		myExtractor.init(line);
		myExtractor.run();
	}
	
	public void init(CommandLine line) {
		; //
	}
	
	public void run() {
		; //
	}
	
	
}
