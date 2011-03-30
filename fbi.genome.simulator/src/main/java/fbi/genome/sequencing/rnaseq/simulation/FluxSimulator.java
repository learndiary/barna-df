package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.Log;
import fbi.genome.io.SpliceGraphIO;
import fbi.genome.model.IntronModel;
import fbi.genome.model.constants.Constants;
import org.apache.commons.math.random.RandomDataImpl;

import java.io.*;
import java.util.Properties;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.concurrent.Callable;

//import gphase.solexa.lp.GraphLPsolver5;
//import gphase.solexa.simulation.Nebulizer;

// TODO check whether selected fragments concord w pro file

public class FluxSimulator implements Callable<Void> {
    /**
     * Default flank length
     */
    private static final int FLANK5_LEN= 50, FLANK3_LEN= 50;
	public static String version= "";

	public static boolean c= false;

	public static final String PROPERTY_BUILD= "simulator.build";
	public static final String PROPERTY_JDK= "simulator.jdk";
	
	//boolean gui= false, help= false, doExpr= false, doLib= false, doSeq= false, doSJ= false, doInstall= false;
	boolean help= false, doExpr= false, doLib= false, doSeq= false, doSJ= false, doInstall= false;
	File file, fileIM, fileGenome;
	FluxSimulatorSettings settings;
	Profiler profiler;
	Fragmenter fragmenter;
	Sequencer sequencer;
	int[] eflanks= null;



	public static void exit(int code) {
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			
		}
	}
	
	private static final String FNAME_PROPERTIES= "simulator.prop";

    // todo reads internal properties file to get simulator version and check java version

    /**
     * Read in properties and return true if properties are valid and everything is fine.
     * Return false otherwise.
     *
     * Currently checks for a program version and checks the current java version. The latter is invalid
     * if it is not new enough.
     *
     * @return valid returns true if properties are valid and everything is fine
     */
	static boolean readProperties() {
		String wrapper= System.getProperty(Constants.PROPERTY_KEY_WRAPPER_BASE);
		if (wrapper== null)
			wrapper= FNAME_PROPERTIES;
		else {
			if (wrapper.startsWith("\'"))	// win double-quoting in bat
				wrapper= wrapper.substring(1);
			if (wrapper.endsWith("\'"))
				wrapper= wrapper.substring(0, wrapper.length()- 1);
			wrapper+= File.separator+ FNAME_PROPERTIES;
		}
		File pFile= new File(wrapper);		
		Properties props= new Properties();
		try {
			props.load(new FileInputStream(pFile));
		} catch (Exception e) {
			; // :)
		}

		if (props.containsKey(PROPERTY_BUILD))
			version= (String) props.get(PROPERTY_BUILD);
		if (props.containsKey(PROPERTY_JDK)) {
			try {
				float v= Float.parseFloat((String) props.get(PROPERTY_JDK));
				String ver= System.getProperty("java.version");
				int p= ver.indexOf('.', 0);
				p= ver.indexOf('.', p+1);
				float v2= Float.parseFloat(ver.substring(0, p));
				if (v2< v) {
					Log.error("[BEHIND] Wrong java version, I need "+v+" but I found "+v2+".");
                    return false;
				}
			} catch (Exception e) {
				; // :)
			}
			
		}
        return true;
		
	}
	
	public static void main(String[] args) {
		
		if(!readProperties()){
            System.exit(-1);
        }
		// prepare options
        FluxOptions options = new FluxOptions();
        try {
            if(!options.parseParameters(args)){
                System.exit(-1);
            }
        } catch (Exception e) {
            Log.error("Error while parsing command line parameters: " + e.getMessage());
            System.exit(-1);
        }


        FluxSimulator sim= new FluxSimulator();
        options.apply(sim);
        try {
            sim.call();
        } catch (Exception e) {
            Log.error(e.getMessage());
        }
    }
	
	
	public static boolean invertTable(File invFile) {
		
		System.err.println("[INFO] inverting .pro file");
		File tmpFile= new File(invFile.getAbsolutePath()+"_inv");
		try {
			Vector<StringTokenizer> lineTokis= new Vector<StringTokenizer>();
			BufferedReader buffy= new BufferedReader(new FileReader(invFile));
			for (String s= buffy.readLine(); s!= null; s= buffy.readLine()) 
				lineTokis.add(new StringTokenizer(s));
			buffy.close();
			
			int c= -1;
			for (int i = 0; i < lineTokis.size(); i++) {
				if (c< 0)
					c= lineTokis.elementAt(i).countTokens();
				else if (c!= lineTokis.elementAt(i).countTokens()) {
					System.err.println("\t[OHNO] inconsistent column count "+c+" <> "+lineTokis.elementAt(i).countTokens());
					return false;
				}
			}
			
			BufferedWriter writer= new BufferedWriter(new FileWriter(tmpFile));
			while (true) {
				for (int i = 0; i < lineTokis.size(); i++) {
					writer.write(lineTokis.elementAt(i).nextToken());
					if (i< lineTokis.size()-1)
						writer.write("\t");
					else
						writer.write("\n");
				}
				if (!lineTokis.elementAt(0).hasMoreTokens()) 
					break;					
			}			
			writer.flush();
			writer.close();

		} catch (Exception e) {
			e.printStackTrace();
			if (tmpFile.exists())
				tmpFile.delete();
			return false;
		}
		
		if (!invFile.delete()) 
			System.err.println("\t[OHNO] failed to remove "+invFile);
		if (!tmpFile.renameTo(invFile)) 
			System.err.println("\t[OHNO] failed to move "+tmpFile+" to "+invFile);
		return true;
	}
	
	public FluxSimulator() {

        // disabled GUI here
		//if ((args== null)|| (args.length== 0)) {
		//	gui= true;
		//}

	}





	public void setFile(String fName) {
		file= new File(fName);
	}
	public void setImodel(String fName) {
		fileIM= new File(fName);
	}
	
	public void setDoSJ(String fName) {
		doSJ= true;
		file= new File(fName);
	}
	public void setDoInstall(){
		doInstall= true;
	}
	public void setDoExpr(){
		doExpr= true;
	}
	public void setDoLib(){
		doLib= true;
	}
	public void setDoSeq(){
		doSeq= true;
	}
	public void setHelp(){
		help= true;
	}
	
	public int[] getEFlanks() {
		if (eflanks == null) {
			eflanks = new int[2];
			eflanks[0]= -1;
			eflanks[1]= -1;
		}

		return eflanks;
	}
	public void setGenomeDir(String path) {
		fbi.genome.model.Graph.overrideSequenceDirPath= path;
		this.fileGenome= new File(path);
	}
	public void set5flank(String length) {
		int x= -1;
		try {
			x= Integer.parseInt(length);
			if (x<= 0)
				throw new RuntimeException("Illegal Value");
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[NOINT] Not a valid length for 5' exon flank: "+length);
		}
		getEFlanks()[0]= x;
	}
	public void set3flank(String length) {
		int x= -1;
		try {
			x= Integer.parseInt(length);
			if (x<= 0)
				throw new RuntimeException("Illegal Value");
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[NOINT] Not a valid length for 3' exon flank: "+length);
		}
		getEFlanks()[1]= x;
	}
	
	

	
	public Void call() throws Exception {
        Log.info("[HELLO] I am the Flux Simulator (build "+version+"), nice to meet you!\n");
        if (file== null || !file.exists()) {
            throw new RuntimeException("I have no parameter file and I want to scream!");
        }

		if (doSJ) {
			if (fileGenome== null|| !fileGenome.exists()) {
				throw new RuntimeException("[AIAIII] I have no directory with the genomic sequences!");
			} else{
                // todo : refoactor this to not use a static variable
				fbi.genome.model.Graph.overrideSequenceDirPath= fileGenome.getAbsolutePath();
            }

			int[] flanks= getEFlanks();
			flanks[0]= flanks[0]< 0? FLANK5_LEN: flanks[0];
			flanks[1]= flanks[1]< 0? FLANK3_LEN: flanks[1];

			IntronModel iModel= new IntronModel();
			if (fileIM!= null){
                Log.info("Reading Model file " + fileIM.getAbsolutePath());
				iModel.read(fileIM);
            }
            Log.info("Extracting splice junctions, 5'sequence "+flanks[0]
                + ", 3'sequence "+eflanks[1]+", intron model "+ iModel.getName());
            SpliceGraphIO.extractSpliceJunctions(flanks[0], flanks[1], iModel, file, null);
            // todo: exit here ?
		}


		// init
		settings= FluxSimulatorSettings.createSettings(file);
		if (settings== null)
			throw new RuntimeException("Unable to create Setting file!");

		Log.info("[INIT] I am collecting information on the run.");
		profiler= new Profiler(settings);
		settings.setProfiler(profiler);
		if (settings.getProFile()!= null&& settings.getProFile().exists()&& settings.getProFile().canRead()) {
			profiler.loadStats();
			if (doExpr&& profiler.isFinishedExpression()) {
				if (!userCLIconfirm("[CAUTION] I overwrite the expression values in file "+settings.getProFile().getName()+", please confirm:\n\t(Yes,No,Don't know)"))
					return null;
				else {
					boolean b= settings.getProFile().delete();	// TODO maybe only remove rfreqs..
					if (settings.getProfiler()!= null)
						settings.getProfiler().status= -1;
				}
			}
            Log.info("");
		}
		fragmenter= new Fragmenter(settings);
		if (settings.getFrgFile()!= null&& settings.getFrgFile().exists()&& settings.getFrgFile().canRead()) {
			if (doExpr|| doLib) {
				//System.err.println("RE-USING library");
				// see doLib()
//				if (!userCLIconfirm("[WARNING] I will overwrite the library file "+settings.getFrgFile().getName()+", please confirm:\n\t(Yes,No,Don't know)"))
//					System.exit(0);
				settings.getFrgFile().delete();
			} else 
				fragmenter.loadStats();
            Log.info("");
		}
		sequencer= new Sequencer(settings);
		if (settings.getErrFile()!= null&& !sequencer.loadErrors())
			throw new RuntimeException(""); // todo: describe the problem
		if (settings.getSeqFile()!= null&& settings.getSeqFile().exists()&& settings.getSeqFile().canRead()) {
			if (doExpr|| doLib|| doSeq) {
				if (!userCLIconfirm("[ATTENTION] I am going to delete the sequencing file "+settings.getSeqFile().getName()+", please confirm:\n\t(Yes,No,Don't know)"))
					return null;
				settings.getSeqFile().delete();
			} else
				sequencer.loadStats();
            Log.info("");
		}
        Log.info("");

		// do
		long t0= System.currentTimeMillis();


        if (doExpr)
			doExpr();
		else{
			Log.info("you did not ask for expression, I skip it.\n");
        }


        if (doLib)
			doLib();
		else{
			Log.info("you did not want me to construct the library, I skip it.");
        }
        Log.info("");


		if (doSeq)
			doSeq();
		else {
			Log.info("sequencing has not been demanded, skipped.");
        }
	

	    Log.info("\n[END] I finished, took me " + (System.currentTimeMillis() - t0) / 1000 + " sec.");
        return null;
	}
	
	void doExpr() {
//		if (!profiler.isFinishedReadAnnotation()) 
//			profiler.readAnnotation();
//
//		if (!profiler.isReady()) {
//			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
//				System.err.println("[MISSING] I lack parameters for performing expression.");
//			System.exit(-1);
//		}
		profiler.run();
//		if (!profiler.profile()) {
//			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
//				System.err.println("[FATAL] Problem during expression, I exit.");
//			System.exit(-1);
//		}
	}
	
	void doLib() {
		if (!fragmenter.isReady()) {
			throw new RuntimeException("[WHATSUP] I am missing parameters for performing fragmentation.");
		}
		// see run()
/*		if (fragmenter.isFinished()&& Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			if (!userCLIconfirm("[CAREFULLY] There are files describing a constructed libary, do you want to overwrite?\n(Yes,No,Don't care)"))
				System.exit(-1);
*/				
		fragmenter.run();
	}
	
	void doSeq() {
		if (!sequencer.isReady()) {
			throw new RuntimeException("[NONONO] I am missing parameters for sequencing.");
		}
		sequencer.run();
	}
	
	public static final String[] user_yes= new String[] {"yes", "y", "si", "yo", "ja", "ok"},
		user_no= new String[] {"no", "n", "nein", "nope", "nix", "noe"};
	

	public static boolean userCLIconfirm(String message) {
		
		if (Constants.verboseLevel== Constants.VERBOSE_SHUTUP)
			return true;
		
		while (true) {
			StringBuffer sb= new StringBuffer();
			int in;

			if (message!= null) {
				System.err.print(message+" ");
				System.err.flush();
			}
			
			try {
				while((in= System.in.read())!= '\n') 
					sb.append((char) in);
			} catch (Exception e) {
				; // :)
			}
			
			String s= sb.toString().toLowerCase().trim();
			for (int i = 0; i < user_yes.length; i++) 
				if (s.equals(user_yes[i]))
					return true;
			for (int i = 0; i < user_no.length; i++) 
				if (s.equals(user_no[i]))
					return false;
		}
	}
	
	public static double sampleLinear(double rank, double par1, double par2) {
		return par1* rank+ par2;
	}
	
	public static double samplePareto(RandomDataImpl rnd, double par1) {
		double d= 0;
		while(d== 0)
			d= rnd.nextUniform(0, 1);
		double val= (1/Math.pow(d, (1d/par1)));
		return val;
	}
	
	/** Generates pseudo-random numbers under a linear distribution
	 * with a<>0.
	 * @param u
	 * @param a
	 * @param len
     * @param more5
	 * @return
	 */
	public static int sampleInversionLinear(double u, double a, int len, boolean more5) {
//		if (a== 0) {
//			System.err.println("[OOPS] linear function failed a==0");
//			return 0;
//		}
//		
//		boolean invert= false;
//		if (a< 0) {
//			a= -a;
//			invert= true;
//		}
//		
		double b= 0;//len* a;
		//u*= len;
		
//		double r= -(b/a)+ Math.sqrt(Math.pow(b/a, 2)+ (2*u/a));
//		double max= -(b/a)+ Math.sqrt(Math.pow(b/a, 2)+ (2/a));

		double r= Math.sqrt(2*u/ a);
		double max= Math.sqrt(2/ a);
		r/= max;
		r*= len;
		if (more5)
			r= len- r;
		
//		if (true/*new Double(r).equals(Double.NaN)*/) {
//			double x1= Math.pow(b/a, 2);
//			double x2= 2*u/a;
//			double x3= Math.sqrt(x1+ x2);
//			double x4= -(b/a);
//			double x5= x4+ x3;
//			System.currentTimeMillis();
//		}
		
		
		return (int) Math.round(r);
	}
	
	public static double sampleLinearUnderTrpt(double u, int len, double raiseAlong) {
		double a= raiseAlong;
		double b= (a>0)?0:len* a;
		double r= sampleInversionLinear(u, a, b);
		return r;
	}
	
	public static double sampleMPVexpr(double rank, double par1, double par2) {
//		double val= par2/ Math.pow(rank, par1)
//			* Math.exp(- (Math.pow(rank, 2)/ 122000000)
//					- (rank/7000));
		double val= par2/ Math.pow(rank, par1);	// par1= 0,6  par2= (41627d/ 2731598d)
		return val;
	}
	
	public static double sampleMPVdecay(double rank, double par1) {
//		double val= Math.exp(- (Math.pow(rank, 2)/ Math.pow(par1, 2))
//				- (rank/par1));
		double val= Math.exp(- (Math.pow(rank, 2)/ 122000000)
				- (rank/7000));

		return val;
	}


    private class FragmentA implements Comparable {
		int start= 0, end= 0;
		public FragmentA(int a, int b) {
			start= a;
			end= b;
		}
		public int length() {
			return end- start+ 1;
		}
		//@Override
		public int compareTo(Object o) {
			int len2= ((FragmentA) o).length();
			return length()- len2;
		}
		@Override
		public String toString() {
			return "["+start+","+end+"]";
		}
	}
	
	/** Generates pseudo-random numbers under a linear distribution
	 * with a<>0.
	 * @param u
	 * @param a
	 * @param b
	 * @return
	 */
	// f^{-1}(x)= (y-b)/a
	public static double sampleInversionLinear(double u, double a, double b) {
		if (a== 0) {
			System.err.println("[OOPS] linear function failed a==0");
			return u;
		}
		
		boolean invert= false;
		if (a< 0) {
			a= -a;
			invert= true;
		}
		
		double r= -(b/a)+ Math.sqrt(Math.pow(b/a, 2)+ (2*u/a));
		
		if (true/*new Double(r).equals(Double.NaN)*/) {
			double x1= Math.pow(b/a, 2);
			double x2= 2*u/a;
			double x3= Math.sqrt(x1+ x2);
			double x4= -(b/a);
			double x5= x4+ x3;
		}
		
		
		return r;
	}

	/**
	 * Returns a gaussian value sampled from a normal (gaussian) 
	 * distribution with >90% of the function covering the range
	 * between the specified boundaries. Subsequently, the mean 
	 * of the distribution is at <code>min+ ((max-min)/ 2)</code> 
	 * and the standard deviation is <code>(max- mean)/ 2.85</code>
	 * since <0.5 of a gaussian distribution exceed 2.85.
	 */
	private static final double CUT_OFF_GAUSSIAN_VAL = 2.85f;
	public static double nextGaussianBetween(Random random, double min, double max) {
		double rdm= 3d;	// gaussian value, stddev 1
		while (rdm< -CUT_OFF_GAUSSIAN_VAL|| rdm> CUT_OFF_GAUSSIAN_VAL)
			rdm= random.nextGaussian();
		double mid= ((double) min)+ (max- min)/ 2f;
		double realValue= mid+ (rdm* (max-mid)/ CUT_OFF_GAUSSIAN_VAL);
		return realValue; 
	}

	
}
