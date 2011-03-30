package fbi.genome.sequencing.rnaseq.simulation;

import commons.Log;
import commons.TableFormatter;
import commons.file.FileHelper;
import commons.system.SystemInspector;
import fbi.genome.io.SpliceGraphIO;
import fbi.genome.model.IntronModel;
import fbi.genome.model.constants.Constants;
import org.apache.commons.math.random.RandomDataImpl;

import java.io.*;
import java.lang.reflect.Method;
import java.util.*;
import java.util.Map.Entry;

//import gphase.solexa.lp.GraphLPsolver5;
//import gphase.solexa.simulation.Nebulizer;

// TODO check whether selected fragments concord w pro file

public class FluxSimulator extends Thread {

	public static String version= "";
	
	public static final String CLI_ABBREV_DO_EXPR= "expr", CLI_ABBREV_DO_RT= "rt", CLI_ABBREV_DO_FRAG= "frag", CLI_ABBREV_DO_FILTER= "filt", CLI_ABBREV_DO_SEQ= "seq";
	public static final char CLI_SHORT_DO_EXPR= 'x', CLI_SHORT_DO_RT= 'r', CLI_SHORT_DO_FRAG= 'f', CLI_SHORT_DO_FILTER= 'i', CLI_SHORT_DO_SEQ= 's';
	public static final int FLANK5_LEN= 50, FLANK3_LEN= 50;
	public static boolean c= false;

	public static final String PROPERTY_BUILD= "simulator.build";
	public static final String PROPERTY_JDK= "simulator.jdk";
	
	private HashMap<String, Method> cliAbbrevMap, cliLongMap;
	private HashMap<Character, Method> cliShortMap;
	private HashMap<Method, String> cliExplMap;
	
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
		
		FluxSimulator sim= new FluxSimulator();
        if(!sim.parseParameters(args)){
            System.exit(-1);
        }

        if (sim.help) {
			sim.printUsage();
			System.exit(0);
		}
		sim.run();
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


    /**
     * Parse the command line parameters and return true if everything is fine.
     * Return false otherwise
     *
     * @param args the command line args
     * @return valid returns true if all command line arguments are valid
     */
    public boolean parseParameters(String[] args) {
        for (int i = 0; args!= null&& i < args.length; i++) {
            if (args[i].startsWith(Constants.CLI_PAR_LONG))
                i= setParameter(getCLIlongMap().get(args[i].substring(Constants.CLI_PAR_LONG.length())), args, i);
            else if (args[i].startsWith(Constants.CLI_PAR_SHORT)) {
                Method m= getCLIabbrevMap().get(args[i].substring(Constants.CLI_PAR_SHORT.length()));
                if (m== null) {
                    int j= Constants.CLI_PAR_SHORT.length();
                    for (; j < args[i].length(); j++)
                        if (getCLIshortMap().get(args[i].charAt(j))== null)	// Character.toLowerCase()
                            break;
                    if (j< args[i].length()) {
                        Log.error("I did not understand the parameter " + args[i] + "!");
                        return false;
                    } // else
                    int x= i;
                    for (j = Constants.CLI_PAR_SHORT.length(); j < args[i].length(); j++)
                        x= setParameter(getCLIshortMap().get(args[i].charAt(j)), args, x);
                    i= x;
                } else
                    i= setParameter(m, args, i);
            } else {
                Log.error("What do you mean by "+ args[i]+ "?\nRunaway argument or bad monday?");
                return false;
            }
        }
        return true;
    }

    private int setParameter(Method m, String[] args, int i) {
		if (m== null) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[PARLEZVOUS] I do not understand parameter "+ args[i]+"!");
			System.exit(-1);
		}
		
		String[] cc= new String[m.getParameterTypes().length];
		if (cc!= null&& cc.length+ i>= args.length) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[WOLLY] Missing arguments for parameter "+args[i]+"!");
			System.exit(-1);
		}
		for (int j = 0; j < cc.length; j++) 
			cc[j]= args[i+ 1+ j];
		try {
			m.invoke(this, cc);
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[FATAL] Failed to set parameter "+args[i]+"!");
			System.exit(-1);
		}
		return (i+ cc.length);
	}
	
	private HashMap<Character, Method> getCLIshortMap() {
		if (cliShortMap == null) {
			cliShortMap = new HashMap<Character, Method>();
			try {
				cliShortMap.put('p', this.getClass().getDeclaredMethod("setFile", new Class[] {String.class}));
				//cliShortMap.put('X', this.getClass().getDeclaredMethod("setGUI", null));
				cliShortMap.put('h', this.getClass().getDeclaredMethod("setHelp", null));
				cliShortMap.put('x', this.getClass().getDeclaredMethod("setDoExpr", null));
				cliShortMap.put('l', this.getClass().getDeclaredMethod("setDoLib", null));
				cliShortMap.put('s', this.getClass().getDeclaredMethod("setDoSeq", null));
				cliShortMap.put('j', this.getClass().getDeclaredMethod("setDoSJ", new Class[] {String.class}));
				cliShortMap.put('g', this.getClass().getDeclaredMethod("setGenomeDir", new Class[] {String.class}));
			} catch (Exception e) {
				e.printStackTrace();
			}			
		}

		return cliShortMap;
	}
	
	private HashMap<String, Method> getCLIabbrevMap() {
		if (cliAbbrevMap == null) {
			cliAbbrevMap = new HashMap<String, Method>();
			try {
				cliAbbrevMap.put("par", this.getClass().getDeclaredMethod("setFile", new Class[] {String.class}));
				//cliAbbrevMap.put("gui", this.getClass().getDeclaredMethod("setGUI", null));
				cliAbbrevMap.put("expr", this.getClass().getDeclaredMethod("setDoExpr", null));
				cliAbbrevMap.put("lib", this.getClass().getDeclaredMethod("setDoLib", null));
				cliAbbrevMap.put("seq", this.getClass().getDeclaredMethod("setDoSeq", null));
				cliAbbrevMap.put("sj", this.getClass().getDeclaredMethod("setDoSJ", new Class[] {String.class}));
				cliAbbrevMap.put("5flank", this.getClass().getDeclaredMethod("set5flank", new Class[] {String.class}));
				cliAbbrevMap.put("3flank", this.getClass().getDeclaredMethod("set3flank", new Class[] {String.class}));
				cliAbbrevMap.put("imodel", this.getClass().getDeclaredMethod("setImodel", new Class[] {String.class}));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		return cliAbbrevMap;
	}
	private HashMap<String, Method> getCLIlongMap() {
		if (cliLongMap == null) {
			cliLongMap = new HashMap<String, Method>();
			try {
				cliLongMap.put("parameter", this.getClass().getDeclaredMethod("setFile", new Class[] {String.class}));
				//cliLongMap.put("graphical", this.getClass().getDeclaredMethod("setGUI", null));
				cliLongMap.put("help", this.getClass().getDeclaredMethod("setHelp", null));
				cliLongMap.put("express", this.getClass().getDeclaredMethod("setDoExpr", null));
				cliLongMap.put("library", this.getClass().getDeclaredMethod("setDoLib", null));
				cliLongMap.put("sequence", this.getClass().getDeclaredMethod("setDoSeq", null));
				cliLongMap.put("junctions", this.getClass().getDeclaredMethod("setDoSJ", new Class[] {String.class}));
				cliLongMap.put("imodel", this.getClass().getDeclaredMethod("setImodel", new Class[] {String.class}));
				cliLongMap.put("genome", this.getClass().getDeclaredMethod("setGenomeDir", new Class[] {String.class}));
				cliLongMap.put("install", this.getClass().getDeclaredMethod("setDoInstall", null));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		return cliLongMap;
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
	
	
	public static void install() {
		String wrapper= System.getProperty(Constants.PROPERTY_KEY_WRAPPER_BASE);
		try {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[INIT] installing..");
			
			File demoDir= new File(wrapper+ File.separator+ ".."+ File.separator+ "demo");
			if ((!demoDir.exists())|| (!demoDir.isDirectory())) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("\t[FAILED] didn't find \'demo\' folder.");
			}
			String[] files= demoDir.list();
			String demoPath= demoDir.getAbsolutePath();
			try {
				demoPath= demoDir.getCanonicalPath();
			} catch (Exception exx) {
				; // :)
			}
			for (int i = 0; i < files.length; i++) {
				if (!files[i].endsWith(FluxSimulatorSettings.DEF_SFX_PAR))
					continue;
				File tmpFile= File.createTempFile("simulator", "install");
				File parFile= new File(demoPath+ File.separator+ files[i]);
				if (!FileHelper.copy(parFile, tmpFile)) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
						System.err.println("\t[FAILED] couldn't copy to "+ System.getProperty("java.io.tmpdir"));
				}
				BufferedReader buffy= new BufferedReader(new FileReader(tmpFile));
				BufferedWriter writer= new BufferedWriter(new FileWriter(parFile));
				int ctr= 0;
				for (String s = null; (s= buffy.readLine())!= null;++ctr) {
					if (!(s.startsWith(FluxSimulatorSettings.PAR_FRG_FNAME)||
							s.startsWith(FluxSimulatorSettings.PAR_PRO_FNAME)||
							s.startsWith(FluxSimulatorSettings.PAR_SEQ_FNAME)||
							s.startsWith(FluxSimulatorSettings.PAR_REF_FNAME)||
							s.startsWith(FluxSimulatorSettings.PAR_ERR_FNAME))) {
						writer.write(s);
						writer.write(System.getProperty("line.separator"));
						continue;
					}
					String[] ss= s.split("\\s");
					writer.write(ss[0]);
					writer.write("\t");
					int p1= ss[1].lastIndexOf("\\"), p2= ss[1].lastIndexOf('/');
					int p= Math.max(p1, p2);
					// IzPack variable substitution, eg ${INSTALL_PATH}${FILE_SEPARATOR}testRun.pro
					int p3= ss[1].lastIndexOf("}");
					p= Math.max(p, p3);
					writer.write(demoPath+ File.separator+ ss[1].substring(p+1));
					writer.write(System.getProperty("line.separator"));
				}
				buffy.close();
				writer.flush();
				writer.close();
				tmpFile.delete();
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("\twrote "+ctr+" lines to "+parFile.getName());
			}
				

			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("\t[OK]");
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("\t[FAILED] "+ e.getMessage());
		}
	}
	
	
	public void run() {
		
		if (doInstall) {
			install();
			System.exit(0);
		}
		
		if (doSJ&& (doExpr|| doLib|| doSeq)) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				if (doExpr|| doLib|| doSeq)
					System.err.println("[TOOMUCH] Cannot mix splice junction extraction with an RNAseq experiment " +
							"(expression, library construction and sequencing).");
			}
			System.exit(-1);
		} else {
			FluxSimulator.c= SystemInspector.checkRuntime();
		}
		if (doSJ) {
			
			if (file== null|| !file.exists()) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[AIIII] I have no parameter file and I want to scream!");
				System.exit(-1);
			}
			if (fileGenome== null|| !fileGenome.exists()) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[AIAIII] I have no directory with the genomic sequences!");
				System.exit(-1);
			} else
				fbi.genome.model.Graph.overrideSequenceDirPath= fileGenome.getAbsolutePath();

			int[] flanks= getEFlanks();
			flanks[0]= flanks[0]< 0? FLANK5_LEN: flanks[0];
			flanks[1]= flanks[1]< 0? FLANK5_LEN: flanks[1];
			
			IntronModel iModel= new IntronModel();
			if (fileIM!= null)
				iModel.read(fileIM);
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("Extracting splice junctions, 5'sequence "+flanks[0]
					+ ", 3'sequence "+eflanks[1]+", intron model "+ iModel.getName());
				SpliceGraphIO.extractSpliceJunctions(flanks[0], flanks[1], iModel, file, null);
		}


        /*
		if (gui) {
			FluxSimulatorGUI gui= FluxSimulatorGUI.createGUI();
			if (file!= null) {
				settings= gui.load(file);
				if (settings!= null)
					gui.loadInit();
			}
			return;
		}
		*/
		
		// CLI
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println("[HELLO] I am the Flux Simulator (build "+version+"), nice to meet you!\n");
		if (file== null) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[AIIII] I have no parameter file and I want to scream!");
			System.exit(-1);
		}
		
		
		// init
		settings= FluxSimulatorSettings.createSettings(file);
		if (settings== null)
			System.exit(-1);
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println("[INIT] I am collecting information on the run.");
		profiler= new Profiler(settings);
		settings.setProfiler(profiler);
		if (settings.getProFile()!= null&& settings.getProFile().exists()&& settings.getProFile().canRead()) {
			profiler.loadStats();
			if (doExpr&& profiler.isFinishedExpression()) {
				if (!userCLIconfirm("[CAUTION] I overwrite the expression values in file "+settings.getProFile().getName()+", please confirm:\n\t(Yes,No,Don't know)"))
					System.exit(0);	
				else {
					boolean b= settings.getProFile().delete();	// TODO maybe only remove rfreqs..
					if (settings.getProfiler()!= null)
						settings.getProfiler().status= -1;
				}
			}
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println();
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
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println();
		}
		sequencer= new Sequencer(settings);
		if (settings.getErrFile()!= null&& !sequencer.loadErrors())
			System.exit(-1);
		if (settings.getSeqFile()!= null&& settings.getSeqFile().exists()&& settings.getSeqFile().canRead()) {
			if (doExpr|| doLib|| doSeq) {
				if (!userCLIconfirm("[ATTENTION] I am going to delete the sequencing file "+settings.getSeqFile().getName()+", please confirm:\n\t(Yes,No,Don't know)"))
					System.exit(0);
				settings.getSeqFile().delete();
			} else
				sequencer.loadStats();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println();
		}
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println();
		
		// do
		long t0= System.currentTimeMillis();
		if (doExpr) 
			doExpr();
		else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println("[NOEXPR] you did not ask for expression, I skip it.\n");
		if (doLib) 
			doLib();
		else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println("[NOLIB] you did not want me to construct the library, I skip it.");
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println();
		if (doSeq)
			doSeq();
		else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println("[NOSEQ] sequencing has not been demanded, skipped.");
	
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println("\n[END] I finished, took me "+(System.currentTimeMillis()- t0)/ 1000+" sec.");
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
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println(); 	// sep line
	}
	
	void doLib() {
		if (!fragmenter.isReady()) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[WHATSUP] I am missing parameters for performing fragmentation.");
			System.exit(-1);
		}
		// see run()
/*		if (fragmenter.isFinished()&& Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			if (!userCLIconfirm("[CAREFULLY] There are files describing a constructed libary, do you want to overwrite?\n(Yes,No,Don't care)"))
				System.exit(-1);
*/				
		fragmenter.run();
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println(); 	// sep line
	}
	
	void doSeq() {
		if (!sequencer.isReady()) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[NONONO] I am missing parameters for sequencing.");
			System.exit(-1);
		}
		sequencer.run();
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println(); 	// sep line
	}
	
	public static final String[] user_yes= new String[] {"yes", "y", "si", "yo", "ja", "ok"},
		user_no= new String[] {"no", "n", "nein", "nope", "nix", "noe"};
	
	public void printUsage() {
		if (Constants.verboseLevel== Constants.VERBOSE_SHUTUP) 
			return;
		
		System.err.println("[HELLO] I am the Flux Simulator (build "+version+"), nice to meet you!\n");
		System.err.println("[NOARGS] So, here is the list of arguments I understand:\n");
		
		TableFormatter tf= new TableFormatter(3);
		tf.setTabRow(true);
		tf.add(new String[] {"Parameter", "Argument", "Description"});
		Iterator<Method> it= getCliExplMap().keySet().iterator();
		while (it.hasNext()) {
			Method m= it.next();
			StringBuilder sb= new StringBuilder("[");
			Object[] oo= getCLIshortMap().entrySet().toArray();
			for (int i = 0; i < oo.length; i++) {
				Entry<Character, Method> en= (Entry<Character, Method>) oo[i];
				if (en.getValue().equals(m)) {
					sb.append(Constants.CLI_PAR_SHORT);
					sb.append(en.getKey());
					break;
				}
			}
			oo= getCLIabbrevMap().entrySet().toArray();
			for (int i = 0; i < oo.length; i++) {
				Entry<String, Method> en= (Entry<String, Method>) oo[i];
				if (en.getValue().equals(m)) {
					if (sb.length()> 1)
						sb.append("|");
					sb.append(Constants.CLI_PAR_SHORT);
					sb.append(en.getKey());
					break;
				}
			}
			oo= getCLIlongMap().entrySet().toArray();
			for (int i = 0; i < oo.length; i++) {
				Entry<String, Method> en= (Entry<String, Method>) oo[i];
				if (en.getValue().equals(m)) {
					if (sb.length()> 1)
						sb.append("|");
					sb.append(Constants.CLI_PAR_LONG);
					sb.append(en.getKey());
					break;
				}
			}
			sb.append("]");
			String s1= sb.toString();
			
			sb= new StringBuilder();
			if (m.getParameterTypes()!= null&& m.getParameterTypes().length> 0) {
				for (int i = 0; i < m.getParameterTypes().length; i++) {
					sb.append(m.getParameterTypes()[i].getSimpleName().toString());
					sb.append(",");
				}
				sb.deleteCharAt(sb.length()- 1);
			}
			String s2= sb.toString();
			String s3= getCliExplMap().get(m);
			
			tf.add(new String[] {s1,s2,s3});
		}
		
		System.err.println(tf.toString());
		System.err.println("\nByebye.");
	}
	
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

	private HashMap<Method, String> getCliExplMap() {
		if (cliExplMap == null) {
			cliExplMap = new HashMap<Method, String>();
			try {
				cliExplMap.put(this.getClass().getDeclaredMethod("setFile", new Class[] {String.class}), "specify parameter file (PAR file)");
				//cliExplMap.put(this.getClass().getDeclaredMethod("setGUI", null), "start graphical user interface (GUI)");
				cliExplMap.put(this.getClass().getDeclaredMethod("setHelp", null), "request command line options");
				cliExplMap.put(this.getClass().getDeclaredMethod("setDoExpr", null), "simulate expression");
				cliExplMap.put(this.getClass().getDeclaredMethod("setDoLib", null), "simulate library construction");
				cliExplMap.put(this.getClass().getDeclaredMethod("setDoSeq", null), "simulate sequencing");
				cliExplMap.put(this.getClass().getDeclaredMethod("setDoSJ", new Class[] {String.class}), "extract splice junctions (GTF file)");
				cliExplMap.put(this.getClass().getDeclaredMethod("set5flank", new Class[] {String.class}), "exonic flank 5' of intron (-sj)");
				cliExplMap.put(this.getClass().getDeclaredMethod("set3flank", new Class[] {String.class}), "exonic flank 3' of intron (-sj)");
				cliExplMap.put(this.getClass().getDeclaredMethod("setImodel", new Class[] {String.class}), "specify the introm model (-sj)");
				cliExplMap.put(this.getClass().getDeclaredMethod("setGenomeDir", new Class[] {String.class}), "set the path to the directory with genomic sequences (-sj)");
				cliExplMap.put(this.getClass().getDeclaredMethod("setDoInstall", null), "install the demonstration (.par) files");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		return cliExplMap;
	}

	
}
