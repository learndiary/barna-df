package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.file.FileHelper;
import fbi.genome.model.constants.Constants;

import java.io.*;
import java.util.Hashtable;
import java.util.Iterator;

public class FluxSimulatorSettings {

	public static final String TMP_SFX= ".tmp";
	
	public static final String NA= "NA";
	
	public static byte DISTR_NI= -1, DISTR_UNIFORM= 0, DISTR_NORMAL= 1, DISTR_GAUSS= 2, DISTR_POISSON= 3;
	public static int[] DISTR_NR_PAR= new int[] {0,0,2,1};
	public static final String PAR_DISTR_UNIFORM= "UNIFORM", PAR_DISTR_NORMAL= "NORMAL", PAR_DISTR_GAUSS= "GAUSS", PAR_DISTR_POISSON= "POISSON", PAR_DISTR_WEIBULL= "WEIBULL";
	public static String[] DISTR_NAMES= new String[] {PAR_DISTR_UNIFORM, PAR_DISTR_NORMAL, PAR_DISTR_GAUSS, PAR_DISTR_POISSON};
	
	public static final char DELIM_FURI= '#';	// flux unique read identifier
	public static final char DELIM_FMOLI= ':';	// flux molecule identifier
	public static final int PRO_COL_NR_MOL= 4, PRO_COL_NR_FRG= 6, PRO_COL_NR_SEQ= 8;
	public static final String PRO_FILE_SEP= "\t";
	public static final ByteArrayCharSequence PRO_FILE_CDS= new ByteArrayCharSequence("CDS"), PRO_FILE_NC= new ByteArrayCharSequence("NC");
	public static final String PAR_TMP_FNAME= "TMP_DIR", PAR_PRO_FNAME= "PRO_FILE_NAME", PAR_FRG_FNAME= "LIB_FILE_NAME", PAR_SEQ_FNAME= "SEQ_FILE_NAME", PAR_REF_FNAME= "REF_FILE_NAME", PAR_ERR_FNAME= "ERR_FILE_NAME";
	public static final String PAR_RUN_NAME = "RUN_NAME";
	public static final String PAR_LOAD_CODING = "LOAD_CODING", PAR_LOAD_NONCODING= "LOAD_NONCODING";
	public static final String PAR_YES = "YES", PAR_NO = "NO";
	public static final String PAR_SEQ_READ_LENGTH= "READ_LENGTH", PAR_SEQ_READ_NUMBER= "READ_NUMBER", PAR_SEQ_PEND= "PAIRED_END", PAR_SEQ_FASTQ= "FASTQ", PAR_SEQ_QTHOLD= "QTHOLD";
	public static final String PAR_FRAG_LAMBDA = "FRAG_LAMBDA";
	public static final String PAR_FRAG_SIGMA = "FRAG_SIGMA";
	public static final String PAR_TSS_MEAN= "TSS_MEAN", PAR_POLYA_SHAPE= "POLYA_SHAPE", PAR_POLYA_SCALE= "POLYA_SCALE";
	public static final String PAR_FRAG_THRESHOLD = "FRAG_THRESHOLD";
	public static final String PAR_RT_FALLOFF = "RT_FALLOFF";
	public static final String PAR_RT_RANDOM = "RT_RANDOM";
	public static final String PAR_FRAG_MODE= "FRAG_MODE", PAR_FRAG_MODE_PHYS= "PHYSICAL", PAR_FRAG_MODE_CHEM= "CHEMICAL";
	public static final String PAR_BP_DISTR= "BP_DISTR", PAR_BP_MOTIF= "BP_MOTIF";
	public static final String PAR_RT_MODE= "RT_PRIMER", PAR_RT_MODE_POLY_DT = "POLY-DT", PAR_RT_MODE_RANDOM = "RANDOM", PAR_RT_MODE_NONE = "NONE"; 
	public static final String PAR_RT_MIN= "RT_MIN", PAR_RT_MAX= "RT_MAX";
	public static final String PAR_FRAG= "FRAGMENTATION", PAR_FILT= "FILTERING";
	public static final String PAR_DEC_P1 = "EXPRESSION_X1";
	public static final String PAR_EXPR_K = "EXPRESSION_K";
	public static final String PAR_EXPR_X0 = "EXPRESSION_X0"; 
	public static final String PAR_NB_CELLS = "NB_CELLS";
	public static final String PAR_MAXTHREAD = "MAXTHREAD";
	public static final String PAR_NB_MOLECULES = "NB_MOLECULES";
	public static final String PAR_FILT_MIN= "FILT_MIN", PAR_FILT_MAX= "FILT_MAX", PAR_FILT_DISTR= "PAR_FILT_DISTR";
	public static final String PAR_FRAG_B4_RT= "FRAG_B4_RT";
	public static final String PAR_GEN_DIR= "GEN_DIR";
	public static String PRO_FILE_CR= "\n", PRO_FILE_TAB= "\t", PAR_COMMENT= "#", TMP_PFX= "sim";
	public static String PAR_PWM_FRAG_SENSE= "PWM_FRAG_SENSE", PAR_PWM_FRAG_ASENSE= "PWM_FRAG_ASENSE", PAR_PWM_RT_ASENSE= "PWM_RT_ASENSE", PAR_PWM_RT_SENSE= "PWM_RT_SENSE";

	public static boolean optDisk= false;
	public static final char SEP_LOC_TID= '@';
	
	public static final int AVG_MOL_CELL= 100000;
	
	public String checkComplete() {
		
		StringBuilder sb= new StringBuilder();
		if (getRefFile()== null|| !getRefFile().exists())
			sb.append("\nmissing: "+ PAR_REF_FNAME);
		if (getProFile()== null)
			sb.append("\nmissing: "+ PAR_PRO_FNAME);
		if (getFrgFile()== null)
			sb.append("\nmissing: "+ PAR_FRG_FNAME);
		if (getSeqFile()== null)
			sb.append("\nmissing: "+ PAR_SEQ_FNAME);

		return sb.toString();
	}
	
	public static boolean checkNewDialogOK(boolean output, String fNameAnn, String fNameTmp, String[] files, boolean checkOverwrite) throws RuntimeException {

		for (int i = 0; i < files.length; i++) {
			File f= new File(files[i]);
			if (files[i]== null) {
				if (output)
					Constants.dialog.showError("Aiaiai, there are file names missing, and I need them!");
				return false;
			}
			if (f.isDirectory()) {
				if (output)
					Constants.dialog.showError("Nonono,\n"+ files[i]+"\nis a directory and I expect a file!");
				return false;
			}
			if (f.exists()) {
				if (checkOverwrite)
					if (!Constants.dialog.checkOverwrite("Ohlala, do you really want to overwrite this file\n"+ files[i]))
						return false;
					else if (!f.delete()) {						
						//throw new RuntimeException("Could not delete\n "+ files[i]);
						if (output)
							Constants.dialog.showError("Ups, I could not delete\n"+ files[i]+"\nremove open file handles or try another file.");
						return false;
					}
			}
			// TODO: linux bug in 
			// (f.exists()&& (!f.canWrite()))|| (!f.exists())&& (f.getParentFile()== null|| (!f.getParentFile().canWrite()))
			if (!FileHelper.canWrite(f)) {
				if (output)
					Constants.dialog.showError("Uuuh, cannot write to file\n"+ files[i]+", try again!");
				return false;
			}
		}
		if (fNameAnn== null) {
			if (output)
				Constants.dialog.showError("Aiiii, I have no annotation file, and I want to scream!");
			return false;
		}
		File fAnn= new File(fNameAnn), fTmp= new File(fNameTmp);
		if (!(fAnn.exists()&& fAnn.canRead())) {
			if (output)
				Constants.dialog.showError("Oooh, problems with annotation file\n"+ fAnn+", nice try!");
			return false;
		}
		// && FileHelper.canWriteToDir(fTmp))
		if (fTmp!= null&& !(fTmp.exists()&& fTmp.isDirectory()&& FileHelper.canWrite(fTmp))) {
			if (output)
				Constants.dialog.showError("Ohoh, problems with the temporary directory\n"+ fTmp+", another time please!");
			return false;
		}
		return true;
	}
	
	public static boolean writeParfile(File f, FluxSimulatorSettings settings) {
		try {
			BufferedWriter writer= new BufferedWriter(new FileWriter(f));
			
			// project
			String s= settings.getRefFile().getAbsolutePath(); 
/*			FileHelper.getRelativePath(f, settings.getRefFile());
			if (!new File(f.getParentFile().getPath()+File.separator+s).exists())
				s= settings.getRefFile().getPath();
			else
				s= s.replace(File.separatorChar, '/');
*/				
			writer.write(PAR_REF_FNAME+ "\t"+ s+ "\n");
			s= settings.getProFile().getAbsolutePath();
			/*FileHelper.getRelativePath(f, settings.getProFile());
			if (!new File(f.getParentFile().getPath()+File.separator+s).getParentFile().exists())
				s= settings.getProFile().getPath();
			else
				s= s.replace(File.separatorChar, '/');*/				
			writer.write(PAR_PRO_FNAME+ "\t"+ s+ "\n");
			s= settings.getFrgFile().getAbsolutePath();
			/*FileHelper.getRelativePath(f, settings.getFrgFile());
			if (!new File(f.getParentFile().getPath()+File.separator+s).getParentFile().exists())
				s= settings.getFrgFile().getPath();
			else
				s= s.replace(File.separatorChar, '/');*/
			writer.write(PAR_FRG_FNAME+ "\t"+ s+ "\n");
			s= settings.getSeqFile().getAbsolutePath();
			/*FileHelper.getRelativePath(f, settings.getSeqFile());
			if (!new File(f.getParentFile().getPath()+File.separator+s).getParentFile().exists())
				s= settings.getSeqFile().getPath();
			else
				s= s.replace(File.separatorChar, '/');*/
			writer.write(PAR_SEQ_FNAME+ "\t"+ s+ "\n");
			if (settings.getGenDir()!= null) {
				s= settings.getGenDir().getAbsolutePath();
				writer.write(PAR_GEN_DIR+ "\t"+ s+ "\n");
			}
			if (settings.getTmpDir()!= null) 
					//&& !settings.getTmpDir().getPath().equalsIgnoreCase(new File(System.getProperty("java.io.tmpdir")).getPath())) 
				{
				s= settings.getTmpDir().getAbsolutePath();
				/*FileHelper.getRelativePath(f, settings.getTmpDir());
				if (!new File(f.getParentFile().getPath()+File.separator+s).exists())
					s= settings.getTmpDir().getPath();
				else
					s= s.replace(File.separatorChar, '/');*/
				writer.write(PAR_TMP_FNAME+ "\t"+ s+ "\n");
			}

			// profiler
			writer.write(PAR_NB_MOLECULES+ "\t"+ settings.getNbMolecules()+ "\n");
			writer.write(PAR_EXPR_K+ "\t"+ settings.getExpDistrP1()+ "\n");
			writer.write(PAR_EXPR_X0+ "\t"+ settings.getExpDistrP2()+ "\n");
			writer.write(PAR_DEC_P1+ "\t"+ settings.getDecDistrP1()+ "\n");
			if (settings.isTssVar()) {
				writer.write(PAR_TSS_MEAN+ "\t");
				writer.write(Double.toString(settings.getTssMean()));
				writer.write("\n");
			}
			if (settings.isPolyAVar()) {
				writer.write(PAR_POLYA_SHAPE+ "\t");
				writer.write(Double.toString(settings.getPolyAshape()));
				writer.write("\n");
				writer.write(PAR_POLYA_SCALE+ "\t");
				writer.write(Double.toString(settings.getPolyAscale()));
				writer.write("\n");
			}
			
			// fragmentation
			writer.write(PAR_RT_MIN+ "\t"+ settings.getMinRTLen()+ "\n");
			writer.write(PAR_RT_MAX+ "\t"+ settings.getMaxRTLen()+ "\n");
			writer.write(PAR_RT_MODE+ "\t"+ settings.getRtMode()+ "\n");
			//writer.write(PAR_RT_RANDOM+ "\t"+ settings.getRtRandom()+ "\n");
			//writer.write(PAR_RT_FALLOFF+ "\t"+ settings.getRtFoff()+ "\n");
			writer.write(PAR_FRAG+ "\t");
			if (settings.isFragment())
				writer.write(PAR_YES+ "\n");
			else
				writer.write(PAR_NO+ "\n");
			writer.write(PAR_FRAG_B4_RT+ "\t");
			if (settings.isFragB4RT())
				writer.write(PAR_YES+ "\n");
			else
				writer.write(PAR_NO+ "\n");
			writer.write(PAR_FRAG_MODE+ "\t"+ settings.getFragMode()+ "\n");
			writer.write(PAR_FRAG_LAMBDA+ "\t"+ settings.getLambda()+ "\n");
			//writer.write(PAR_FRAG_SIGMA+ "\t"+ settings.getSigma()+ "\n");
			writer.write(PAR_FRAG_THRESHOLD+ "\t"+ settings.getThold()+ "\n");
			writer.write(PAR_FILT+ "\t");
			if (settings.isFilter())
				writer.write(PAR_YES+ "\n");
			else
				writer.write(PAR_NO+ "\n");
			writer.write(PAR_LOAD_CODING+ "\t");
			if (settings.isLoadCoding())
				writer.write(PAR_YES+ "\n");
			else
				writer.write(PAR_NO+ "\n");
			writer.write(PAR_LOAD_NONCODING+ "\t");
			if (settings.isLoadNoncoding())
				writer.write(PAR_YES+ "\n");
			else
				writer.write(PAR_NO+ "\n");
			writer.write(PAR_FILT_MIN+ "\t"+ settings.getFiltMin()+ "\n");
			writer.write(PAR_FILT_MAX+ "\t"+ settings.getFiltMax()+ "\n");

			// sequencing
			writer.write(PAR_SEQ_READ_NUMBER+ "\t"+ settings.getReadNr()+ "\n");
			writer.write(PAR_SEQ_READ_LENGTH+ "\t"+ settings.getReadLength()+ "\n");
			writer.write(PAR_SEQ_PEND+ "\t"+ (settings.isPairedEnd()?PAR_YES:PAR_NO)+ "\n");
			writer.write(PAR_SEQ_FASTQ+ "\t"+ (settings.isFastQ()?PAR_YES:PAR_NO)+ "\n");
			if (settings.isFastQ()) {
				if (!Float.isNaN(settings.getQthold()))
					writer.write(PAR_SEQ_QTHOLD+ "\t"+ Float.toString(settings.qthold)+ "\n");
				if (settings.getErrFile()!= null) {
					/*s= FileHelper.getRelativePath(f, settings.getErrFile());
					if (!new File(f.getParentFile().getPath()+File.separator+s).getParentFile().exists())
						s= settings.getErrFile().getPath();
					else
						s= s.replace(File.separatorChar, '/');*/
					s= settings.getErrFile().getAbsolutePath();
					try {
						s= settings.getErrFile().getCanonicalPath();
					} catch (Exception e) {
						; // :)
					}
					writer.write(PAR_ERR_FNAME+ "\t"+ s+ "\n");
				}
			}
			writer.flush();
			writer.close();
			return true;
			
		} catch (Exception e) {
			return false;
		}
	}
	public static boolean appendProfile(FluxSimulatorSettings settings, int colNr, Hashtable<CharSequence,Long> mapFrags) {
		try {
			if (Constants.progress!= null)
				Constants.progress.setString("Updating .pro file ");
			
			long total= 0;
			Iterator<Long> iter= mapFrags.values().iterator();
			while(iter.hasNext())
				total+= iter.next();
			
			BufferedReader buffy= new BufferedReader(new FileReader(settings.proFile));
			File tmpF= new File(System.getProperty(Constants.PROPERTY_TMPDIR)+ File.separator+settings.proFile.getName()+TMP_SFX);
			BufferedWriter wright= new BufferedWriter(new FileWriter(tmpF));
			String[] token;
			long bytesRead= 0, bytesTotal= settings.proFile.length();
			int perc= 0, lineCtr= 0;
			String nullStr= Double.toString(0d)+PRO_FILE_TAB+Long.toString(0);
			for (String s= null; (s= buffy.readLine())!= null;++lineCtr) {
				
				bytesRead+= s.length()+ PRO_FILE_CR.length();
				if (lineCtr% 1000== 0&& bytesRead* 10d/ bytesTotal> perc) {
					++perc;
					if (Constants.progress!= null)
						Constants.progress.progress();	// setValue(perc)
				}
				
				token= s.split(FluxSimulatorSettings.PRO_FILE_SEP);
				if (token.length== colNr) {
					wright.write(s);
					wright.write(FluxSimulatorSettings.PRO_FILE_SEP);
				} else
					for (int i = 0; i < colNr; i++) {
						wright.write(token[i]);
						wright.write(FluxSimulatorSettings.PRO_FILE_SEP);
					}
				String id= token[0]+ "@"+ token[1];
				if (mapFrags.containsKey(id)) {
					long absCnt= mapFrags.get(id);
					if (FluxSimulator.c&& false) {
						Iterator<CharSequence> kk= mapFrags.keySet().iterator();
						for(CharSequence kx= null,k;kk.hasNext();kx=k) {
							k= kk.next();
							if (token[1].equals(k)) {
								if (kk.hasNext()) {
									absCnt= mapFrags.get(kk.next());
									break;
								}
								if (kx!= null) {
									absCnt= mapFrags.get(kx);
									break;
								}
								absCnt= mapFrags.get(k);
								break;
							}
						}
					}
					double relFreq= absCnt/ (double) total;
					if (Double.isNaN(relFreq))
						System.currentTimeMillis();
					wright.write(Double.toString(relFreq));
					wright.write(FluxSimulatorSettings.PRO_FILE_SEP);
					wright.write(Long.toString(absCnt));
				} else 
					wright.write(nullStr);
				
				wright.write(PRO_FILE_CR);
				if (lineCtr%1000== 0)
					wright.flush();
			}
			buffy.close();
			wright.flush();
			wright.close();
			
			settings.proFile.delete();
			FileHelper.move(tmpF, settings.proFile, null);
			
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println(" OK");
			return true;
			
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
	}
	
	public static FluxSimulatorSettings createSettings(File f) {
	try {
			
		if (f== null|| !f.exists()) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				if (f== null)
					System.err.println("\t[UHOH] I have no parameter file and I want to scream!");
				else
					System.err.println("\t[UHOH] parameter file does not exist: "+f.getCanonicalPath());
			}
			return null;
		} else {
			f= new File(f.getCanonicalPath());	// kill Win32ShellFolder instances, they fuck up relative path conversion
		}
		
		FluxSimulatorSettings settings= new FluxSimulatorSettings();
		fillDefaults(settings);
		settings.setParFile(f);
		
		BufferedReader buffy= new BufferedReader(new FileReader(f));
		String refname= null, proname= null, frgname= null, errname= null, seqname= null, tmpname= null, genname= null;
		for (String s; (s= buffy.readLine())!= null; ) {
			s= s.trim();
			if (s.startsWith(PAR_COMMENT))
				continue;
			else if (s.startsWith(PAR_RUN_NAME))
				settings.name= s.substring(PAR_RUN_NAME.length()).trim();
			else if (s.startsWith(PAR_SEQ_READ_LENGTH)) 
				settings.readLength= Integer.parseInt(s.substring(PAR_SEQ_READ_LENGTH.length()).trim());
			else if (s.startsWith(PAR_PRO_FNAME)) 
				proname= s.substring(PAR_PRO_FNAME.length()).trim();
			else if (s.startsWith(PAR_REF_FNAME))
				refname= s.substring(PAR_REF_FNAME.length()).trim();
			else if (s.startsWith(PAR_FRG_FNAME))
				frgname= s.substring(PAR_FRG_FNAME.length()).trim();
			else if (s.startsWith(PAR_ERR_FNAME))
				errname= s.substring(PAR_ERR_FNAME.length()).trim();
			else if (s.startsWith(PAR_SEQ_FNAME))
				seqname= s.substring(PAR_SEQ_FNAME.length()).trim();
			else if (s.startsWith(PAR_TMP_FNAME)) {
				tmpname= s.substring(PAR_TMP_FNAME.length()).trim();
				System.setProperty(Constants.PROPERTY_TMPDIR, tmpname);
			} else if (s.startsWith(PAR_GEN_DIR))
				genname= s.substring(PAR_GEN_DIR.length()).trim();
			if (s.startsWith(PAR_COMMENT))
				continue;
			
			// EXPRESSION
			else if (s.startsWith(PAR_EXPR_K))
				try {
					s= s.substring(PAR_EXPR_K.length()).trim();
					settings.expDistrP1= Double.parseDouble(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_EXPR_K+": "+s);
					return null;
				}									
			else if (s.startsWith(PAR_EXPR_X0)) 
				try {
					s= s.substring(PAR_EXPR_X0.length()).trim();
					settings.expDistrP2= Double.parseDouble(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_EXPR_X0+": "+s);
					return null;
				}				
			else if (s.startsWith(PAR_DEC_P1))
				try {
					s= s.substring(PAR_DEC_P1.length()).trim();
					settings.decDistrP1= Double.parseDouble(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_DEC_P1+": "+s);
					return null;
				}				
//				else if (s.startsWith(PAR_NB_CELLS))
//					try {
//						s= s.substring(PAR_NB_CELLS.length()).trim();
//						settings.nbCells= Long.parseLong(s);
//					} catch (NumberFormatException e) {
//						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
//							System.err.println("\tNo valid value for "+PAR_NB_CELLS+": "+s);
//						return null;
//					}				
			else if (s.startsWith(PAR_NB_MOLECULES))
				try {
					s= s.substring(PAR_NB_MOLECULES.length()).trim();
					settings.nbMolecules= Long.parseLong(s);
					settings.nbCells= settings.nbMolecules/ AVG_MOL_CELL;
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_NB_MOLECULES+": "+s);
					return null;
				}							
			else if (s.startsWith(PAR_TSS_MEAN))
				try {
					s= s.substring(PAR_TSS_MEAN.length()).trim();
					if (!s.toUpperCase().contains(NA))
						settings.tssMean= Double.parseDouble(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_TSS_MEAN+": "+s);
					return null;
				}							
			else if (s.startsWith(PAR_POLYA_SHAPE))
				try {
					s= s.substring(PAR_POLYA_SHAPE.length()).trim();
					if (!s.toUpperCase().contains(NA))
						settings.polyAshape= Double.parseDouble(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_POLYA_SHAPE+": "+s);
					return null;
				}							
			else if (s.startsWith(PAR_POLYA_SCALE))
				try {
					s= s.substring(PAR_POLYA_SCALE.length()).trim();
					//if (!s.toUpperCase().contains(NA))	// 110108: allow NaN to override defaults
						settings.polyAscale= Double.parseDouble(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_POLYA_SCALE+": "+s);
					return null;
				}							
			
		
				// RT / FRAGMENTER				
			else if (s.startsWith(PAR_FRAG_LAMBDA))
				try {
					s= s.substring(PAR_FRAG_LAMBDA.length()).trim();
					settings.lambda= Double.parseDouble(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_FRAG_LAMBDA+": "+s);
					return null;
				}
			else if (s.startsWith(PAR_FRAG_SIGMA))
				try {
					s= s.substring(PAR_FRAG_SIGMA.length()).trim();
					settings.sigma= Double.parseDouble(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_FRAG_SIGMA+": "+s);
					return null;
				}				
			else if (s.startsWith(PAR_FRAG_THRESHOLD))
				try {
					s= s.substring(PAR_FRAG_THRESHOLD.length()).trim();
					settings.thold= Double.parseDouble(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_FRAG_THRESHOLD+": "+s);
					return null;
				}
			else if (s.startsWith(PAR_RT_MIN))
				try {
					s= s.substring(PAR_RT_MIN.length()).trim();
					settings.minRTLen= Integer.parseInt(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_RT_MIN+": "+s);
					return null;
				}				
			else if (s.startsWith(PAR_RT_MAX))
				try {
					s= s.substring(PAR_RT_MAX.length()).trim();
					settings.maxRTLen= Integer.parseInt(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_RT_MAX+": "+s);
					return null;
				}				
			else if (s.startsWith(PAR_FILT_DISTR)) {
				s= s.substring(PAR_FILT_DISTR.length()).trim();
				settings.fileFilterDistr= new File(s);
			} else if (s.startsWith(PAR_FILT_MIN))
				try {
					s= s.substring(PAR_FILT_MIN.length()).trim();
					settings.filtMin= Integer.parseInt(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_FILT_MIN+": "+s);
					return null;
				}				
			else if (s.startsWith(PAR_FILT_MAX))
				try {
					s= s.substring(PAR_FILT_MAX.length()).trim();
					settings.filtMax= Integer.parseInt(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_FILT_MAX+": "+s);
					return null;
				}
			else if (s.startsWith(PAR_PWM_FRAG_SENSE)) {
				s= s.substring(PAR_PWM_FRAG_SENSE.length()).trim();
				settings.filePWMfragSense= new File(s);
			}
			else if (s.startsWith(PAR_PWM_FRAG_ASENSE)) {
				s= s.substring(PAR_PWM_FRAG_ASENSE.length()).trim();
				settings.filePWMfragAsense= new File(s);
			}
			else if (s.startsWith(PAR_PWM_RT_SENSE)) {
				s= s.substring(PAR_PWM_RT_SENSE.length()).trim();
				settings.filePWMrtSense= new File(s);
			}
			else if (s.startsWith(PAR_PWM_RT_ASENSE)) {
				s= s.substring(PAR_PWM_RT_ASENSE.length()).trim();
				settings.filePWMrtAsense= new File(s);
			}
				
				
				
			else if (s.startsWith(PAR_LOAD_CODING)) {
				s= s.substring(PAR_LOAD_CODING.length()).trim();
				settings.setLoadCoding(s.equals(PAR_YES));
			} else if (s.startsWith(PAR_LOAD_NONCODING)) {
				s= s.substring(PAR_LOAD_NONCODING.length()).trim();
				settings.setLoadNoncoding(s.equals(PAR_YES));
			
			} else if (s.startsWith(PAR_FRAG)) {
				s= s.substring(PAR_FRAG.length()).trim();
				settings.setFragment(s.equals(PAR_YES));
			} else if (s.startsWith(PAR_BP_DISTR)) {
				s= s.substring(PAR_FILT.length()).trim();
				
				if (s.startsWith(PAR_DISTR_UNIFORM)) {
					settings.bpDistr= DISTR_UNIFORM;
					s= s.substring(PAR_DISTR_UNIFORM.length()).trim();
				} else if (s.equals(PAR_DISTR_NORMAL)) {
					settings.bpDistr= DISTR_NORMAL;
					s= s.substring(PAR_DISTR_NORMAL.length()).trim();
				} else if (s.equals(PAR_DISTR_GAUSS)) {
					settings.bpDistr= DISTR_GAUSS;
					s= s.substring(PAR_DISTR_GAUSS.length()).trim();
				} else {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_BP_DISTR+": "+s);
					return null;
				}
				settings.bpDistrPar= new float[DISTR_NR_PAR[settings.bpDistr]];
				String[] ss= s.split("\\s");
				if (ss.length!= settings.bpDistrPar.length) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tInvalid number of arguments for "+DISTR_NAMES[settings.bpDistr]+": "+ss.length+", expected "+ settings.bpDistrPar.length);
					return null;
				}
				for (int i = 0; i < ss.length; i++) {
					try {
						settings.bpDistrPar[i]= Float.parseFloat(ss[i]);
					} catch (NumberFormatException e) {
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
							System.err.println("\tNot a float number "+ss[i]);
						return null;
					}
				}
			} else if (s.startsWith(PAR_BP_MOTIF)) {
				s= s.substring(PAR_BP_MOTIF.length()).trim();
				settings.bpMotif= new File(s);
			} else if (s.startsWith(PAR_FILT)) {
				s= s.substring(PAR_FILT.length()).trim();
				settings.setFilter(s.equals(PAR_YES));
			} else if (s.startsWith(PAR_FILT)) {
				s= s.substring(PAR_FILT.length()).trim();
				settings.setFilter(s.equals(PAR_YES));
			} else if (s.startsWith(PAR_FRAG_B4_RT)) {
				s= s.substring(PAR_FRAG_B4_RT.length()).trim();
				settings.setFragB4RT(s.equals(PAR_YES));
			} else if (s.startsWith(PAR_TMP_FNAME))
				settings.tmpDir= new File(s.substring(PAR_TMP_FNAME.length()).trim());
			else if (s.startsWith(PAR_FRAG_MODE)) {
				s= s.substring(PAR_FRAG_MODE.length()).trim();
				settings.setFragMode(s);
			} else if (s.startsWith(PAR_RT_MODE)) {
				s= s.substring(PAR_RT_MODE.length()).trim();
				settings.setRtMode(s);
			} else if (s.startsWith(PAR_SEQ_READ_LENGTH))
				try {
					s= s.substring(PAR_SEQ_READ_LENGTH.length()).trim();
					settings.readLength= Integer.parseInt(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_SEQ_READ_LENGTH+": "+s);
					return null;
				}
			else if (s.startsWith(PAR_SEQ_READ_NUMBER))
				try {
					s= s.substring(PAR_SEQ_READ_NUMBER.length()).trim();
					settings.readNr= Long.parseLong(s);
				} catch (NumberFormatException e) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("\tNo valid value for "+PAR_SEQ_READ_NUMBER+": "+s);
					return null;
				}
			else if (s.startsWith(PAR_SEQ_PEND)) {
				s= s.substring(PAR_SEQ_PEND.length()).trim();
				settings.pairedEnd= s.equals(PAR_YES); 
			} else if (s.startsWith(PAR_SEQ_FASTQ)) {
				s= s.substring(PAR_SEQ_FASTQ.length()).trim();
				settings.fastQ= s.equals(PAR_YES); 
			} else if (s.startsWith(PAR_SEQ_QTHOLD)) {
				s= s.substring(PAR_SEQ_QTHOLD.length()).trim();
				settings.qthold= Float.parseFloat(s); 
			} else if (s.startsWith(PAR_MAXTHREAD)) {
				s= s.substring(PAR_MAXTHREAD.length()).trim();
				try {
					int x= Integer.parseInt(s);
					if (x< 1|| x> 20)
						System.err.println("\tNot a valid number of threads for "+ PAR_MAXTHREAD+ ": "+ s);
					else
						settings.maxThreads= x;
				} catch (NumberFormatException e) {
					System.err.println("\tNot a valid number of threads for "+ PAR_MAXTHREAD+ ": "+ s);
				}
			}

			
		}
		buffy.close();

			// init rel file names
/*			if (refname!= null) {
				try {
					settings.refFile= new File(f.getParentFile().getPath()+ File.separator+ refname.replace('/', File.separatorChar)).getAbsoluteFile();
					if (!settings.refFile.exists()) 
						settings.refFile= new File(refname);
				} catch (Exception e) {
					//if (!settings.refFile.exists())
					settings.refFile= new File(refname);
				}
			}
			if (proname!= null) {
				try {
					settings.proFile= new File(f.getParentFile().getPath()+ File.separator+ proname.replace('/', File.separatorChar)).getAbsoluteFile();
					if (!settings.proFile.getParentFile().exists())
						settings.proFile= new File(proname);
				} catch (Exception e) {
					//if (!settings.proFile.getParentFile().exists())
					settings.proFile= new File(proname);
				}
			}
			if (frgname!= null) {
				try {
					settings.frgFile= new File(f.getParentFile().getPath()+ File.separator+ frgname.replace('/', File.separatorChar)).getAbsoluteFile();
					if (!settings.frgFile.getParentFile().exists())
						settings.frgFile= new File(frgname);
				} catch (Exception e) { 
					//if (!settings.frgFile.getParentFile().exists())
					settings.frgFile= new File(frgname);
				}
			}
			if (errname!= null) {
				try {
					settings.errFile= new File(f.getParentFile().getPath()+ File.separator+ errname.replace('/', File.separatorChar)).getAbsoluteFile();
				} catch (Exception e) { 
					//if (!settings.seqFile.getParentFile().exists())
					settings.errFile= new File(errname);
				}
			}
			if (seqname!= null) {
				try {
					settings.seqFile= new File(f.getParentFile().getPath()+ File.separator+ seqname.replace('/', File.separatorChar)).getAbsoluteFile();
					if (!settings.seqFile.getParentFile().exists())
						settings.seqFile= new File(seqname);
				} catch (Exception e) { 
					//if (!settings.seqFile.getParentFile().exists())
					settings.seqFile= new File(seqname);
				}
			}
			if (tmpname!= null) {
				try {
					settings.tmpDir= new File(f.getParentFile().getPath()+ File.separator+ tmpname.replace('/', File.separatorChar)).getAbsoluteFile();
					if (!settings.tmpDir.getParentFile().exists())
						settings.tmpDir= new File(tmpname);
				} catch (Exception e) { 
					//if (!settings.tmpDir.exists())
					settings.tmpDir= new File(tmpname);
				}
			} else {
				settings.tmpDir= new File(System.getProperty("java.io.tmpdir"));
			}
*/
			if (refname!= null) 
				settings.refFile= new File(refname);
			if (proname!= null) 
				settings.proFile= new File(proname);
			if (frgname!= null) 
				settings.frgFile= new File(frgname);
			if (seqname!= null) 
				settings.seqFile= new File(seqname);
			if (errname!= null) {
				File fx= new File(errname);
				if (fx.exists())
					settings.errFile= fx;
			}
			if (tmpname!= null) {
				File fx= new File(tmpname);
				if (fx.exists()&& fx.isDirectory())
					settings.tmpDir= fx;
			}
			if (genname!= null) {
				File fx= new File(genname);
				if (fx.exists()&& fx.isDirectory())
					settings.genDir= fx;
			}
			
			return settings;
			
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
		
	}

	/**
	 * Fills all not inited attributes with default values
	 * @param settings
	 * @return
	 */
	public static String fillDefaults(FluxSimulatorSettings settings) {
		
		StringBuilder sb= new StringBuilder();
		if (Double.isNaN(settings.expDistrP1)) {
			settings.expDistrP1= DEF_EXP_DISTR_K;
			sb.append("\n\t"+ PAR_EXPR_K+"   "+ DEF_EXP_DISTR_K);
		}
		if (Double.isNaN(settings.expDistrP2)) {
			settings.expDistrP2= DEF_EXP_DISTR_X0;
			sb.append("\n\t"+ PAR_EXPR_X0+"   "+ DEF_EXP_DISTR_X0);
		}
		if (Double.isNaN(settings.decDistrP1)) {
			settings.decDistrP1= DEF_EXP_DECAY_P1;
			sb.append("\n\t"+ PAR_DEC_P1+"   "+ DEF_EXP_DECAY_P1);
		}


		if (Double.isNaN(settings.polyAshape)) {
			settings.polyAshape= DEF_POLYA_SHAPE;
			sb.append("\n\t"+ PAR_POLYA_SHAPE+"   "+ DEF_POLYA_SHAPE);
		}
		if (Double.isNaN(settings.polyAscale)) {
			settings.polyAscale= DEF_POLYA_SCALE;
			sb.append("\n\t"+ PAR_POLYA_SCALE+"   "+ DEF_POLYA_SCALE);
		}
		if (Double.isNaN(settings.tssMean)) {
			settings.tssMean= DEF_TSS_MEAN;
			sb.append("\n\t"+ PAR_TSS_MEAN+"   "+ DEF_TSS_MEAN);
		}

		if (settings.rtMode== null) {
			//settings.rtMode= DEF_RT_MODE;
			//sb.append("\n\t"+ PAR_RT_MODE+"   "+ DEF_RT_MODE);
		}
		if (settings.minRTLen< 0) {
			settings.minRTLen= DEF_RT_MIN_LEN;
			sb.append("\n\t"+ PAR_RT_MIN+"   "+ DEF_RT_MIN_LEN);
		}
		if (settings.maxRTLen< 0) {
			settings.maxRTLen= DEF_RT_MAX_LEN;
			sb.append("\n\t"+ PAR_RT_MAX+"   "+ DEF_RT_MAX_LEN);
		}
		if (settings.fragment) {
			if (settings.fragMode== null) {
				settings.fragMode= DEF_FRAG_MODE;
				sb.append("\n\t"+ PAR_FRAG_MODE+"   "+ DEF_FRAG_MODE);
			}
			if (Double.isNaN(settings.lambda)) {
				settings.lambda= DEF_FRAG_LAMBDA;
				sb.append("\n\t"+ PAR_FRAG_LAMBDA+"   "+ DEF_FRAG_LAMBDA);
			}
			if (Double.isNaN(settings.sigma)) {
				settings.sigma= DEF_FRAG_SIGMA;
				sb.append("\n\t"+ PAR_FRAG_SIGMA+"   "+ DEF_FRAG_SIGMA);
			}
			if (Double.isNaN(settings.thold)) {
				settings.thold= DEF_FRAG_THRESHOLD;
				sb.append("\n\t"+ PAR_FRAG_THRESHOLD+"   "+ DEF_FRAG_THRESHOLD);
			}
		}
		if (settings.filter) {
			if (settings.filtMin< 0) {
				settings.filtMin= DEF_FILT_MIN;
				sb.append("\n\t"+ PAR_FILT_MIN+"   "+ DEF_FILT_MIN);
			}
			if (settings.filtMax< 0) {
				settings.filtMax= DEF_FILT_MAX;
				sb.append("\n\t"+ PAR_FILT_MAX+"   "+ DEF_FILT_MAX);
			}
		}

		if (settings.readNr< 0) {
			settings.readNr= DEF_READ_NR;
			sb.append("\n\t"+ PAR_SEQ_READ_NUMBER+"   "+ DEF_READ_NR);
		}
		if (settings.readLength< 0) {
			settings.readLength= DEF_READ_LENGTH;
			sb.append("\n\t"+ PAR_SEQ_READ_LENGTH+"   "+ DEF_READ_LENGTH);
		}
		if (settings.filtMax< 0) {
			settings.filtMax= DEF_FILT_MAX;
			sb.append("\n\t"+ PAR_FILT_MAX+"   "+ DEF_FILT_MAX);
		}
		//settings.pairedEnd= DEF_PAIRED_END;
		if (settings.errFile!= null) {
			if (Float.isNaN(settings.qthold)) {
				settings.qthold= DEF_QTHOLD;
				sb.append("\n\t"+ PAR_SEQ_QTHOLD+"   "+ DEF_QTHOLD);
			}
		}
		
		return sb.toString();
	}
	
	
	public static final int DEF_READ_LENGTH= 36, DEF_RT_RANDOM= 2000, DEF_RT_FALLOFF= 2000;
	public static final float DEF_QTHOLD= 33f;
	public static final long DEF_READ_NR= 5000000;
	public static final String DEF_RT_MODE= PAR_RT_MODE_RANDOM;
	public static final int DEF_RT_MIN_LEN= 500, DEF_RT_MAX_LEN= 5500;
	public static final String DEF_FRAG_MODE= PAR_FRAG_MODE_PHYS;
	public static final long DEF_NB_CELLS= 50, DEF_NB_MOLECULES= 5000000;
	public static final double DEF_TSS_MEAN= 25d, DEF_POLYA_SHAPE= 2d, DEF_POLYA_SCALE= 300d;
	public static final boolean DEF_POLYA_VAR= true, DEF_TSS_VAR= true;
	public static final double DEF_EXP_DISTR_K= -0.6, DEF_EXP_DISTR_X0= 50000000/*0.01523906519187669635*/, DEF_EXP_DECAY_P1= 9500;
	public static final double DEF_FRAG_SIGMA= 0.05, DEF_FRAG_LAMBDA= 900, DEF_FRG_LAMBDA= 700, DEF_FRAG_THRESHOLD= 0.1;
	public static final String DEF_SFX_PAR= ".par", DEF_SFX_PRO= ".pro", DEF_SFX_LIB= ".lib", DEF_SFX_SEQ= ".bed", DEF_SFX_QFASTA= ".qfasta", DEF_SFX_ERR= ".err", DEF_SFX_GEM= ".map";
	public static final int DEF_FILT_MIN = 200, DEF_FILT_MAX= 250;
	public static final int SF_CELL_HI = 50; 	// 500;
	public static final int SF_CELL_MED = 15;
	public static final int SF_SHORT = 1500;
	public static final int SF_MEDIUM = 3000;
	public static final boolean DEF_FRAG_B4_RT= false, DEF_PAIRED_END= false, DEF_FASTQ= false;
	
	public static FluxSimulatorSettings createDefaults() {
		FluxSimulatorSettings settings= new FluxSimulatorSettings();
		
		settings.readLength= DEF_READ_LENGTH;
		
		settings.nbCells= DEF_NB_CELLS;
		settings.nbMolecules= DEF_NB_MOLECULES;
		
		settings.expDistrP1= DEF_EXP_DISTR_K;
		settings.expDistrP2= DEF_EXP_DISTR_X0;
		settings.decDistrP1= DEF_EXP_DECAY_P1;
		
		settings.tssMean= DEF_TSS_MEAN;
		settings.polyAshape= DEF_POLYA_SHAPE;
		settings.polyAscale= DEF_POLYA_SCALE;
		
		settings.rtMode= DEF_RT_MODE;
		settings.minRTLen= DEF_RT_MIN_LEN;
		settings.maxRTLen= DEF_RT_MAX_LEN;
		
		settings.sigma= DEF_FRAG_SIGMA;
		settings.lambda= DEF_FRAG_LAMBDA;
		settings.thold= DEF_FRAG_THRESHOLD;
		settings.fragMode= DEF_FRAG_MODE;
		
		settings.filtMin= DEF_FILT_MIN;
		settings.filtMax= DEF_FILT_MAX;
		
		settings.readNr= DEF_READ_NR;
		settings.readLength= DEF_READ_LENGTH;
		settings.qthold= DEF_QTHOLD;
		settings.pairedEnd= DEF_PAIRED_END;
		
		return settings;
	}
	

	
	File parFile, refFile, proFile, frgFile, seqFile, genDir, errFile;
	File filePWMfragSense, filePWMfragAsense, filePWMrtSense, filePWMrtAsense;
	File fileFilterDistr;
	float qthold= Float.NaN;
	int readLength= -1;
	long readNr= -1;
	double decDistrP1 = Double.NaN;	// was: 0
	double expDistrP1 = Double.NaN;
	double expDistrP2 = Double.NaN;
	long nbCells= -1;
	long nbMolecules= -1;
	int minRTLen= -1, maxRTLen= -1;
	File tmpDir;
	double thold = Double.NaN;
	double sigma = Double.NaN;
	double lambda= Double.NaN;
	double tssMean= Double.NaN, polyAshape= Double.NaN, polyAscale= Double.NaN;
	String rtMode= null, fragMode= null;
	int filtMin= -1, filtMax= -1;
	String name;
	Profiler profiler= null;
	boolean fragB4RT= DEF_FRAG_B4_RT;
	boolean pairedEnd= DEF_PAIRED_END, fastQ= DEF_FASTQ;
	boolean fragment= false, filter= false;	
	boolean loadCoding= true, loadNoncoding= true;
	int maxThreads= 1;
	byte bpDistr= DISTR_NI;
	float[] bpDistrPar= null;
	/**
	 * 
	 */
	File bpMotif= null;

	public int getMaxThreads() {
		return maxThreads;
	}

	public void setMaxThreads(int maxThreads) {
		this.maxThreads = maxThreads;
	}

	public File getRefFile() {
		return refFile;
	}


	public File getProFile() {
		return proFile;
	}


	public File getFrgFile() {
		return frgFile;
	}


	public File getSeqFile() {
		return seqFile;
	}


	public int getReadLength() {
		return readLength;
	}


	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}


	public double getDecDistrP1() {
		return decDistrP1;
	}


	public void setDecDistrP1(double decDistrP1) {
		this.decDistrP1 = decDistrP1;
	}


	public double getExpDistrP1() {
		return expDistrP1;
	}


	public void setExpDistrP1(double expDistrP1) {
		this.expDistrP1 = expDistrP1;
	}


	public double getExpDistrP2() {
		return expDistrP2;
	}


	public void setExpDistrP2(double expDistrP2) {
		this.expDistrP2 = expDistrP2;
	}


	public long getNbCells() {
		return nbCells;
	}


	public void setNbCells(long nbCells) {
		this.nbCells = nbCells;
	}


	public long getNbMolecules() {
		return nbMolecules;
	}


	public void setNbMolecules(long nbMolecules) {
		this.nbMolecules = nbMolecules;
	}


	public File getTmpDir() {
		return tmpDir;
	}

	public void save() {
		if (getParFile()!= null) {
			writeParfile(getParFile(), this);
		}
	}

	public void setTmpDir(File tmpDir) {
		this.tmpDir = tmpDir;
	}


	public double getThold() {
		return thold;
	}


	public void setThold(double thold) {
		this.thold = thold;
	}


	/**
	 * @deprecated not used
	 * @return
	 */
	public double getSigma() {
		return sigma;
	}


	public void setSigma(double sigma) {
		this.sigma = sigma;
	}


	public double getLambda() {
		return lambda;
	}


	public void setLambda(double lambda) {
		this.lambda = lambda;
	}


	public void setRefFile(File refFile) {
		this.refFile = refFile;
	}


	public void setProFile(File proFile) {
		this.proFile = proFile;
	}


	public void setFrgFile(File frgFile) {
		this.frgFile = frgFile;
	}


	public void setSeqFile(File seqFile) {
		this.seqFile = seqFile;
	}


	public String getName() {
		return name;
	}


	public void setName(String name) {
		this.name = name;
	}


	public Profiler getProfiler() {
		return profiler;
	}


	public void setProfiler(Profiler profiler) {
		this.profiler = profiler;
	}


	public int getFiltMin() {
		return filtMin;
	}


	public void setFiltMin(int filtMin) {
		this.filtMin = filtMin;
	}


	public int getFiltMax() {
		return filtMax;
	}


	public void setFiltMax(int filtMax) {
		this.filtMax = filtMax;
	}


	public boolean isFragB4RT() {
		return fragB4RT;
	}


	public void setFragB4RT(boolean fragB4RT) {
		this.fragB4RT = fragB4RT;
	}


	public String getRtMode() {
		return rtMode;
	}

	public String getFragMode() {
		return fragMode;
	}

	public void setRtMode(String rtMode) {
		this.rtMode = rtMode;
	}
	
	public void setFragMode(String fragMode) {
		this.fragMode = fragMode;
	}


	public boolean isPairedEnd() {
		return pairedEnd;
	}


	public void setPairedEnd(boolean pairedEnd) {
		this.pairedEnd = pairedEnd;
	}


	public long getReadNr() {
		return readNr;
	}


	public void setReadNr(long readNr) {
		this.readNr = readNr;
	}

	public File getParFile() {
		return parFile;
	}

	public void setParFile(File parFile) {
		this.parFile = parFile;
	}

	public int getMinRTLen() {
		return minRTLen;
	}

	public void setMinRTLen(int minRTLen) {
		this.minRTLen = minRTLen;
	}

	public int getMaxRTLen() {
		return maxRTLen;
	}

	public void setMaxRTLen(int maxRTLen) {
		this.maxRTLen = maxRTLen;
	}

	public boolean isFragment() {
		return fragment;
	}

	public void setFragment(boolean fragment) {
		this.fragment = fragment;
	}

	public boolean isFilter() {
		return filter;
	}

	public void setFilter(boolean filter) {
		this.filter = filter;
	}

	public boolean isLoadCoding() {
		return loadCoding;
	}

	public void setLoadCoding(boolean loadCoding) {
		this.loadCoding = loadCoding;
	}

	public boolean isLoadNoncoding() {
		return loadNoncoding;
	}

	public void setLoadNoncoding(boolean loadNoncoding) {
		this.loadNoncoding = loadNoncoding;
	}

	public File getGenDir() {
		return genDir;
	}

	public void setGenDir(File genDir) {
		this.genDir = genDir;
	}

	public File getErrFile() {
		return errFile;
	}

	public void setErrFile(File errFile) {
		this.errFile = errFile;
	}

	public float getQthold() {
		return qthold;
	}

	public void setQthold(float qthold) {
		this.qthold = qthold;
	}

	public boolean isFastQ() {
		return fastQ;
	}

	public void setFastQ(boolean fastQ) {
		this.fastQ = fastQ;
	}

	public double getTssMean() {
		return tssMean;
	}

	public void setTssMean(double tssMean) {
		this.tssMean = tssMean;
	}

	public double getPolyAshape() {
		return polyAshape;
	}

	public void setPolyAshape(double polyAshape) {
		this.polyAshape = polyAshape;
	}

	public double getPolyAscale() {
		return polyAscale;
	}

	public void setPolyAscale(double polyAscale) {
		this.polyAscale = polyAscale;
	}
	
	@Override
	public Object clone() throws CloneNotSupportedException {
		FluxSimulatorSettings settings= new FluxSimulatorSettings();
		settings.parFile= parFile;
		settings.decDistrP1= decDistrP1;
		settings.errFile= errFile;
		settings.expDistrP1= expDistrP1;
		settings.expDistrP2= expDistrP2;
		settings.fastQ= fastQ;
		settings.filter= filter;
		settings.filtMax= filtMax;
		settings.filtMin= filtMin;
		settings.fragB4RT= fragB4RT;
		settings.fragment= fragment;
		settings.fragMode= fragMode;
		settings.frgFile= frgFile;
		settings.genDir= genDir;
		settings.lambda= lambda;
		settings.loadCoding= loadCoding;
		settings.loadNoncoding= loadNoncoding;
		settings.maxRTLen= maxRTLen;
		settings.minRTLen= minRTLen;
		settings.name= name;
		settings.nbCells= nbCells;
		settings.nbMolecules= nbMolecules;
		settings.pairedEnd= pairedEnd;
		settings.parFile= parFile;
		settings.polyAscale= polyAscale;
		settings.polyAshape= polyAshape;
		settings.proFile= proFile;
		settings.profiler= profiler;
		settings.qthold= qthold;
		settings.readLength= readLength;
		settings.readNr= readNr;
		settings.refFile= refFile;
		settings.rtMode= rtMode;
		settings.seqFile= seqFile;
		settings.sigma= sigma;
		settings.thold= thold;
		settings.tmpDir= tmpDir;
		settings.tssMean= tssMean;
		
		return settings;
	}

	public boolean isTssVar() {
		return (!(Double.isNaN(tssMean)));
	}

	public boolean isPolyAVar() {
		return (!(Double.isNaN(polyAshape)|| Double.isNaN(polyAscale)));
	}

	public void setTssVar(boolean tssVar) {
		if (!tssVar)
			this.tssMean = Double.NaN;
	}

	public void setPolyAVar(boolean polyAVar) {
		if (!polyAVar) {
			this.polyAscale = Double.NaN;
			this.polyAshape = Double.NaN;
		}
	}

	public File getFilePWMfragSense() {
		return filePWMfragSense;
	}
 
	public File getFilePWMrtSense() {
		return filePWMrtSense;
	}

	public byte getBpDistr() {
		return bpDistr;
	}

	public void setBpDistr(byte bpDistr) {
		this.bpDistr = bpDistr;
	}

	public float[] getBpDistrPar() {
		return bpDistrPar;
	}

	public void setBpDistrPar(float[] bpDistrPar) {
		this.bpDistrPar = bpDistrPar;
	}

	public File getFilePWMfragAsense() {
		return filePWMfragAsense;
	}

	public void setFilePWMfragAsense(File filePWMfragAsense) {
		this.filePWMfragAsense = filePWMfragAsense;
	}

	public File getFilePWMrtAsense() {
		return filePWMrtAsense;
	}

	public void setFilePWMrtAsense(File filePWMrtAsense) {
		this.filePWMrtAsense = filePWMrtAsense;
	}

	public void setFilePWMfragSense(File filePWMfragSense) {
		this.filePWMfragSense = filePWMfragSense;
	}

	public void setFilePWMrtSense(File filePWMrtSense) {
		this.filePWMrtSense = filePWMrtSense;
	}

	public File getFileFilterDistr() {
		return fileFilterDistr;
	}

	public void setFileFilterDistr(File fileFilterDistr) {
		this.fileFilterDistr = fileFilterDistr;
	}

}
