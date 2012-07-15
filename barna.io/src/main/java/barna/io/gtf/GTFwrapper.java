/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.io.gtf;

import barna.commons.ByteArrayCharSequence;
import barna.commons.io.IOHandler;
import barna.commons.io.IOHandlerFactory;
import barna.commons.log.Log;
import barna.commons.utils.ArrayUtils;
import barna.commons.utils.StringUtils;
import barna.io.AbstractFileIOWrapper;
import barna.io.AnnotationWrapper;
import barna.io.FileHelper;
import barna.io.Sorter;
import barna.model.*;
import barna.model.commons.IntVector;
import barna.model.commons.MyArrayHashMap;
import barna.model.constants.Constants;
import barna.model.gff.GFFObject;

import java.io.*;
import java.text.Collator;
import java.util.*;
import java.util.concurrent.Future;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//import genome.tools.MyArray;
//import org.apache.commons.collections.BidiMap;


/**
 * <u>07/11/21</u>: bugfix in read() for resetting last line at chromosome change
 * <u>07/11/22</u>: can now handle comment lines ('# ...').
 * <u>07/11/30</u>: bugfix for not adding last gene twice when abort condition in read()
 * <u>08/08/18</u>: bugfix in sweepChromosome(): whole subsequent chr was skipped when the 
 * 										chr that was intented to be skipped consisted only of 1 line.
 * <u>08/11/14</u>: bugfix for windows line separators, <code>bytesRead</code> adapted
 * TODO: allow for custom gene clustering
 * 
 * @author micha
 *
 */
public class GTFwrapper extends AbstractFileIOWrapper implements AnnotationWrapper {

    /**
     * Line splitter pattern
     */
    static final Pattern SPLITTER_PATTERN = Pattern.compile("\\s");
    /**
     * Pattern to findTranscriptID the transcript ID
     */
    static final Pattern TRANSCRIPTID_PATTERN = Pattern.compile(GFFObject.TRANSCRIPT_ID_TAG);

    /**
     * Helper to compare entries based on their global position
     */
    class GlobalTranscriptPositionComparator implements Comparator<String>{
        Map<byte[], Integer> map;
        int[] fieldNrs;
        //ByteArrayCharSequence cs;

        GlobalTranscriptPositionComparator(Map<byte[], Integer> map, int[] fieldNrs) {
            this.map = map;
            this.fieldNrs = fieldNrs;
        }

        public int compare(String o1, String o2) {
            ByteArrayCharSequence cs = new ByteArrayCharSequence(o1);
            ByteArrayCharSequence[] fields= ByteArrayCharSequence.split(cs, SPLITTER_PATTERN, fieldNrs);	// 1,4,-1
            byte[] key1= ByteArrayCharSequence.merge(fields[1], fields[0]);
            Integer mapValue1 = map.get(key1);
            int v1= Math.abs(mapValue1);

            cs.clear();
            cs.append(o2);
            fields= ByteArrayCharSequence.split(cs, SPLITTER_PATTERN, fieldNrs);	// 1,4,-1
            byte[] key2= ByteArrayCharSequence.merge(fields[1], fields[0]);
            Integer mapValue2 = map.get(key2);
            int v2= Math.abs(mapValue2);
            int result = v1 - v2;
            return result;
        }
    }
	
	static boolean sort = false;

	static String usage = "GTFChrReader [options] <inputFile>\n\n"
			+ "where options may be\n"
			+ "-sort\t sort according to chr, transcriptID, start\n\n"
			+ "micha, 2007";

	public static final String NORMED_FILE_TAG = "norman";

	public static final int MAX_GTF_LINE_LENGTH = 1000;

	public static final String[] CHROMOSOME_FILETER_MITO= new String[] {
		"^CHRM.*$", "^M$", "^MT$"
	};
	
	public static final String[] CHROMOSOME_FILETER_UNKNOWN= new String[] {
		"^U$", "^UN$", "^CHRU$", "^CHRUN$", ".*UNKNOWN.*" 
	};
	
	public static final String[] CHROMOSOME_FILETER_HAPLO_HETERO= new String[] {
		".*_HAP.*", // chr5_h2_hap1
		"^\\d{1,}H$", // droso, hetero?
		"^CHR\\d{1,}H$", "^X{1}H$", "^CHRXH$", "^Y{1}H$", "^CHRYH$"
	};
	
	public static final String[] CHROMOSOME_FILETER_RANDOM= new String[] {
		".*_RANDOM.*"
	};	
	
	public static String[] compose(String[] a, String[] b) {
		String[] c= new String[a.length+ b.length];
		for (int i = 0; i < a.length; i++) 
			c[i]= a[i];
		for (int i = 0; i < b.length; i++) 
			c[a.length+ i]= b[i];
		return c;
	}
	
	
	public static String[] DEFAULT_CHROMOSOME_FILTER = new String[] {
			"^NT_.+$" 
	}; // default chromosomes filtered off
	static {
		DEFAULT_CHROMOSOME_FILTER= compose(CHROMOSOME_FILETER_MITO, DEFAULT_CHROMOSOME_FILTER);
		DEFAULT_CHROMOSOME_FILTER= compose(CHROMOSOME_FILETER_UNKNOWN, DEFAULT_CHROMOSOME_FILTER);
		DEFAULT_CHROMOSOME_FILTER= compose(CHROMOSOME_FILETER_HAPLO_HETERO, DEFAULT_CHROMOSOME_FILTER);
		DEFAULT_CHROMOSOME_FILTER= compose(CHROMOSOME_FILETER_RANDOM, DEFAULT_CHROMOSOME_FILTER);
	}

	protected InputStream inputStream = null;

	Species species = null;

	int limitGTFObs= 0;
	// long bytesRead = 0, size = 0l; // from superclass

	int readAheadLimit= -1, readAheadTranscripts= 50000;

	Vector<String> readChr = null, skippedChr = new Vector<String>(), skippedFeatures = new Vector<String>(), skippedTranscripts = new Vector<String>();
	int readTranscripts = 0, readExons = 0, readGenes = 0, readChrs = 0;
	int skippedObjects = 0;

	HashMap<String,Integer> filtSomeIDs = null;
	boolean[] filtSomeIDSuccess = null;
	DirectedRegion[] filtRegs = null;
	String[] filtChrIDs = null, noIDs = DEFAULT_CHROMOSOME_FILTER, filtGeneIDs = null, filtTrptIDs = null, 
		readFeatures = new String[] { "exon", "CDS", "start_codon", "stop_codon" }, allowSources= null;

	boolean silent = false, stars= true, outputWarnings= false, chromosomeWise = true, strandWise= true, geneWise = true,
			clusterGenes = true;
	boolean clustered = false;
	boolean readGTF = false, readGene = true, printStatistics= true;
	
	int nrGenes= 0, nrTranscripts= 0, nrExons= 0;
	int nrLinesRead = 0;
	
	/**
	 * reuse the file handle.
	 */
	boolean reuse= false;
	String lastLine= null;
	
	Gene[] genes = null;

	GFFObject[] gtfObj;
	
	boolean keepOriginalLines= false;
	public boolean isKeepOriginalLines() {
		return keepOriginalLines;
	}

	public void setKeepOriginalLines(boolean keepOriginalLines) {
		this.keepOriginalLines = keepOriginalLines;
	}


	Vector<String> vLines= null;

	public Vector<String> getVLines() {
		return vLines;
	}

	public static void main(String[] args) {

        System.out.println(SPLITTER_PATTERN);
	}

	void initSpecies() {
		StringTokenizer st = new StringTokenizer(getInputFile().getName(), "_.");
		String speName = null, annotationVersion = null, genomeVersion = null;
		int speNb = -1;
		while (st.hasMoreTokens()) {
			String s = st.nextToken();
			int tmpSpNb = Species.getSpeciesNumber(s);
			if (tmpSpNb >= 0) {
				if (speName == null)
					speName = s;
				else {
					if (!speName.equalsIgnoreCase(s))
						System.err
								.println("Conflicting info on species in file name.");
				}
				if (speNb >= 0 && tmpSpNb != speNb)
					System.err.println("Conflicting info in file name: " + s
							+ " <> " + Species.SP_NAMES_COMMON[speNb]);
				continue;
			}

			tmpSpNb = Species.getSpeciesNumberForGenomeVersion(s);
			if (tmpSpNb >= 0) {
				if (speNb >= 0 && tmpSpNb != speNb)
					System.err.println("Conflicting info in file name: " + s
							+ " <> " + Species.SP_NAMES_COMMON[speNb]);
				speNb = tmpSpNb;
				speName = Species.SP_NAMES_COMMON[speNb];
				if (genomeVersion == null)
					genomeVersion = s;
				else {
					if (Species.getGenomeVerNb(genomeVersion) != Species
							.getGenomeVerNb(s))
						System.err
								.println("Conflicting info on genome in file name: "
										+ genomeVersion + " <> " + s);
				}
				continue;
			}
		}

		if (speName == null) {
			if (!silent)
				System.out
						.println("No valid species name found in file name, guessing \'human\'.");
			species = new Species("human");
		} else
			species = new Species(speName);

		if (genomeVersion == null) {
			genomeVersion = species.getDefaultGenomeVersion();
			if (!silent)
				System.out
						.println("No valid build found in file name, guessing "
								+ genomeVersion + ".");
			species.setGenomeVersion(genomeVersion);
		} else
			species.setGenomeVersion(genomeVersion);

		if (annotationVersion != null)
			species.setAnnotationVersion(annotationVersion);
	}

	/**
	 * Default constructor requires a <code>File</code> instance, 
	 * because of sweeping while reading and 2-pass reads of the
	 * source during sorting require to open multiple streams for
	 * reading from the source.<br>
	 * <b>Note</b>: For the same reasoning the wrapper does not 
	 * permit an <code>InputStream</code>-based constructor.
	 * @param inputFile
	 */
	public GTFwrapper(File inputFile) {
		super(inputFile);
		readChr = new Vector();
		// initSpecies();
	}

	public GTFwrapper(String absFName) {
		//super(absFName);
		this(new File(absFName));
	}



	boolean overlap(Vector trans1, Vector trans2) {
		// for (int i = 0; i < trans1.size(); i++) {
		// France f1= (France) trans1.elementAt(i);
		// if (!f1.isExon())
		// continue;
		// for (int j = 0; j < trans2.size(); j++) {
		// France f2= (France) trans2.elementAt(j);
		// if (!f2.isExon())
		// continue;
		// DefaultDirectedRegion d1= new DefaultDirectedRegion(f1.getStrand(),
		// f1.getStart(), f1.getEnd());
		// DefaultDirectedRegion d2= new DefaultDirectedRegion(f1.getStrand(),
		// f2.getStart(), f2.getEnd());
		// if (d1.overlaps(d2))
		// return true;
		// }
		// }
		return false;
	}

	static HashMap getChromosomes(HashMap transGTF) {
		Iterator iter = transGTF.values().iterator();
		HashMap chrHash = new HashMap();
		while (iter.hasNext()) {
			Vector gtfsVec = (Vector) iter.next();
			GFFObject o = (GFFObject) (gtfsVec).elementAt(0);
			HashMap tHash = (HashMap) chrHash.remove(o.getChromosome());
			if (tHash == null)
				tHash = new HashMap();
			Vector v = (Vector) tHash.remove(o.getTranscriptID());
			if (v == null)
				v = new Vector();
			for (int i = 0; i < gtfsVec.size(); i++)
				v.add(gtfsVec.elementAt(i));
			tHash.put(o.getTranscriptID(), v);
			chrHash.put(o.getChromosome(), tHash);
		}
		return chrHash;
	}

	private void checkMegaClusters(Transcript[] trans) {
		java.util.Arrays.sort(trans, new DirectedRegion.PositionComparator()); // sort
																				// ascending
		for (int i = trans.length - 1; i >= 0; --i) {
			System.out.println(trans[i].getTranscriptID() + "\t"
					+ trans[i].getLength() + "\t" + trans[i].getExons().length);
		}
	}

	public boolean isEOF() {
		return (bytesRead == getInputFile().length());
	}

	/**
	 * depends on the order of gtf-objects from the input for creating
	 * transcripts to create from the outside
	 * 
	 * @param encode
	 * @return
	 */
	public static Species assemble(GFFObject[] gtfObs, String speName,
			String buildVersion) {

		Species spe = new Species(speName);
		// Gene[] ge= assemble(gtfObs);
		// spe.setAnnotationVersion(buildVersion);
		// for (int i = 0; i < ge.length; i++) {
		// spe.addGene(ge[i]);
		// ge[i].setSpecies(spe);
		// }

		return spe;
	}

	public String getLastReadChr() {
		if (readChr == null)
			return null;
		return (String) readChr.elementAt(readChr.size() - 1);
	}

	static protected Transcript[][] clusterTranscripts(DirectedRegion[] regions) {

		int maxPlus = Integer.MIN_VALUE, maxMin = Integer.MIN_VALUE;
		Vector clusters = new Vector();
		Vector vPlus = null, vMinus = null;
		for (int i = 0; i < regions.length; i++) {

			DirectedRegion r = regions[i];
			if (regions[i].isForward()) {
				if (regions[i].getStart() > maxPlus) {
					if (vPlus != null)
						clusters.add(vPlus);
					vPlus = new Vector();
				}
				vPlus.add(regions[i]);
				if (regions[i].getEnd() > maxPlus)
					maxPlus = regions[i].getEnd();
			} else {
				if (Math.abs(regions[i].getStart()) > maxMin) {
					if (vMinus != null)
						clusters.add(vMinus);
					vMinus = new Vector();
				}
				vMinus.add(regions[i]);
				if (Math.abs(regions[i].getEnd()) > maxMin)
					maxMin = Math.abs(regions[i].getEnd());
			}
		}

		if (vPlus != null)
			clusters.add(vPlus);
		if (vMinus != null)
			clusters.add(vMinus);

		System.gc();

		return (Transcript[][]) ArrayUtils.toField(clusters);
	}

	static protected HashMap getGroups(String id, GFFObject[] obj) {

		HashMap hash = new HashMap();
		for (int i = 0; i < obj.length; i++) {
			String s = obj[i].getAttribute(id);
			Vector tAttrib = (Vector) hash.get(s);
			if (tAttrib == null) {
				tAttrib = new Vector();
				hash.put(s, tAttrib);
			}
			tAttrib.add(obj[i]);
		}

		return hash;
	}

	public String[] getGeneIDs(String[] protIDs) {
		Vector v = new Vector();
		Pattern[] patterns = new Pattern[protIDs.length];
		for (int i = 0; i < protIDs.length; i++)
			patterns[i] = Pattern.compile(protIDs[i]); // case sensitive
														// .toUpperCase()
		try {
			BufferedReader buffy;
			if (getInputFile() != null)
				buffy = new BufferedReader(getReader());
			else
				buffy = new BufferedReader(new InputStreamReader(inputStream));
			String line;
			int lineCtr = 0;
			Vector gtfVec = new Vector();
			System.out.println("Scanning for geneIDs ");
			long t0 = System.currentTimeMillis();
			while (buffy.ready()) {
				lineCtr++;
				if (lineCtr % 1000 == 0) {
					System.out.print("*");
					System.out.flush();
				}
				line = buffy.readLine();
				int x = 0;
				for (int i = 0; i < patterns.length; i++) {
					// Matcher m= patterns[i].matcher(line); // case sensitive
					// .toUpperCase()
					// if (m.find())
					if (line.contains(protIDs[i]))
						break;
				}
				if (x == protIDs.length)
					continue;

				String[] tokens = line.split("\t");
				for (int i = 8; i < tokens.length; i++) {
					if (tokens[i].equals(GFFObject.GENE_ID_TAG)) {
						String id = tokens[i + 1];
						if (id.startsWith("\"") || id.startsWith("\'"))
							id = id.substring(1, id.length() - 1);
						v.add(id);
						break;
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		System.out.println("\nPrimary list " + v.size());
		Comparator compi = Collator.getInstance();
		Vector vv = new Vector();
		for (int i = 0; i < v.size(); i++) {
			ArrayUtils.addUniqueSorted(vv, v.elementAt(i), compi);
		}
		System.out.println("Found " + vv.size() + " gene ids for "
				+ protIDs.length + " proteins.");
		return (String[]) ArrayUtils.toField(vv);
	}

	public void reset() {
		bytesRead = 0l;
		nrLinesRead = 0;
		nrGenes= 0;
		nrTranscripts= 0;
		nrExons= 0;
		vLines= null;
		genes= null;
		gtfObj= null;
		readChr = new Vector();
		skippedChr = new Vector();
		skippedTranscripts = new Vector();
		skippedTranscripts = new Vector();
		buffy= null;	// has to be inited
		// TODO: region and ID-filtering?? or do it method-level..
	}

	public void sweepToChromosome(String chrom) {
		BufferedReader buffy = null;
		long size = 0l;
	
		if (getInputFile() != null) {
			try {
				size = getInputSize();
				inputStream = getInputStream();
				inputStream.skip(bytesRead);
				buffy = new BufferedReader(new InputStreamReader(inputStream));
				// buffy= new BufferedReader(new FileReader(file));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	
		// reset(); // start at the beginning
		// int lastPerc= (int) ((bytesRead* 10d)/ size);
		try {
			// raf.seek(bytesRead);
			boolean reset = false;
			while (true) {
				String line = buffy.readLine();
				if (line == null) {
					if (reset) {
						System.out.println("WARNING: chromosome " + chrom
								+ " not found for sweeping in "+ getInputFile().getName());
						break;
					} else {
						buffy = new BufferedReader(getReader());
						reset();
						line = buffy.readLine();
						reset = true; // only once
					}
	
				}
	
				String[] tokens = line.split("\t");
				String chr = null;
				if (tokens[0].length() <= 3)
					chr = "chr" + tokens[0].trim();
				else
					chr = tokens[0].trim();
				if (chr.equalsIgnoreCase(chrom))
					break;
	
				bytesRead += (line.length() + 1); // else
				++nrLinesRead;
				int perc = (int) ((bytesRead * 10d) / size);
				// if (perc> lastPerc&& !silent) {
				// ++lastPerc;
				// System.out.print("*");
				// System.out.flush();
				// }
			}
			buffy.close();
	
		} catch (Exception e) {
			e.printStackTrace();
		}
		return;
	
	}

	public void write() {
		// TODO implement/copy
	}
	
	public int[] getTxPerLocus() {
		if (txPerLocus== null)
			return null;
		return txPerLocus.toIntArray();
	}
	
	public int[] getTxLengths() {
		if (txLengths== null)
			return null;
		return txLengths.toIntArray();
	}
	
	IntVector txPerLocus= null, txLengths= null;
	public void scanFile() {

        if(!silent && stars){
            Log.progressStart("scanning");
        }
		reset();
		
		setReadGTF(false);
		setReadGene(true);
		setReadAheadLimit(1);
		setReadAheadTranscripts(1); // keep mem low
		txPerLocus= new IntVector();
		txLengths= new IntVector();

		try {
			read();
			int chkGeneCnt= 0, chkTxCnt= 0, chkExonCnt= 0;
			while (getGenes()!= null) {				
				for (int i = 0; i < getGenes().length; i++) {
					++chkGeneCnt;
					txPerLocus.add(getGenes()[i].getTranscriptCount());
					for (int j = 0; j < getGenes()[i].getTranscripts().length; j++) {
						++chkTxCnt;
						chkExonCnt+= getGenes()[i].getTranscripts()[j].getExons().length;
						txLengths.add(getGenes()[i].getTranscripts()[j].getExonicLength());
					}
				}
				read();
			}
			//System.err.println("[CHECK] read "+chkGeneCnt+" genes, "+chkTxCnt+" tx, "+chkExonCnt+" exons.");
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		boolean b= close();
        if(!silent && stars){
            Log.progressFinish(StringUtils.OK, true);
        }

	}
	
	
	
	
	public void sweepToChromosomeStrand(String chrom, byte strand) {
		BufferedReader buffy = null;
		long size = 0l;
		if (getInputFile() != null) {
			try {
				if (lineSeparator== null)
					lineSeparator= FileHelper.guessFileSep(inputFile);
				size = inputFile.length();
				inputStream = new FileInputStream(inputFile);
				inputStream.skip(bytesRead);
				buffy = new BufferedReader(new InputStreamReader(inputStream));
				// buffy= new BufferedReader(new FileReader(file));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		// reset(); // start at the beginning
		// int lastPerc= (int) ((bytesRead* 10d)/ size);
		try {
			// raf.seek(bytesRead);
			boolean reset = false;
			while (true) {
				String line = buffy.readLine();
				if (line == null) {
					if (reset) {
						System.out.println("WARNING: chromosome " + chrom+ " strand "+ strand
								+ " not found for sweeping "+ getInputFile().getName());
						break;
					} else {
						buffy = new BufferedReader(getReader());
						reset();
						line = buffy.readLine();
						reset = true; // only once
					}

				}

				String[] tokens = line.split("\t");
				String chr = null;
				if (tokens[0].length() <= 3)
					chr = "chr" + tokens[0].trim();
				else
					chr = tokens[0].trim();
				
				if (chr.equalsIgnoreCase(chrom)&& (strand== 0|| strand== GFFObject.parseStrand(tokens[6]))) {
					break;
				}

				bytesRead += (line.length() + lineSeparator.length()); // else
				++nrLinesRead;
				int perc = (int) ((bytesRead * 10d) / size);
				// if (perc> lastPerc&& !silent) {
				// ++lastPerc;
				// System.out.print("*");
				// System.out.flush();
				// }
			}
			buffy.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return;

	}

	void skipToNextChromosome_raf(RandomAccessFile raf) {
		try {
			// raf.seek(bytesRead);
			String lastChrID = null;
			int lastPerc = 0;
			long t0 = System.currentTimeMillis();
			while (true) {
				int perc = (int) ((bytesRead * 10d) / raf.length());
				if (perc > lastPerc && !silent) {
					++lastPerc;
					System.out.print("*");
					System.out.flush();
				}
				String line = raf.readLine();
				if (line == null) {
					break;
				}
				// String[] tokens= line.split("\t");
				// String chr= tokens[0];
				StringTokenizer toki = new StringTokenizer(line, "\t");
				String chr = toki.nextToken();
				if (chr.length() <= 3)
					chr = "chr" + chr.trim();
				else
					chr = chr.trim();

				if (lastChrID == null) {
					lastChrID = chr;
					System.out.print("skip " + lastChrID + " -> ");
					System.out.flush();
				}
				if (!chr.equals(lastChrID)) {
					System.out.println(chr + " "
							+ ((System.currentTimeMillis() - t0) / 1000)
							+ " sec");
					break;
				}
				bytesRead += (line.length() + lineSeparator.length()); // else
				++nrLinesRead;
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	BufferedReader skipToNextChromosome(BufferedReader buffy, long size, String lastChrID) {
		try {
			// raf.seek(bytesRead);
			int lastPerc = (int) ((bytesRead * 10d) / size);
			// long t0= System.currentTimeMillis();
			while (true) {
                if(!silent && stars){
                    Log.progress(bytesRead, size);
                }

				// buffy.mark(MAX_GTF_LINE_LENGTH); // does not work as planned, see reset() below
				long saveBytesRead= bytesRead;
				String line = buffy.readLine();
				if (line == null) {
					break;
				}
				String[] tokens = line.split("\t");
				String chr = tokens[0];
				// StringTokenizer toki= new StringTokenizer(line, "\t");
				// String chr= toki.nextToken();
				if (chr.length() <= 3)
					chr = "chr" + chr.trim();
				else
					chr = chr.trim();

				if (lastChrID == null) {
					lastChrID = chr;
					// System.out.print("skip "+lastChrID+" -> ");
					// System.out.flush();
				}
				if (!chr.equals(lastChrID)) {
					
					// buffy.reset();
					/*
					 * java.io.IOException: Mark invalid
	at java.io.BufferedReader.reset(BufferedReader.java:485)
	at gphase.io.gtf.GTFReader.skipToNextChromosome(GTFReader.java:1134)
	at gphase.io.gtf.GTFReader.read(GTFReader.java:1533)
	at gphase.tools.IntronExonSeqRetriever.read(IntronExonSeqRetriever.java:86)
	at gphase.tools.IntronExonSeqRetriever.main(IntronExonSeqRetriever.java:23)
					 */
					if (reuse)
						lastLine= line;
					else {
						buffy.close();
						bytesRead= saveBytesRead;
						buffy = getBuffy();
					}
						
					return buffy;
				}
				bytesRead += (line.length() + lineSeparator.length()); // else
				++nrLinesRead;
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		return null;
	}

	public String getNextChromosome() {
		RandomAccessFile raf = null;
		long size = 0l;
		if (getInputFile() != null) {
			try {
				raf = new RandomAccessFile(getInputFile(), "r");
				size = raf.length();
			} catch (Exception e) {
				e.printStackTrace();
			}
		} else if (inputStream != null)
			System.err.println(this.getClass()
					+ " only supports input from Files.");

		Gene[] geneReg = null;
		Comparator startCompi = new AbstractRegion.StartComparator();
		try {
			raf.seek(bytesRead);
			String line = raf.readLine();
			if (line == null)
				return null;
			raf.close();

			String[] tokens = line.split("\t");
			if (tokens[0].length() <= 3)
				return "chr" + tokens[0];
			else
				return tokens[0];
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;

	}

	public boolean checkChr(Gene[] ge) {
		int posStrand = 0;
		int negStrand = 0;
		boolean ok = true;
		for (int i = 0; i < ge.length; i++) {
			if (ge[i].getStrand() > 0)
				++posStrand;
			else if (ge[i].getStrand() < 0)
				++negStrand;
			if (ge[i].getStart() == 0 || ge[i].getEnd() == 0
					|| ge[i].getChromosome() == null || ge[i].getStrand() == 0) {
				System.err.println("Gene not inited " + ge[i].getGeneID());
				ok = false;
			}
			Transcript[] trpt = ge[i].getTranscripts();
			if (trpt.length == 0 || trpt.length > 10) {
				System.err.println("Too many transcripts per gene "
						+ ge[i].getGeneID());
				ok = false;
			}
			Exon[] ex = ge[i].getExons();
			if (ex.length == 0 || ex.length > 100) {
				System.err.println("Too many exons per gene "
						+ ge[i].getGeneID());
				ok = false;
			}

			for (int j = 0; j < trpt.length; j++) {
				if (trpt[j].getStart() == 0 || trpt[j].getEnd() == 0
						|| trpt[j].getChromosome() == null
						|| trpt[j].getStrand() == 0) {
					System.err.println("Transcript not inited "
							+ trpt[j].getTranscriptID());
					ok = false;
				}

			}
			for (int j = 0; j < ex.length; j++) {
				if (ex[j].getStart() == 0 || ex[j].getEnd() == 0
						|| ex[j].getChromosome() == null
						|| ex[j].getStrand() == 0) {
					System.err.println("Exon not inited "
							+ ex[j].getTranscripts()[0].getTranscriptID());
					ok = false;
				}

			}
		}

		if (Math.min(((double) posStrand) / negStrand, ((double) negStrand)
				/ posStrand) > 0.75d) {
			System.err.println("Unequal strand distribution, pos " + posStrand
					+ ", neg " + negStrand);
			ok = false;
		}
		return ok;
	}

	void printReadStatistics() {
		if (printStatistics) {
			System.err.println("\nRead:\t" + readChrs + " chromosomes, "
					+ readGenes + " genes, " + readTranscripts + " transcripts, "
					+ readExons + " exons.");
			System.err.println("Skipped:\t" + skippedObjects+ " lines.");
			System.err.print("Chromosomes:\t");
			for (int i = 0; i < skippedChr.size(); i++)
				System.err.print(skippedChr.elementAt(i) + " ");
			System.err.print("\nTranscripts:\t");
			for (int i = 0; i < skippedTranscripts.size(); i++)
				System.err.print(skippedTranscripts.elementAt(i) + " ");
			System.err.print("\nFeatures:\t");
			for (int i = 0; i < skippedFeatures.size(); i++)
				System.err.print(skippedFeatures.elementAt(i) + " ");
			System.err.println();
		}
	}

	boolean checkFeature(String feature) {
		int x; // skip mRNA, gene...
		for (x = 0; readFeatures != null && x < readFeatures.length; x++) {
			if (readFeatures[x].equalsIgnoreCase(feature))
				break;
		}
		if (readFeatures != null && x == readFeatures.length) {
			ArrayUtils.addUnique(skippedFeatures, feature);
			return false;
		}
		return true;
	}

	boolean checkSource(String source) {
		int x; // skip mRNA, gene...
		for (x = 0; allowSources != null && x < allowSources.length; x++) {
			if (allowSources[x].equalsIgnoreCase(source))
				break;
		}
		if (allowSources != null && x == allowSources.length) 
			return false;		
		return true;
	}

	boolean checkChromosome(String chrID) {
		
		if (getInputFile().getName().contains("cow") && chrID.equals("30"))
			chrID = "X";
		
		// neg chromosome filtering
		if (noIDs != null) {
			int u;
			String chrIDu = chrID.toUpperCase();
			for (u = 0; u < noIDs.length; u++) {
				if (Pattern.matches(noIDs[u], chrIDu))
					break;
			}
			if (u < noIDs.length)
				return false;

		}
		
		// pos chromosome filtering
		if (filtChrIDs != null) {
			int u;
			for (u = 0; u < filtChrIDs.length; u++) {
				if (filtChrIDs[u].equalsIgnoreCase(chrID))
					break;
			}
			if (u == filtChrIDs.length)
				return false;
		}

		for (int i = 0; readChr != null && i < readChr.size(); i++)
			if (readChr.elementAt(i).equals(chrID)) {
				if (!silent) {
					System.err.println("Line " + nrLinesRead+ ": chromosome " + chrID+ " already read!");
				}
				return false;
			}

		return true;

	}

	boolean checkOverlapRegion(Transcript trpt) {
		if (filtRegs != null) {
			int i;
			for (i = 0; i < filtRegs.length; i++)
				if (filtRegs[i].overlaps(trpt))
					break;
			if (i == filtRegs.length)
				return false;
		}
		return true;
	}

	boolean checkObject(GFFObject obj)  {
		if (obj.getStart()!= GFFObject.POSITION_UNDEFINED&& 
				obj.getEnd()!= GFFObject.POSITION_UNDEFINED&&
				obj.getStrand()!= GFFObject.STRAND_UNDEFINED&&
				checkFeature(obj.getFeature())&&
				checkSource(obj.getSource())&&
				checkAttributes(obj))
			return true;
		return false;
	}
	
	boolean checkAttributes(GFFObject obj) {
		if (filtSomeIDs != null) {
			
			Iterator<String> iter= obj.getAttributes().values().iterator();
			while (iter.hasNext()) {
				String[] t = iter.next().split(" ");
				for (int j = 0; j < t.length; j++) {
					Integer in = (Integer) filtSomeIDs.get(t[j]);
					if (in == null)
						continue; // ID not in filter list
					filtSomeIDs.put(t[j], new Integer(in.intValue() + 1));
					return true;
				}
			}
		}
		
		return true;
	}

	Transcript readFeature(GFFObject obj, Transcript trpt) {

		trpt.addAttribute(GFFObject.GENE_ID_TAG, obj.getAttribute(GFFObject.GENE_ID_TAG));
		
		if (obj.getFeature().equalsIgnoreCase("exon")) {
			String exonID = obj.getAttribute(GFFObject.EXON_ID_TAG);
			Exon e = new Exon(trpt, exonID, obj.getStart(), obj.getEnd());
			e.setChromosome(obj.getSeqname());
			e.setStrand((byte) obj.getStrand());
			++readExons;

		} else if (obj.getFeature().equalsIgnoreCase("CDS")
				|| obj.getFeature().equals("start_codon")
				|| obj.getFeature().equals("stop_codon")) {
			trpt.addCDS(obj.getStart(), obj.getEnd());
			Translation trans= trpt.getTranslations()[0];
			Vector<Transcript> vtx= new Vector<Transcript>(1); // TODO very inefficient 
			if (obj.getFeature().equals("start_codon")) {
				SpliceSite ss= new SpliceSite(trans.get5PrimeEdge(), SpliceSite.TYPE_CODON_START, trpt.getGene());
				trpt.getGene().addSpliceSite(ss, vtx);
				trans.setCodonStart(ss);
			} else if (obj.getFeature().equals("stop_codon")) {
				SpliceSite ss= new SpliceSite(trans.get3PrimeEdge(), SpliceSite.TYPE_CODON_STOP, trpt.getGene());
				trpt.getGene().addSpliceSite(ss, vtx);
				trans.setCodonStop(ss);
			}
			Object[] keys = obj.getAttributes().keySet().toArray();
			for (int i = 0; i < keys.length; i++) {
				if ((!GFFObject.GENE_ID_TAG.equals(keys[i]))
						&& (!GFFObject.TRANSCRIPT_ID_TAG.equals(keys[i]))
						&& (!GFFObject.EXON_ID_TAG.equals(keys[i]))) { // TODO
																		// make
																		// it
																		// nice:
																		// protein_id
																		// "...",
																		// other_id
																		// "..."
					String putProtID = obj.getAttribute((String) keys[i]);
					String[] idTokens = putProtID.split(" ");
					for (int j = 0; j < idTokens.length; j++)
						trpt.getTranslations()[0].addProteinID(idTokens[j]);
				}
			}
		} else if(obj.getFeature().equals("start_codon")|| obj.getFeature().equals("stop_codon")) {
			
		} else {
			int y;
			for (y = 0; readFeatures != null && y < readFeatures.length; y++)
				if (readFeatures[y].equalsIgnoreCase(obj.getFeature()))
					break;
			if (readFeatures == null || y < readFeatures.length) {
				String regID = obj.getFeature();
				DirectedRegion newReg = new DirectedRegion(obj.getStart(), obj
						.getEnd(), trpt.getStrand());
				newReg.setID(obj.getAttribute(obj.getFeature() + "_id"));
				newReg.setScore(obj.getScore());
				newReg.setChromosome(trpt.getChromosome());
				DirectedRegion[] reg = null;
				if (trpt.getAttributes() != null)
					reg = (DirectedRegion[]) trpt.getAttributes().remove(regID);
				if (reg == null)
					reg = new DirectedRegion[] { newReg };
				else {
					int i;
					for (i = 0; i < reg.length; i++) {
						if (reg[i].getID().equals(newReg.getID())) {
							if (reg[i].getScore() != newReg.getScore())
								System.out
										.println("WARNING: regions with non-matching scores "
												+ newReg + ", " + reg[i]);
							reg[i].set5PrimeEdge(Math.min(newReg
									.get5PrimeEdge(), reg[i].get5PrimeEdge()));
							reg[i].set3PrimeEdge(Math.max(newReg
									.get3PrimeEdge(), reg[i].get3PrimeEdge()));
							break;
						}
					}
					if (i == reg.length)
						reg = (DirectedRegion[]) ArrayUtils.add(reg, newReg);
				}

				trpt.addAttribute(regID, reg); // "domain",
			}

		}

		return trpt;
	}

	/**
	 * Requires the gtf file to be ordered by chromosomes AND that data (CDS,
	 * exon) from the same transcript is in adjacent lines (which is good for
	 * RefSeq in worm, where the same tID/gID exists on 4 diferent chromosomes).
	 * It does NOT require any order of the transcript starts within the
	 * chromosome, clustering is performed exhaustively.
	 * 
	 * 
	 * retrieves every transcript that is overlapping a certain region of
	 * <code>ovlReg</code>
	 * 
	 * @return
	 */
	/*
	 * Note that the clustering cannot be separated from the dynamic gene
	 * building, since reduncancy checks for splice sites and exons within the
	 * genes depend on it.
	 */
	BufferedReader buffy= null;
	private BufferedReader getBuffy() {

		if (buffy== null) {
			if (getInputFile() != null) {
				try {
					if (size== bytesRead) {
						genes= null;
						gtfObj= null;
						return null;
					}
					inputStream = new FileInputStream(getInputFile());
					inputStream.skip(bytesRead);
					buffy = new BufferedReader(new InputStreamReader(inputStream));
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		
		return buffy;
	}
	
	boolean readAll= false;
	boolean warnFirstSkip= true;
	public void read() {

		BufferedReader buffy = getBuffy();
		if (buffy== null)
			return;
		
		if (bytesRead== 0) {
			skippedObjects= 0;
			warnFirstSkip= true;
		}
		
		clustered = false;
		
		Vector gtfV = null;
		if (readGTF)
			gtfV = new Vector();
		Vector<Gene> geneV= null;
		int cnt= 0, cntTrpt= 0;
		if (readGene)
			geneV= new Vector<Gene>();
		String line = null;
		try {
			Transcript lastTrpt = null;
			String lastChrID = null;
			byte lastStrand= 0;
			Transcript trpt = null;

			if (keepOriginalLines) {
				if (vLines== null)
					vLines= new Vector<String>();
				else
					vLines.removeAllElements();
			} else
				vLines= null;
			boolean inited= false;
			if (reuse&& lastLine!= null) {
				if (keepOriginalLines)
					vLines.add(lastLine);					
				line= lastLine;
				lastLine= null;
				inited= true;
			}
			
			while (true) {

				// read line
				if (inited) 
					inited= false;
				else {
					line= buffy.readLine();
					if (line == null) {
						if (trpt!= null) {	// necessary
							if (readGene&& (geneV.size()== 0|| geneV.lastElement()!= trpt.getGene()))
								geneV.add(trpt.getGene());
						} 
						
						if (!silent)
							printReadStatistics();
						break;
					} else {
						if (vLines!= null)
							vLines.add(line);
						++nrLinesRead;
						bytesRead += line.length() + getLineSeparator().length();				
                        if(!silent && stars){
                            Log.progress(bytesRead, size);
                        }
					}
				}
				
				// check line
				if (line.startsWith("#"))
					continue;
				GFFObject obj= readBuildObject(line);
				if (!checkObject(obj)) {	// object based criteria
					++skippedObjects;
					if (warnFirstSkip) {
						Log.warn("skipped line "+ obj);
						warnFirstSkip= false;
					}
					if (trpt!= null && geneWise
							&& ((readAheadLimit> 0&& cnt== (readAheadLimit+1))|| (readAheadTranscripts> 0&& cntTrpt>= readAheadTranscripts))) {
							// NO: wait for end of LOCUS !!! -> no overlap
							//&& (!obj.getAttribute(GFFObject.TRANSCRIPT_ID_TAG).equals(trpt.getTranscriptID()))) 

						Gene currGene= geneV.elementAt(geneV.size()- 1);
						boolean breakit= true;
						for (int i = 0; i < currGene.getTranscripts().length; i++) {
							Transcript tx= currGene.getTranscripts()[i];
							int tstart= Math.abs(tx.getStart()), tend= Math.abs(tx.getEnd());
							if (tstart> tend) {int h= tstart; tstart= tend; tend= h;}
							if (((obj.getStart()>= tstart&& obj.getStart()<= tend)|| (obj.getEnd()>= tstart&& obj.getEnd()<= tend))
									|| (obj.getAttribute(GFFObject.TRANSCRIPT_ID_TAG).equals(trpt.getTranscriptID()))) {
								breakit= false;
								break;
							}
						}
						if (breakit) { 
							if (reuse)
								lastLine= line;
							if (vLines!= null)
								vLines.remove(vLines.size()- 1);
							break;
						} else
							continue;
					} else {
						continue;
					}
				}
				String chrID = obj.getSeqname();
				if (!chrID.equals(lastChrID)) {		// chromosome
					//++readChrs;
					if (lastChrID != null) { 						
						if (!checkChromosome(chrID)) {
							++skippedObjects;
							if (warnFirstSkip) {
								Log.warn("skipped chromosome "+ chrID);
								warnFirstSkip= false;
							}
							getReadChr().add(lastChrID);
							ArrayUtils.addUnique(getSkippedChr(), chrID);
							buffy= skipToNextChromosome(buffy, size, chrID);
							if (buffy== null)
								break;
						} else {	// only if not sweeped
							bytesRead-= line.length()+ lineSeparator.length();
							--nrLinesRead;
						}
						if ((chromosomeWise|| geneWise)&& (!readAll)) {
							if (vLines!= null&& vLines.size()> 0)
								vLines.remove(vLines.size()-1);
							if (reuse)
								lastLine= line;
							if (readGene&& (geneV.size()== 0|| geneV.lastElement()!= trpt.getGene())) 
								add(geneV, trpt);
							break;
						}
					}
					lastChrID= chrID;
				}
				
				if (lastStrand!= obj.getStrand()) {	// strand
					if (lastStrand== 0)
						lastStrand= (byte) obj.getStrand();
					else if (chromosomeWise&& (!readAll)) {
						if (readGene&& (geneV.size()== 0|| geneV.lastElement()!= trpt.getGene())) 
							add(geneV, trpt);
						if (reuse) 
							lastLine= line;						
						if (vLines!= null&& vLines.size()> 0)
							vLines.remove(vLines.size()-1);
						--nrLinesRead;
						bytesRead-= line.length()+ lineSeparator.length();
						break;
					}
				}

				// gtf objects are finished here
				if (readGTF) {
					gtfV.add(obj);
					if (gtfV.size()== getLimitGTFObs()) {
						if (reuse)
							lastLine= line;
						break;
					}
				}
				if (!readGene)
					continue;

				// build gene
				String tid = obj.getAttribute(GFFObject.TRANSCRIPT_ID_TAG);
                if (tid== null) {
                    throw new RuntimeException(
                            "I have no transcript ID, and I want to scream!\n" +
                            line+ "\n"+
                            obj.getAttribute(GFFObject.TRANSCRIPT_ID_TAG)
                    );
                }
				if (lastTrpt== null|| (!tid.equals(lastTrpt.getTranscriptID()))) {
					boolean overlap = true;
					if (!checkOverlapRegion(trpt)) // should be ok for the
													// moment. maybe also check
													// in merging step
						overlap = false;

					// save gene
					if (overlap) {
						if (trpt != null) {
							if (readGene&& (!geneWise)&& (geneV.size()== 0|| geneV.lastElement()!= trpt.getGene())) 
								add(geneV, trpt);
						}
					} else {
						skippedTranscripts.add(lastTrpt.getTranscriptID());
						Log.warn("skipped transcript "+ lastTrpt.getTranscriptID());
					}

					Gene ge = new Gene(Gene.getUniqueID());
					ge.setStrand(obj.getStrand());
					ge.setChromosome(obj.getSeqname());
					ge.setSpecies(species);
					trpt = new Transcript(tid);
					lastTrpt = trpt;
					trpt.setStrand(obj.getStrand());
					trpt.setChromosome(obj.getChromosome());
					trpt.setSource(obj.getSource());
					ge.addTranscript(trpt);
					++readTranscripts;
					++cntTrpt;
				}

					// read the feature
				trpt = readFeature(obj, trpt);
				// pull identical exon attributes up to transcript level
//				if (trpt.getExons() != null && trpt.getExons()[0].getAttributes() != null) {
//					Object[] keys = trpt.getExons()[0].getAttributes().keySet().toArray();
//					for (int i = 0; i < keys.length; i++) {
//						int j;
//						for (j = 1; j < trpt.getExons().length; j++) {
//							if (!trpt.getExons()[j].getAttribute(keys[i])
//									.equals(
//											trpt.getExons()[0]
//													.getAttribute(keys[i])))
//								break;
//						}
//						if (j == trpt.getExons().length) {
//							Object val = trpt.getExons()[0]
//									.getAttribute(keys[i]);
				Features:	
//							for (int k = 0; k < trpt.getExons().length; k++)
//								trpt.getExons()[k].getAttributes().remove(
//										keys[i]);
//							trpt.addAttribute(keys[i], val);
//						}
//					}
//				}
				
					// gene clustering
				if (geneWise) {					
						
					boolean newGene= false;
					if (geneV.size()== 0) {
						newGene= true;
					} else {
						if (trpt.getGene()!= geneV.elementAt(geneV.size()-1)) {
							if (clusterGenes) {
								if (trpt.overlaps((DirectedRegion) geneV.elementAt(geneV.size()- 1)))
									geneV.elementAt(geneV.size()-1).merge(trpt.getGene());
								else
									newGene= true;
							} else {
//                                Object a1 = trpt.getGene().getAttribute(GFFObject.GENE_ID_TAG);
//                                Object a2 = geneV.elementAt(geneV.size() - 1).getAttribute(GFFObject.GENE_ID_TAG);
                                Object a1 = trpt.getGene().getGeneID();
                                Object a2 = geneV.elementAt(geneV.size() - 1).getGeneID();
                                if (a1.equals(a2))
									geneV.elementAt(geneV.size()-1).merge(trpt.getGene());
								else
									newGene= true;
							}
						}
					}
					if (newGene) {
						++cnt;
						if ((readAheadLimit> 0&& cnt== (readAheadLimit+1))|| (readAheadTranscripts> 0&& cntTrpt> readAheadTranscripts)) {
							if (reuse) {
								lastLine= line;
							}
							--nrLinesRead;
							if (vLines!= null&& vLines.size()> 0)
								vLines.remove(vLines.size()-1);
							bytesRead -= line.length() + lineSeparator.length();
							break;
						} else {
							++readGenes;
							add(geneV, trpt);
						}
					}
						
				}
				
			} // end while(true)
			
			if (getInputFile()!= null&& !reuse) {
				if (buffy!= null)
					buffy.close();
				this.buffy= null;
			}
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
			//
		if (readGene&& geneV!= null) {	// happens on "empty lines" (pb with terminator, single line files)
			this.genes = (Gene[]) ArrayUtils.toField(geneV); // doesnt matter which sorting, no?!
			if (this.genes!= null)
				nrGenes+= this.genes.length;
//			BufferedWriter www= new BufferedWriter(new FileWriter("N:\\txmap.reader", true));
			for (int i = 0; this.genes!= null&& i < this.genes.length; i++) {
				Gene g= this.genes[i];
				nrTranscripts+= g.getTranscriptCount();
				for (int j = 0; j < g.getTranscripts().length; j++) {
					nrExons+= g.getTranscripts()[j].getExons().length;
//					www.write(g.getTranscripts()[j].getTranscriptID()
//							+"\t"+g.getTranscripts()[j].getExons().length+ "\n");					
					if (g.getTranscripts()[j].getExons().length== 0)
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
							System.err.println("[SOS] transcript "
									+g.getTranscripts()[j].getTranscriptID()
									+" has no exons!");
				}
			}
//			www.flush();
//			www.close();
			if (chromosomeWise)
				clustered = false;
		}
		if (readGTF)
			this.gtfObj = (GFFObject[]) ArrayUtils.toField(gtfV);

//		if (bytesRead == size && Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
//			System.err.println();

	}
	
	private void add(Vector<Gene> geneV, Transcript trpt) {
		geneV.add(trpt.getGene());
	}

	/**
	 * Builds up a GTF object from a GTF line.
	 * @param line a GTF line of the format <code>\<seqname\> \<source\> \<feature\>  \<start\> 
	 * \<end\> \<score\> \<strand\> \<frame\> [attributes] [comments]</code>
	 * @return an instance of <code>GTFObject</code> representing that line 
	 */
	protected GFFObject readBuildObject(String line) {

		GFFObject obj = new GFFObject(); 
		
		String[] tokens = line.split("\\s"); 	// actually it should be a tab (specification), but for awk users..
		if (tokens.length < 8) {
			if (!silent)
				System.err.println("line " + nrLinesRead
						+ ": skipped (<8 token)!\n\t" + line);
			return obj;
		}
		
			// init attributes
		obj.setSeqname(GFFObject.filterRedundantStrings(tokens[0]));
		obj.setSource(GFFObject.filterRedundantStrings(tokens[1]));
		obj.setFeature(GFFObject.filterRedundantStrings(tokens[2]));

		try {
			obj.setStart(Integer.parseInt(tokens[3]));
		} catch (NumberFormatException e) {
			if (!silent)
				System.err.println("Error in line "+nrLinesRead+": invalid start position \'"+tokens[3]+"\'.");
		}
		try {
			obj.setEnd(Integer.parseInt(tokens[4]));
		} catch (NumberFormatException e) {
			if (!silent)
				System.err.println("Error in line "+nrLinesRead+": invalid end position \'"+tokens[4]+"\'.");
		}

		try {
			obj.setScore(tokens[5]);
		} catch (NumberFormatException e) {
			//if (!silent)
				System.err.println("Error in line "+nrLinesRead+": invalid score \'"+tokens[5]+"\'.");
		}
		obj.setStrand(GFFObject.parseStrand(tokens[6].trim()));
		
		if (tokens[7].charAt(0)!= GFFObject.SYMBOL_NOT_INITED)	
			try {
				obj.setFrame(Byte.parseByte(tokens[7]));
			} catch (NumberFormatException e) {
				if (!silent)
					System.err.println("Error line "+nrLinesRead+": invalid frame assignment \'"+tokens[7]+"\'.");
			}
			
		String[] attrTokens = new String[0];
		if (tokens.length > 8) {
			String h = line.substring(line.indexOf(tokens[8]),
					line.length()).trim(); // attributes, comments
			attrTokens = h.split(";\\s+");
			for (int i = 0; i < attrTokens.length; i++) {
				h = attrTokens[i].trim();
				h = h.replaceAll("\\s+", " ");
				int sep = Math.max(h.indexOf(' '), h.indexOf('=')); // = in ASD gff
				if (sep < 0) // comments
					break;
				if (sep >= 0) {
					String id = h.substring(0, sep);
					String val = h.substring(sep + 1, h.length());
					int ofStart= 0, ofEnd= 0;
					if (val.charAt(0) == '\"')
						ofStart= 1;
					if (val.charAt(val.length()- 1) == '\"')
						ofEnd= 1;
					if (val.length()>= 2&& val.charAt(val.length()- 2) == '\"')
						ofEnd= 2;
					if (val.length()- ofEnd<= ofStart)
						val= "";	// empty fields
					else
						val = val.substring(ofStart, val.length()- ofEnd);
					obj.addAttribute(id, val);
				}
			}
		}
		
		return obj;

	}
	

	private Gene[] clusterLoci(Gene[] g) {
		if (g == null)
			return null;

		Comparator compi = new DirectedRegion.StartComparator();
		java.util.Arrays.sort(g, compi);

		Vector v = new Vector();
		for (int i = 0; i < g.length; i++) {
			v.add(g[i]);
			int j;
			for (j = i + 1; j < g.length; j++) {
				if (g[i].overlaps(g[j]))
					g[i].merge(g[j]);
				else
					break;
			}
			i = j - 1;
		}

		return (Gene[]) ArrayUtils.toField(v);
	}
	
	public int countTranscripts() {
		HashMap<String, String> hash= null;

		try {
			BufferedReader buffy = new BufferedReader(getReader());
			String s= null; StringTokenizer st;
			while (buffy.ready()) {
				s= buffy.readLine();
				String[] a= s.split("\\s");
				if (a.length< 8)
					continue;
				if (hash== null) 
					hash= new HashMap<String, String>((int) (getInputSize()/ (20f *s.length())));	// 10 exons, 2x for "CDS" and "exon"
				for (int i = 8; i < a.length; i++) {
					if (a[i].equals(GFFObject.TRANSCRIPT_ID_TAG)) {
						hash.put(a[0]+a[i+1], a[i+1]);
						break;
					}
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return hash.size();
	}

	/**
	 * per definition there cannot 2 genes with same start in the arrays
	 * (merged!) therefore binary search works
	 * 
	 * not so sure, during the multiple merges there can happen the case that
	 * two have the same
	 * 
	 * @param geneReg
	 * @param gene
	 * @param startCompi
	 * @param endCompi
	 * @return
	 */
	private boolean mergeGene(Gene[][] geneReg, Gene gene,
			Comparator startCompi, Comparator endCompi) {
		boolean merged = false;
		if (geneReg[0] == null) {
			geneReg[0] = new Gene[] { gene };
			geneReg[1] = geneReg[0]; // not needed here, see efficient
										// version..
		} else {
			Vector vRest = new Vector(geneReg[0].length);
			for (int i = 0; i < geneReg[0].length; i++) {
				if (geneReg[0][i].overlaps(gene))
					gene.merge(geneReg[0][i]);
				else
					vRest.add(geneReg[0][i]);
			}
			vRest.add(gene);
			geneReg[0] = (Gene[]) ArrayUtils.toField(vRest);
			geneReg[1] = geneReg[0]; // not needed here, see efficient
										// version..
		}

		return merged;
	}

	/**
	 * Sorts data read from <code>inputFile</code> to a default
	 * target file.
	 * @return a file handle to which data has been written
	 */
	public File sort() {
		File sortedFile= new File(FileHelper.append(inputFile.getAbsolutePath(), "_sorted", false, null));
        // this alwasy creates unzipped files, so
        // we have to remove any .gz extensions!
        if(sortedFile.getName().toLowerCase().endsWith(".gz")){
            String fn = sortedFile.getAbsolutePath();
            sortedFile = new File(fn.substring(0, fn.length()-3));
        }
		sort(sortedFile);
        return sortedFile;
	}

	/**
	 * Sorts the <code>inputFile</code> to a certain target
	 * file.
	 * @param targetFile the file to which sorted data is 
	 * written
	 */
	public void sort(File targetFile) {
        if(!silent && stars){
            Log.progressStart("sorting GTF file");
        }
        try {
			OutputStream oStream= new FileOutputStream(targetFile);
			sort(oStream);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		} finally {
	        if(!silent && stars){
	            Log.progressFinish(StringUtils.OK, true);
	        }
		}
	}

	/**
	 * Sorts the source data read from <code>inputFile</code>
	 * and writes the sorted result to the given 
	 * <code>OutputStream</code>.
	 * @param outputStream stream to which sorted data is 
	 * written
	 */
	public void sort(OutputStream outputStream) {
		// attentionAttention, the nrs have to be sorted - look in find()
        int[] fieldNrs= new int[]{0,3,6,-1};	// ,-1 for transcriptid

        Map<byte[], Integer> tssMap = null;
        try {
            tssMap = getTSS(fieldNrs);  
		} catch (Exception error) {
			throw new RuntimeException(error);
		} finally {
			System.gc();
		}

        try {
			sort(outputStream, fieldNrs, tssMap);
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			System.gc();
		}

        
	}
	
    /**
     * Sort data read from an <code>InputStream</code> and 
     * write it to an <code>OutputStream</code>.
     *
     * @param istream the source stream
     * @param ostream the target stream
     * @throws Exception in case of any errors
     */
	private void sort(OutputStream outStr, int[] fieldNrs, Map<byte[], Integer> transcriptPositions) throws Exception {
        IOHandler io = IOHandlerFactory.createDefaultHandler();
        PipedOutputStream out = null;
        PipedInputStream in = null;
        BufferedWriter writer = null;
        Future sorterThread = null;

        ByteArrayCharSequence cs = null;
        long lineCount = 0l;
        ByteArrayCharSequence[] fields = null;
		try {
			cs= new ByteArrayCharSequence(1000);

            
            out = new PipedOutputStream();
            in = new PipedInputStream(out);
            writer= new BufferedWriter(new OutputStreamWriter(out));

            sorterThread = Sorter.create(in, outStr, true, "\\s")
                    .field(0, false) // chr
                    .field(6, false) // strand
                    .field(new GlobalTranscriptPositionComparator(transcriptPositions, new int[]{fieldNrs[0], fieldNrs[3]}))
                            //.field(1, true) // tpos
                    .field(fieldNrs[3], false) // tid
                    .field(3, true).sortInBackground();// pos
            fieldNrs= new int[]{fieldNrs[0], fieldNrs[3]};	// 090901: start pos from map
			HashSet<String> setInvalidTx= new HashSet<String>();
			int nrInvalidLines= 0;

			InputStream fileInput = getInputStream();
            io.addStream(fileInput);

            while (io.readLine(fileInput,cs) != -1){
                lineCount++;
                // issue #56 make sure we skip empty lines and comment lines
                if (cs.length() == 0 || cs.charAt(0)== '#'){
                    continue;
                }

                fields= ByteArrayCharSequence.split(cs, SPLITTER_PATTERN,fieldNrs);	// 1,4,-1
                byte[] key= ByteArrayCharSequence.merge(fields[1], fields[0]);
                int val= transcriptPositions.get(key);
                if (val== Integer.MIN_VALUE) {
                    ++nrInvalidLines;
                    setInvalidTx.add(fields[1].toString());
                }else{
				    writer.write(cs.toString());
				    writer.write(lineSeparator);
                }
			}
            writer.flush();
            writer.close();
            sorterThread.get();

			transcriptPositions = null;

			System.gc();

			
			if (nrInvalidLines> 0){
                Log.message("\tskipped "+nrInvalidLines+" lines in " + setInvalidTx.size()+ " transcripts");
                Log.message("");
                for (String s : setInvalidTx) {
                    Log.message(s + " ");
                }
                Log.message("");
            }
        }catch (Exception error){
            logError(fieldNrs, cs, lineCount, fields, error);
            throw error; // rethrow
		} finally {
            io.close();
            if(out != null) try {out.close();} catch (IOException e) {}
            if(in != null) try {in.close();} catch (IOException e) {}
            if(writer != null) try {writer.close();} catch (IOException e) {}
            if(sorterThread != null){
                sorterThread.cancel(true);
            }
		}
	}

    /**
     * Helper method to print more usefull infomation to log error if the gtf can not be parsed.
     *
     * @param fieldNrs the field numbers inspected by the parser
     * @param cs the current line
     * @param lineCount the line count
     * @param fields the parsed fields
     * @param error the exception
     */
    private void logError(int[] fieldNrs, ByteArrayCharSequence cs, long lineCount, ByteArrayCharSequence[] fields, Exception error) {
        Log.error("Error while sorting GTF File : " + error.getMessage());
        if(cs != null){
            Log.error("Line "+lineCount+" caused the problem : " + cs.toString());
            Log.error("Please check your GTF file and watch out for too much tabs or spaces!");
            Log.error("We currently use '"+SPLITTER_PATTERN+"'. This single tab or space characters.");
            Log.error("Note that the GTF definition states that the first eight fields are separated by tab (\\t).");
            Log.error("Our parser also allows a single space, but not multiple spaces !");
            Log.error("We look for fields : "+ Arrays.toString(fieldNrs) + ".");
            if(fields  != null){
                Log.error("The field content for the problematic line : ");
                for (int i = 0; i < fields.length; i++) {
                    switch (i){
                        case 0: Log.error("Field " + fieldNrs[i] + " : "+fields[i] + " - The chromosome, parsed as Text");break;
                        case 1: Log.error("Field " + fieldNrs[i] + " : "+fields[i] + " - Start Position, parsed as Number");break;
                        case 2: Log.error("Field " + fieldNrs[i] + " : "+fields[i] + " - Strand, parsed as Text");break;
                        case 3: Log.error("Field " + fieldNrs[i] + " : "+fields[i] + " - Transcript ID, parsed as Text");break;
                    }

                }
            }
        }
    }

    /**
     * Find the field index of the transcript id or return -1
     *
     * @param input the line
     * @return transcriptIDindex the field index if the transcript ID
     */
	private int findTranscriptID(CharSequence input) {
    	int index = 0;
		Matcher m= SPLITTER_PATTERN.matcher(input);
	    for(int mCtr= 0;m.find();index= m.end(), mCtr++) {
	    	if (mCtr< 8)
	    		continue;	    	
	    	CharSequence match = input.subSequence(index, m.start());
	    	Matcher m2= TRANSCRIPTID_PATTERN.matcher(match);
	    	if (m2.find())
	    		return mCtr+1;
	    }
	    return -1;
	}
	
    /**
     * Create a map using the concatenated transcript ID and the chromosome name as key
     * and the minimal transcript position as value. This is later used to sort the
     * entries relative to their global position. We pass the field[] here, where the last position
     * is filled with the identified index of the transcript IDs. An exception is thrown if
     * a problem occurs while parsing the file.
     *
     * @param fieldNrs the fields to check. Position 3 is set the to the transcriptID index
     * @return transcriptMap map from the concatenated transcriptID+chromosome as key and the minimal global position as value
     * @throws Exception in case of any error
     */
    protected Map<byte[], Integer> getTSS(int[] fieldNrs) throws Exception{
    	
        IOHandler io = IOHandlerFactory.createDefaultHandler();
        ByteArrayCharSequence cs= new ByteArrayCharSequence(1000);
        ByteArrayCharSequence[] fields = null;
        int lineCounter = 0;
        try{
            /*
             estimate line count to estimate map size
             */
			long estLineCount = (getInputSize() / 100);
			int estIDCount = 0;
			if (estLineCount <= 50000){
                /*
                reference annotation
                10 exons per transcript, 2 for CDS and exon line
                 */
				estIDCount = (int) (estLineCount / 40);
            }else if (estLineCount <= 500000){
                /*
                mRNA collection
                7 exons per transcript, no CDS
                 */
				estIDCount = (int) (estLineCount / 6); //
            }else{
				// ESTs, 25 mio lines -> 4 mio IDs
				estIDCount = (int) (estLineCount / 6);
            }


            // the map that maps from chr+transcriptID -> min transcript position
			MyArrayHashMap<byte[], Integer> transcriptPositions = new MyArrayHashMap<byte[], Integer>(estIDCount);
			transcriptPositions.setIncrementSize((int) (estIDCount * 0.3));

            InputStream reader = getInputStream();
            io.addStream(reader);


            //int[] fieldNrs= new int[]{1,4,7,-1};	// ,-1 for gene
            int bytesRead = 0;

            while(io.readLine(reader, cs) != -1){
            	
                lineCounter++;
				bytesRead += cs.length() + getLineSeparator().length();

				if (cs.charAt(0)== '#') {
					continue;
				}
                // also ignore empty lines
                if(cs.length() == 0){
                    continue;
                }

				if (fieldNrs[3]< 0) {
					fieldNrs[3]= findTranscriptID(cs);
					if (fieldNrs[3]<0) {
                        throw new RuntimeException("No transcript ID found.\n"+cs);
					}
				}

				fields= ByteArrayCharSequence.split(cs, SPLITTER_PATTERN, fieldNrs);	// 1,4,tid,gid
				for (int i = 0; i < fields.length; i++) {
					if (fields[i]== null) {
                        throw new RuntimeException("I could not find field number "+fieldNrs[i]+" in line " + lineCounter
                        +"\n"
                        +"\tline skipped check format of GTF file"+"\n"
                        +"\t(first 8 fields and transcript_id in the same column!)");
					}
				}
					// transcript start
				int start = fields[1].parseInt();
				byte[] key= ByteArrayCharSequence.merge(fields[3], fields[0]);
				byte strand= GFFObject.parseStrand(fields[2]);
				strand=(strand== 0)?1:strand;	// convert unknown to 1, for multiplication afterwards map.put()

				Integer oldVal = transcriptPositions.get(key);
//				String DEBUG= key.toString();
				// 090901: consistency check of aligned transcripts,
				// must be to the same strand on the same chromosome
				if (oldVal!= null&& !oldVal.equals(Integer.MIN_VALUE)) {
					if (strand* oldVal< 0) {
						oldVal= Integer.MIN_VALUE;
						transcriptPositions.put(key, oldVal);
						continue;
					}
				}

				if (oldVal == null || Math.abs(oldVal) > start) {
					transcriptPositions.put(key, (strand* start));
				}
				
			}
            return transcriptPositions;
        } catch (Exception error) {
            logError(fieldNrs, cs, lineCounter, fields, error);
        	throw error; // rethrow
        }finally {
            io.close();
        }
    }

	
	private int fieldTID= -1, fieldGID= -1;
	public int getFieldTID() {
		if (fieldTID< 0) {
			guessFieldTIDGID();
		}
		return fieldTID;
	}
	public int getFieldGID() {
		if (fieldGID< 0) {
			guessFieldTIDGID();
		}
		return fieldGID;
	}
	private void guessFieldTIDGID() {
		if (inputFile== null)
			return;
		try {
			BufferedReader buffy= new BufferedReader(getReader());
			String s= buffy.readLine();
			String[] ss= s.split("\\s");
			for (int i = 9; i < ss.length; i+=2) {
				if (ss[i].equals(GFFObject.TRANSCRIPT_ID_TAG))
					fieldTID= i+ 1;
				else if (ss[i].equals(GFFObject.GENE_ID_TAG))
					fieldGID= i+ 1;
				if (fieldTID>= 0&& fieldGID>= 0)
					break;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

    /**
     * Checks that the current file is a valid GTF file and checks
     * that the file is sorted!
     *
     * @return sorted true if file is valid and sorted
     * @throws RuntimeException in case the file is not valid
     */
	public boolean isApplicable() {
		reset();
		long lines= isApplicable(getInputFile(), clusterGenes);
		if (lines< 0)
			return false;
		return true;
	}
	
	/**
	 * Checks for correct sorting, returns number of lines read (<0 if not applicable).
	 * @param inputFile the file from which is read
	 * @param clusterGenes indicates whether native gene clustering is applied
	 * @return number of lines read, or -(number of lines read) up to the unsorted
	 * entry
	 */
	protected long isApplicable(File inputFile, boolean clusterGenes) {
		try {
			InputStream fis= getInputStream();
			long linesOK= isApplicable(fis, FileHelper.getSize(inputFile));
			fis.close();
			return linesOK;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	/**
	 * Checks for correct sorting, returns number of lines read (<0 if not applicable).
	 * <b>Note:</b> does not close the given stream.
	 * @param inputStream stream from which is read
	 * @param clusterGenes indicates whether native gene clustering is applied
	 * @param size total size of data in the stream, if known, otherwise <= 0
	 * @return number of lines read, or -(number of lines read) up to the unsorted
	 * entry
	 */
	protected long isApplicable(InputStream inputStream, long size) {
		
		long bytesRead = 0l;
		BufferedReader buffy = new BufferedReader(new InputStreamReader(inputStream));

		int lineCtr = 0;
		try {
			String lastChrID = null, lastGID = null, lastTID = null, lastStrand = null;
			HashMap<String, String> chrMap = new HashMap<String, String>(25, 1f), 
				tidMap = new HashMap<String, String>(), gidMap = new HashMap<String, String>();
			int tidField = -1, gidField = -1, lastStart = -1;
			String line;
            if(size > 0){
                Log.progressStart("Checking GTF");
            }
			while ((line = buffy.readLine())!= null) {	// TODO && !isStop()
	
				if (line.startsWith("#"))
					continue;
				if (line == null)
					break;
	
				++lineCtr;
				bytesRead += line.length() + 1;
                if(size> 0){	// TODO !silent && stars
                    Log.progress(bytesRead, size);
                }
				String[] tokens = line.split("\\s");
														
					// trim attribute strings
				for (int i = 0; i < tokens.length; i++) {
					if (tokens[i].endsWith(";")) {
						if (tokens[i].startsWith("\""))
							tokens[i] = tokens[i].substring(1, tokens[i]
									.length() - 2);
						else
							tokens[i] = tokens[i].substring(0, tokens[i]
									.length() - 1);
					}
				}
				
				if (tokens.length < 8) {
				    throw new RuntimeException("GTF Error line " + lineCtr+ " - less than 8 tokens.");
				}
	
				// changes chr or strand
				if ((!tokens[0].equals(lastChrID))
						|| (!tokens[6].equals(lastStrand))) {
					String chrStrand = tokens[0] + tokens[6];
					if (chrMap.get(chrStrand) != null) {
							Log.warn("Unsorted in line " + lineCtr
									+ " - chr/strand " + tokens[0] + " "
									+ tokens[6] + " already read.");
						return (-lineCtr);
					}
					chrMap.put(chrStrand, chrStrand);
					lastChrID = tokens[0];
					lastStrand = tokens[6];
					tidMap = new HashMap<String, String>(); 
					gidMap = new HashMap<String, String>();
					lastTID= null;	// re-init per chromosome
					lastGID= null;
					lastStart= -1;
				} else {
					tokens[0]= lastChrID;
					tokens[6]= lastStrand;
				}
	
				// changes transcript ID
				if (tidField < 0
						|| (!tokens[tidField - 1]
								.equals(GFFObject.TRANSCRIPT_ID_TAG))) {
						// try to find TID field
					for (int i = 8; i < tokens.length; i++) {
						if (tokens[i].equals(GFFObject.TRANSCRIPT_ID_TAG)) {
							tidField = i + 1;
							break;
						}
					}
					if (tidField < 0) {
                        throw new RuntimeException("Error line " + lineCtr+ " - no TID:\n"+line);
					}
				}
				if (tidField>= 0&& !tokens[tidField].equals(lastTID)) {
					if (tidMap.get(tokens[tidField]) != null) {
//                            if(!silent && stars){
//                                Log.progressFailed("Unsorted in line " + lineCtr
//									+ " transcript id " + tokens[tidField]
//									+ " used twice, on: " + tokens[0] + ","
//									+ tidMap.get(tokens[tidField]));
//                            }else{
                                Log.warn("Unsorted in line " + lineCtr
                                        + " transcript id " + tokens[tidField]
                                        + " used twice, on: " + tokens[0] + ","
                                        + tidMap.get(tokens[tidField]));
//                            }
                        return (-lineCtr);
					}
					tidMap.put(tokens[tidField], tokens[0]);
					if (clusterGenes) {
						int newStart = Integer.parseInt(tokens[3]);
						if (lastTID != null) {
							if (lastStart > newStart) {
//	                            if(!silent && stars){
//	                                Log.progressFailed("Unsorted in line " + lineCtr
//											+ " - cannot perform gene clustering: "
//											+ tokens[0] + " " + tokens[6] + " "
//											+ tokens[tidField] + " @ " + tokens[3]
//											+ " after " + lastTID + " @ " + lastStart);
//	                            }else{
									Log.warn("Unsorted in line " + lineCtr
											+ " - cannot perform gene clustering: "
											+ tokens[0] + " " + tokens[6] + " "
											+ tokens[tidField] + " @ " + tokens[3]
											+ " after " + lastTID + " @ " + lastStart);
//	                            }
								buffy.close();
								return (-lineCtr);
							}
						}
						lastStart = newStart;
					}
					lastTID = tokens[tidField];
				} else {
					tokens[tidField]= lastTID;
				}
	
				// changes in gene ID
				if (!clusterGenes) {
					if (gidField < 0
							|| (!tokens[gidField - 1]
									.equals(GFFObject.GENE_ID_TAG))) {
						for (int i = 8; i < tokens.length; i++) {
							if (tokens[i].equals(GFFObject.GENE_ID_TAG)) {
								gidField = i + 1;
								break;
							}
						}
						if (tidField < 0) {
                            if(size> 0){	// TODO !silent && stars
                                Log.progressFailed("Error line " + lineCtr+ " - no GID.");
                            }else{
								Log.warn("Error line " + lineCtr+ " - no GID.");
                            }
							return (-lineCtr);
						}
					}
					if (!tokens[gidField].equals(lastGID)) {
						if (gidMap.get(tokens[gidField]) != null)	// TODO && !silent
							Log.warn("Warning line " + lineCtr+ " gene id " + tokens[gidField]
									+ " used twice, on: " + tokens[0] + ","
									+ gidMap.get(tokens[gidField]));
						gidMap.put(tokens[gidField], tokens[0]);
						lastGID = tokens[gidField];
					}
				}
			}
	
			buffy.close();

            if(size> 0){	// TODO !silent && stars
                Log.progressFinish(StringUtils.OK, true);
            }

		} catch (Exception e) {
            if(size> 0){	// TODO !silent && stars
                Log.progressFailed("ERROR");
            }
            Log.error("Error while checking GTF file!", e);
		}finally {
            if(buffy != null) try {buffy.close();} catch (IOException e) {}
        }

		return lineCtr;
	}

	public void setFiltTrptIDs(String[] filtTrptIDs) {
		this.filtTrptIDs = filtTrptIDs;
	}

	public void setSilent(boolean silent) {
		this.silent = silent;
	}

	public String[] getFiltSomeIDs() {
		if (filtSomeIDs == null)
			return null;

		String[] result = new String[filtSomeIDs.size()];
		Object[] keys = filtSomeIDs.keySet().toArray();
		for (int i = 0; i < result.length; i++)
			result[i] = (String) keys[i];
		return result;
	}

	public void setFiltSomeIDs(String[] filtSomeIDs) {
		this.filtSomeIDs = new HashMap(filtSomeIDs.length);
		for (int i = 0; i < filtSomeIDs.length; i++)
			this.filtSomeIDs.put(filtSomeIDs[i].toUpperCase(), new Integer(0));

	}

	public void setReadAllLines() {
		filtChrIDs = null;
		filtGeneIDs = null;
		filtRegs = null;
		filtSomeIDs = null;
		filtTrptIDs = null;
		readFeatures = null;
		noIDs = null;
	}

	public String[] getFiltSomeIDsNotFound() {
		if (filtSomeIDs == null)
			return null;

		Vector v = new Vector();
		Object[] keys = filtSomeIDs.keySet().toArray();
		for (int i = 0; i < keys.length; i++) {
			Integer cnt = (Integer) filtSomeIDs.get(keys[i]);
			if (cnt.intValue() == 0)
				v.add(keys[i]);
		}
		return (String[]) ArrayUtils.toField(v);
	}

	public boolean isChromosomeWise() {
		return chromosomeWise;
	}

	public void setChromosomeWise(boolean chromosomeWise) {
		this.chromosomeWise = chromosomeWise;
		geneWise= !chromosomeWise;
		strandWise= !chromosomeWise;
	}

	public String[] getFiltChrIDs() {
		return filtChrIDs;
	}

	public void setFiltChrIDs(String[] filtChrIDs) {
		this.filtChrIDs = filtChrIDs;
	}

	public Gene[] getGenes() {
		if (genes == null)
			return null;
		if ((!geneWise)&& chromosomeWise&& (!clustered)) {
			for (int i = 0; i < genes.length; i++)
				genes[i].setConstruct(false);
			genes = clusterLoci(genes);
			clustered = true;
			readGenes += this.genes.length;
		}

		return genes;
	}

	public void setGenes(Gene[] genes) {
		this.genes = genes;
	}

	public DirectedRegion[] getFiltRegs() {
		return filtRegs;
	}

	public void setFiltRegs(DirectedRegion[] filtRegs) {
		this.filtRegs = filtRegs;
	}

	public boolean isReadGene() {
		return readGene;
	}

	public void setReadGene(boolean readGene) {
		this.readGene = readGene;
	}

	public boolean isReadGTF() {
		return readGTF;
	}

	public void setReadGTF(boolean readGTF) {
		this.readGTF = readGTF;
	}

	public String[] getReadFeatures() {
		return readFeatures;
	}

	public void setReadFeatures(String[] readFeatures) {
		this.readFeatures = readFeatures;
	}

	public void addReadFeature(String newFeature) {
		String[] r = new String[readFeatures.length + 1];
		for (int i = 0; i < readFeatures.length; i++) {
			r[i] = readFeatures[i];
		}
		r[r.length - 1] = newFeature;
		readFeatures = r;
	}

	public int getReadGenes() {
		return readGenes;
	}
	
	public int getReadTranscripts() {
		return readTranscripts;
	}
	
	public Vector getReadChr() {
		if (readChr == null) {
			readChr = new Vector<String>();
		}
	
		return readChr;
	}

	public Vector getSkippedChr() {
		if (skippedChr == null) {
			skippedChr = new Vector<String>();
		}

		return skippedChr;
	}

	public long getBytesRead() {
		return bytesRead;
	}

	public int getUnclusteredGeneNb() {
		if (this.genes == null)
			return 0;
		return this.genes.length;
	}
	
	public static int getGTFLineToken(String line, String token) {
		String[] tokens= line.split("\\s");
		for (int i = 0; i < tokens.length; i++) {
			if (tokens[i].startsWith("\"")) 
				tokens[i]= tokens[i].substring(1);
			if (tokens[i].startsWith(token))
				return i;	
		}
		
		return -1;
	}

	public boolean isClusterGenes() {
		return clusterGenes;
	}

	public void setClusterGenes(boolean clusterGenes) {
		this.clusterGenes = clusterGenes;
	}

	public int getReadAheadLimit() {
		return readAheadLimit;
	}

	public void setReadAheadLimit(int readAheadLimit) {
		this.readAheadLimit = readAheadLimit;
	}

	public int getLimitGTFObs() {
		return limitGTFObs;
	}

	public void setLimitGTFObs(int limitGTFObs) {
		this.limitGTFObs = limitGTFObs;
	}

	public boolean isGeneWise() {
		return geneWise;
	}

	public void setGeneWise(boolean geneWise) {
		this.geneWise = geneWise;
		if (geneWise) {
			clusterGenes= true;
			chromosomeWise= false;
			strandWise= false;
		}
		
	}

	public boolean isStrandWise() {
		return strandWise;
	}

	public void setStrandWise(boolean strandWise) {
		this.strandWise = strandWise;
	}

	public String[] getAllowSources() {
		return allowSources;
	}

	public void setAllowSources(String[] allowSources) {
		this.allowSources = allowSources;
	}

	public int getReadAheadTranscripts() {
		return readAheadTranscripts;
	}

	public void setReadAheadTranscripts(int readAheadTranscripts) {
		this.readAheadTranscripts = readAheadTranscripts;
	}

	public long getInputSize() {
		return size;
	}

	public String[] getNoIDs() {
		return noIDs;
	}

	public void setNoIDs(String[] noIDs) {
		this.noIDs = noIDs;
	}

	public boolean isStars() {
		return stars;
	}

	public void setStars(boolean stars) {
		this.stars = stars;
	}

	public boolean isPrintStatistics() {
		return printStatistics;
	}

	public void setPrintStatistics(boolean printStatistics) {
		this.printStatistics = printStatistics;
	}

	public GFFObject[] getGtfObj() {
		return gtfObj;
	}

	boolean stop= false;
	public boolean isStop() {		
		return stop;
	}

	public boolean setStop() {
		stop= true;
		return true;
	}

	public boolean setStop(boolean stop) {
		if (stop)
			return setStop();
		this.stop= stop;
		return true;
	}

	public void run() {
		; // :)
		
	}
	
	public boolean close() {
		try {
			getBuffy().close();
			return true;
		} catch (Exception e) {
			return false;
		}
	}

	public boolean isReuse() {
		return reuse;
	}

	public void setReuse(boolean reuse) {
		this.reuse = reuse;
	}

	public boolean isReadAll() {
		return readAll;
	}

	public void setReadAll(boolean readAll) {
		this.readAll = readAll;
	}

	public int getNrGenes() {
		return nrGenes;
	}

	public int getNrTranscripts() {
		return nrTranscripts;
	}

	public int getNrExons() {
		return nrExons;
	}

	public int getNrLinesRead() {
		return nrLinesRead;
	}

	@Override
	public int getNrInvalidLines() {
		return skippedObjects;
	}

}
