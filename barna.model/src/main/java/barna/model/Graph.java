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

/*
 * Created on Mar 8, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package barna.model;

//import gphase.NMDSimulator;

import barna.commons.ByteArrayCharSequence;
import barna.commons.utils.ArrayUtils;
import barna.commons.utils.StringUtils;
import barna.model.constants.Constants2;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//import gphase.algo.AlignmentGenerator;
//import gphase.algo.AlignmentWrapper;
//import gphase.db.EnsemblDBAdaptor;
//import gphase.ext.ClustalWrapper;
//import gphase.graph.SpliceGraph;

/**
 * 
 * 
 * @author micha
 */
public class Graph implements Serializable {
	
	static final long serialVersionUID = 6025857021729053707L;
	HashMap speciesHash= null;	// maps EnsemblPfx to Species
	public static final String GRAPH_ENCODE_REF_GENES= "check_enc.fasta";
	public static String overrideSequenceDirPath= null;  
	
	static HashMap buildMap= new HashMap(3);	// maps species builds (standard source) to date
	/*
	 * The Ensembl database names have the format:
	 * 
	 * <species><database type><release number><data version>
	 * 
	 * As an example, the name of the Human Core database for the Ensembl 32 release would be:
	 * homo_sapiens_core_32_35e
	 * A letter suffix to the data version indicates a change in data without a change in assembly; e.g. a new gene build.
	 * 
	 * --
	 * 
	 * Basically there are 3 identifiers for a genome version:
	 * NCBI build (e.g., 35), trivial name (hg17), date (200411)
	 * (which then are mapped to a Ensembl-version in ensemblDB)
	 */
//	static {
//		
//		HashMap map;
//		for (int i = 0; i < EnsemblDBAdaptor.SPECIES.length; i++) {
//			if (EnsemblDBAdaptor.SPECIES[i].contains("homo_sapiens")) {
//				map= new HashMap(1);
//				map.put(new Integer(35), new Integer(200405));	// NCBI build 35 = golden path hg17
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			}else if (EnsemblDBAdaptor.SPECIES[i].contains("mus_musculus")) {
//				map= new HashMap(2);
//				map.put(new Integer(33), new Integer(200405));	// NCBI build 33 = golden path mm5 (200405)
//				map.put(new Integer(34), new Integer(200503));	// NCBI build 34 = golden path mm6
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			} else if (EnsemblDBAdaptor.SPECIES[i].contains("rattus_norvegicus")) {
//				map= new HashMap(1);
//				map.put(new Integer(34), new Integer(200411));	// RGSC build 3.45 (no golden path)
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			} else if (EnsemblDBAdaptor.SPECIES[i].contains("canis_familiaris")) {
//				map= new HashMap(1);
//				map.put(new Integer(1), new Integer(200407));	// UCSC
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			} else if (EnsemblDBAdaptor.SPECIES[i].contains("gallus_gallus")) {
//				map= new HashMap(1);
//				map.put(new Integer(1), new Integer(200402));	// UCSC
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			} else if (EnsemblDBAdaptor.SPECIES[i].contains("pan_troglodytes")) {
//				map= new HashMap(1);
//				map.put(new Integer(3), new Integer(200311));	// UCSC
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			}
//		}
//		
//	}

	public void initTU() {
		Iterator iter= speciesHash.values().iterator();
		while (iter.hasNext()) {
			Species spe= (Species) iter.next();
			spe.initTU();
		}
	}
	
	public static String[] decodeStableID(String stableID) {

		Pattern patty= Pattern.compile("(\\D*)(\\d*)");	// speciesID, stableID
		Matcher matty= patty.matcher(stableID);
		if (!matty.matches())
			System.out.println(stableID);
		
		return new String[] {matty.group(1).substring(0,matty.group(1).length()- 1), matty.group(2)};
	}		
		
	public static String decodeStableID(String stableID, boolean number) {
	
		String[] both= decodeStableID(stableID);
		
		if (number) 
			return both[1];
		else
			return both[0];
	}

	public Graph() {
		speciesHash= new HashMap();
	}

	public int countNonCodingLoci() {
		Iterator iter= speciesHash.values().iterator();
		int cnt= 0;
		while (iter.hasNext()) {
			Species spe= (Species) iter.next();
			cnt+= spe.countNonProteinCodingLoci();
		}
		return cnt;
	}
	public Graph(int nbSpecies) {
		speciesHash= new HashMap(nbSpecies);
	}
	
	public Graph(Species[] newSpecies) {
		this(newSpecies.length);
		addSpecies(newSpecies);
	}


	public static void main(String[] args) {

	}
	
	public Gene getGene(String stableGeneID, String speciesID) {
	
		Species spec= (Species) speciesHash.get(speciesID);		// get species
		if (spec== null)
			return null;

		return spec.getGene(stableGeneID);	
	}

	public Gene[] getGenes() {
		
		Iterator iter= speciesHash.values().iterator();
		Gene[] result= new Gene[0];
		while (iter.hasNext()) {
			Gene[] tmp= ((Species) iter.next()).getGenes();
			Gene[] old= result;
			result= new Gene[old.length+ tmp.length];
			for (int i = 0; i < old.length; i++) 
				result[i]= old[i];
			for (int i = 0; i < tmp.length; i++) 
				result[old.length+ i]= tmp[i];
		}
		
		return result;	
	}
	
	public Exon[] getExons(int region) {
		Vector<Exon> v= new Vector<Exon>();
		Gene[] ge= getGenes();
		for (int i = 0; i < ge.length; i++){
            //v= (Vector) ArrayUtils.addAll(v, ge[i].getExons(region));
            Collections.addAll(v, ge[i].getExons(region));
        }

		return (Exon[]) ArrayUtils.toField(v);
	}

	public Species getSpeciesByEnsemblPrefix(String speciesPfx) {
		
		if (Species.isValidEnsemblPrefix(speciesPfx))
			return (Species) speciesHash.get(speciesPfx);
		
		return getSpeciesByName(Species.SP_NAMES_BINOMIAL[11]);	// tetraodon, genoscan
	}
	
	
	/**
	 * Inefficient - traverses hash via loop. 
	 * Species names are to be given in the format "mus_musculus". 
	 * (Warning, contains an ineficient iteration over all species.)
	 * @param binomialName
	 * @return
	 */
	public Species getSpeciesByName(String binomialName) {
		
		if (speciesHash== null)
			return null;
		
		Iterator iter= speciesHash.values().iterator();
		Species spec;
		while (iter.hasNext()) {		// iterate hashmap
			spec= (Species) iter.next();
			String s= spec.getBinomialName();
			if (s.equalsIgnoreCase(binomialName))	// "mus_musculus"
				return spec;
		}
		
		return null;
	}	
	
	public int getTranscriptCount() {
		if (speciesHash== null) 
			return 0;
		Species[] ge= getSpecies();
		int cnt= 0;
		for (int i = 0; i < ge.length; i++) 
			cnt+= ge[i].getTranscriptCount();
		return cnt;
		
	}
	
	public Species getSpeciesByGeneID(String stableGeneID) {
		
			String[] idDec= decodeStableID(stableGeneID);		
			
			return getSpeciesByEnsemblPrefix(idDec[0]);		// get species
	}	
	
	
	public String countGenesTranscriptsExons() {
		
		String result= "[G,T,E]: ";
		Iterator specIter= speciesHash.values().iterator();
		while(specIter.hasNext()) {
			Iterator geneIter= ((Species) specIter.next()).geneHash.values().iterator();
			int geneCount= 0, transCount= 0, exonCount= 0;
			while(geneIter.hasNext()) {
				Gene gene= (Gene) geneIter.next();
				if (gene.getTranscripts()!= null)
					transCount+= gene.getTranscripts().length;
				if (gene.getExons()!= null)
					exonCount+= gene.getExons().length;
				++geneCount;
			}
			result+= "("+geneCount+","+transCount+","+exonCount+") ";
		}
		
		return result;
	}
	
	
	public boolean addGene(Gene newGene) {
		
			// get speciesHash
		Species spec= (Species) speciesHash.get(newGene.getSpecies().getCommonName()); 
		return spec.addGene(newGene);
	}
	
	/**
	 * 
	 * @param newSpecies
	 * @return <code>true</code> if new <code>Species[]</code> has been added successfully
	 */
	public boolean addSpecies(Species newSpecies) {
		
		if (newSpecies== null)
			return false;
		
		if (speciesHash== null)
			speciesHash= new HashMap();
		
		String key= newSpecies.getCommonName();
		if (speciesHash.get(key)!= null)
			return false;	// return false if key already occupied
		
		speciesHash.put(key, newSpecies);
		return true;
	}
	public Exon[] getExons() {

		Vector v= new Vector();
		Iterator iter= speciesHash.values().iterator();
		while(iter.hasNext()) {
			Species sp= (Species) iter.next();
			Exon[] ex= sp.getExons();
			for (int i = 0; i < ex.length; i++) {
				v.add(ex[i]);
			}
		}
		return (Exon[]) ArrayUtils.toField(v);
	}

	/**
	 * @return
	 */
	public Species[] getSpecies() {

		Collection c= speciesHash.values();
		Species[] sp= new Species[c.size()];
		Iterator it= c.iterator();
		int i= 0;
		while (it.hasNext())
			sp[i++]= (Species) it.next();
		return sp;
	}
	
	public static String readSequence(Gene g) {
		
		String seq= readSequence(
				g.getSpecies(),
				g.getChromosome(),
				g.isForward(),
				g.getStart(),
				g.getEnd()
		);
		
		if (!g.isForward()) 
			seq= invertSequence(seq);
		return seq;
	}
	
	public static String readSequence(SpliceSite s) {
		String seq= readSequence(
				s.getGene().getSpecies(),
				s.getGene().getChromosome(),
				s.getGene().isForward(),
				s.getPos()- SpliceSite.SPLICE_SITE_FLANK_SIZE,
				s.getPos()+ SpliceSite.SPLICE_SITE_FLANK_SIZE
		);
		
		return seq;
	}
	
	public static String readSequence(SpliceSite s, int flank5, int flank3) {
		int posStart= s.getPos();
		int posEnd= s.getPos();
		if (s.isDonor()) {
			++posStart;
			posEnd+= 2;
		} else if (s.isAcceptor()) {
			posStart-= 2;
			--posEnd;
		}
		
			// make all here
//		posStart-= flank5;
//		posEnd+= flank3;
//		
//			// normalize here already for correctly extending end to 1st not read position
//		if (!s.getGene().isForward()) {
//			int h= posStart;
//			posStart= posEnd;
//			posEnd= h;
//		}
//		++posEnd;	// extend to first nonread position
			
		String seq= readSequence(
				s.getGene().getSpecies(),
				s.getGene().getChromosome(),
				s.getGene().isForward(),
				posStart- flank5,
				posEnd+ flank3
		);
		return seq;
	}
	
	public static String readSequence(Exon e) {
		
		String seq= readSequence(
				e.getGene().getSpecies(),
				e.getGene().getChromosome(),
				e.getGene().isForward(),
				e.getStart(),
				e.getEnd()
		);
		
		if (!e.getGene().isForward()) 
			seq= invertSequence(seq);
		return seq;
	}

	public static String readSequence(Transcript t) {
		
		String seq= readSequence(
				t.getGene().getSpecies(),
				t.getGene().getChromosome(),
				t.getGene().isForward(),
				t.getStart(),
				t.getEnd()
		);
		
		if (!t.getGene().isForward()) 
			seq= invertSequence(seq);
		return seq;
	}
	
	public static String readSequence(DirectedRegion reg) {
		
		String seq= readSequence(
				reg.getSpecies(),
				reg.getChromosome(),
				reg.isForward(),
				reg.getStart(),
				reg.getEnd()
		);
		
			// neg strand already reversed in subroutine
		return seq;
	}

	public static String invertSequence(String seq) {
		
		seq= reverseSequence(seq);
		seq= complementarySequence(seq); 
		return seq;
	}
	
	public static String reverseSequence(String seq) {
		
		StringBuffer sb= new StringBuffer(seq.length());
		for (int i = (seq.length()-1); i >= 0; i--) 
			sb.append(seq.charAt(i));
		return sb.toString();
	}

    public static char complementaryCharacter(char c) {
        boolean wasLow= Character.isLowerCase(c);
        c= Constants2.NA_COMPL_IUPAC[Character.toUpperCase(c)- 65];
        if (wasLow)
            c= Character.toLowerCase(c);
        return c;
    }

	public static String complementarySequence(String seq) {
		
		StringBuffer sb= new StringBuffer(seq.length());
		for (int i = 0; i < seq.length(); i++) {
			
			char c= seq.charAt(i);
            c= complementaryCharacter(c);
			sb.append(c);
		}
		return sb.toString();
	}

    public static String readSequence(Species spe, CharSequence chromosome, boolean forwardStrand,
            long start, long end){
        return readSequence(spe, chromosome, forwardStrand, start, end, false);
    }

	/**
	 * 
	 * @param spe
	 * @param chromosome
	 * @param forwardStrand
	 * @param start 1st position to be read
	 * @param end	1st position to be read
     * @param isCircular circular chromosome
	 * @return
	 */
	public static String readSequence(Species spe, CharSequence chromosome, boolean forwardStrand, 
			long start, long end, boolean isCircular)
				throws RuntimeException {
			
			if (!forwardStrand) {	// WAS: (start< 0), neg strand genes
				start= -start;
				end= -end;
				if (start> end) {	// WAS: for both strands, but shouldnt incur on + strand
					long h= start;
					start= end;
					end= h;
				}
			}
            /*
            Thasso: 25.5.11 reset start if we are outside of the bounds
             */
		    if(start < 0){
                start = 0;
            }else{
			    start--;	// this is ok
            }
			byte[] seq= new byte[(int) (end- start)];
			String s= null;
			long p= -1;
			try {
//				System.out.println(getSequenceDirectory(speRealName)+ File.separator+ "chr"+ chromosome+ Constants2.CHROMOSOME_EXT);
//				System.out.println(start+"-"+end);				
//				if (spe== null)
//					return null; // no sequence available
				
				RandomAccessFile raf= getRAF(spe, chromosome);
				
				// this has to be done properly, with the file linebreak-length x
				// offset+x+start+(start/line)*x
				// p= offset+ 1+ start+ (start/line);
				// 100215: should be proper now
				String pfx= null, sfx= null;
				if (start< 0) {		// circular genomes
					//System.err.println("Neg seek: "+forwardStrand+", "+start+", "+end);
					pfx= readSequence(spe, chromosome, forwardStrand, chrLen+ start, chrLen);	// start< 0
					start= 0;
				}
				if (end> chrLen) {
					sfx= readSequence(spe, chromosome, forwardStrand, 1, (end- chrLen)+ 1);
					end= chrLen;
				}
				
				p= headerOffset+ fileSep.length()+ start+ ((start/lineLen)* fileSep.length());
				assert(p>= 0);
				raf.seek(p);	// fpointer is different from reading point!
				int pos= 0;
				int nextN= (int) (lineLen- (start% lineLen));				// read (end of) first line
				while (pos+ nextN< seq.length&& p+pos+nextN< raf.length()) {		// full lines
					raf.readFully(seq,pos,nextN);
					raf.skipBytes(1);
					pos+= nextN;
					nextN= lineLen;
				}
				int a= seq.length- pos;
				int b= (int) (raf.length()- p- pos- 2);
				int rest= Math.min(a, b);	// catch EOF (when reading range larger than file)
				if (a < 0)
					rest= b;
				else if (b< 0)
					rest= a;
				if (a< 0&& b< 0)
					rest= 0;
				try {
                    if(pos+p < raf.length())
					    raf.readFully(seq,pos,rest);	// read start of last line
				} catch (Exception e) {	//EOFException, IndexOutOfBoundsException
					System.err.println("Problems reading "+chromosome+": "+(p+pos)+", "+rest+"> "+raf.length()+" into "+seq.length+": "+e.getMessage());
					System.err.println("check for the right species/genome version!");
					e.printStackTrace();
					return s;
				}

				s= new String(seq);
				if (seq.length- pos> rest)
					s= s.substring(0, s.length()- ((seq.length- pos)- rest));
				if (!forwardStrand) 
					s= StringUtils.reverseComplement(s);

                s = s.trim();
				
				if (isCircular && pfx!= null) {
					if (forwardStrand)
						s= pfx+ s;
					else
						s+= pfx;
				}
				if (isCircular && sfx!= null) {
					if (forwardStrand)
						s+= sfx;
					else
						s= sfx+ s;
				}
				
			} catch (Exception e) {
				throw new RuntimeException("Problems reading sequence " +
                        chromosome+": pos "+p+", len "+seq.length+ ",\n" +
                        "check whether chromosomal sequence exists / has the correct size",
						e);
				//e.printStackTrace();
			}
			
			return s;
		}

	private static CharSequence lastChr= null;
	private static int lineLen= -1, headerOffset= -1;
	private static long chrLen= -1;	// make int !!!
	private static String fileSep= null;
	private static RandomAccessFile raf= null;
	private static RandomAccessFile readChromosome_old(CharSequence chromosome) {
		
		String dirPath= getSequenceDirectory(null); 
		String fName= dirPath+File.separator+ chromosome+ Constants2.CHROMOSOME_EXT;
		if (!new File(fName).exists()) {	// try to find file case insensitive
			String[] list= new File(dirPath).list();
			String chrFile= chromosome+ Constants2.CHROMOSOME_EXT;
			for (int i = 0; list!= null&& i < list.length; i++) 
				if (list[i].equalsIgnoreCase(chrFile)) {
					chrFile= list[i];
					break;
				}
			fName= dirPath+ File.separator+ chrFile;
		}
		File f= new File(fName);
		try {
			BufferedInputStream buffy= new BufferedInputStream(new FileInputStream(f));
			
			// init file characteristics
			byte[] buf= new byte[1024];
			fileSep= null;
			boolean tag= false;
			int totByte= 0, readByte= 0;
			while (fileSep== null|| !tag) {
				readByte= buffy.read(buf, 0, buf.length);
				totByte+= readByte;
				for (int i = 0; i < readByte; i++) {
					if (buf[i]== '>')
						tag= true;
					else if (buf[i]== '\n'|| buf[i]== '\r') {
						fileSep= Character.toString((char) buf[i]);
						if (i+ 1< readByte) {
							if (buf[i+ 1]== '\n'|| buf[i+ 1]== '\r') 
								fileSep+= Character.toString((char) buf[++i]);
						} else {
							buf[0]= (byte) buffy.read();
							if (buf[0]== '\n'|| buf[0]== '\r')
								readByte= 0;
							else {
								readByte= 1;
							}
							break;
						}
					} else if (tag&& fileSep!= null) {
						readByte-= i- 1;
						totByte-= i- 1;
						System.arraycopy(buf, i, buf, 0, readByte);
						break;
					}
				}
			}
			headerOffset= totByte;
			
			// read first fasta line
			int mark= readByte;
			lineLen= -1;
			readByte= buffy.read(buf, readByte, buf.length- readByte);
			for (int i = mark; i < buf.length; ++i) {
				if (buf[i]== '\n'|| buf[i]== '\r') {
					lineLen= i- mark;
					break;
				}
			}
			assert(lineLen> 0);
			chrLen= f.length()- headerOffset- 2* fileSep.length();	// fsep 1st and last line
			chrLen= chrLen- ((chrLen/ (lineLen+ fileSep.length()))* fileSep.length());
			chrSeq= new byte[(int) chrLen];
			for (int i = 0; i < buf.length; i++) {
				
			}
			System.arraycopy(buf, 0, chrSeq, 0, lineLen);
			mark= lineLen;
			
			// read all lines
			
			
			

			
			lastChr= chromosome;
			return raf;
		} catch (Exception e) {
			e.printStackTrace();
		}

		return raf;
	}

	public static String getSequenceDirectory(Species spe) {

		if (overrideSequenceDirPath!= null)
			return overrideSequenceDirPath;
		
		// Species spe= getSpeciesByName(realName);
		if (spe== null)
			return null;
		String bin= spe.getBinomialName();
		int binSep= bin.indexOf('_');
		String bin1= bin.substring(0, binSep);
		String bin2= bin.substring(binSep+1, bin.length());
//		HashMap hm= (HashMap) buildMap.get(realName);
//		Object o= hm.get(new Integer(spe.getBuildVersion()));
//		int dateID= ((Integer) ((HashMap) buildMap.get(realName))
//						.get(new Integer(spe.getBuildVersion()))).intValue();		// extract ID (date of build)
		
		
		String seqDirName= Character.toUpperCase(bin.charAt(0))+ "."+ bin2;
		File speciesGenome=  new File(Constants2.DATA_DIR+ File.separator
				+ Constants2.SEQUENCES_SUBDIR+ File.separator
				+ seqDirName);
		String[] list= speciesGenome.list();
		if (list== null)
			System.err.println("\nSequence Directory not found: "+ speciesGenome.getAbsolutePath()+"\nplease provide build version: "+spe.getAnnotationVersion());
		
//		String oneLetterPFX= ""+ bin1.charAt(0)+ bin2.charAt(0);
//		String threeLetterPFX= ""+ bin1.substring(0,3)+ bin2.substring(0,3);
//		String oneThreeLetterPFX= ""+ bin1.charAt(0)+ bin2.substring(0,3);
//		Vector combiVec= new Vector();
//		String x= spe.getAnnotationVersion();
//		for (int i = 0; i< Species.ASSEMBLY_PFX.length; i++) 
//			combiVec.add(Species.ASSEMBLY_PFX[i].toUpperCase()+x);
//		combiVec.add(oneLetterPFX.toUpperCase()+x);
//		combiVec.add(threeLetterPFX.toUpperCase()+x);
//		combiVec.add(oneThreeLetterPFX.toUpperCase()+x);
		int i;
		String subdirName= null;
		for (i = 0; i < list.length; i++) {
//			int j;
//			for (j = 0; j < combiVec.size(); j++) {
//				String s= (String) combiVec.elementAt(j);
//				int p= list[i].toUpperCase().indexOf(s);
//				if (p>= 0&& (p+ s.length()== list[i].length()|| list[i].charAt(p+s.length())== '_'))	// dateID
//					break;
//			}
//			if (j< combiVec.size())
//				break;
			StringTokenizer st= new StringTokenizer(list[i], "_.");
			while (st.hasMoreTokens()) {
				if (st.nextToken().toUpperCase().contains(spe.getGenomeVersion().toUpperCase())) {
					subdirName= list[i];
					break;
				}
			}
			if (subdirName!= null)
				break;
		}
//		if (i> list.length- 1) {
		if (subdirName== null) {
			System.err.println("Not found genomeVer "+spe.getGenomeVersion()+".");
//			System.err.print("Build not found: "+ bin+" ( ");
//			for (int j = 0; j < combiVec.size(); j++) 
//				System.out.print(combiVec.elementAt(j));
//			System.out.println(")");
		}		
		
		return speciesGenome.getAbsolutePath()+ File.separator
				+ list[i]+ File.separator+ Constants2.CHROMOSOME_DIR;
	}

	public static String getSequenceDirectory_old(String realName) {
	
		String seqDirName= Character.toUpperCase(realName.charAt(0))+ "."
				+ realName.substring(realName.indexOf('_')+ 1);
		File speciesGenome= new File(Constants2.DATA_DIR+ File.separator
				+ Constants2.SEQUENCES_SUBDIR+ File.separator
				+ seqDirName);
		Pattern patty= Pattern.compile("^(\\D+)(\\d+).*");	// assuming that eg mouse "...mm5" is not relevant
		Matcher matty;
		int highestVersion= -1;
		String goldenPath= null;
		String[] list= speciesGenome.list();
		for (int i = 0; i < list.length; i++) {
			matty= patty.matcher(list[i]);
			if (!matty.matches())
				continue;
			if (!matty.group(1).equals("golden_path_"))
				continue;
			int ver= Integer.parseInt(matty.group(2));
			if (ver> highestVersion) {
				highestVersion= ver;
				goldenPath= list[i];
			}
		}
		
		return speciesGenome.getAbsolutePath()+ File.separator
				+ goldenPath+ File.separator+ Constants2.CHROMOSOME_DIR;
	}

	
	public void init() {
		//initExonSplicePattern();
///		initExonHomology();
//		filter(g);
//		System.out.println("Graph filtered ---");
//		
//		g.init();
//		System.out.println("Graph inited ---");
	}
	
	void initExonSplicePattern() {
//		for (int i = 0; i < getSpecies().length; i++) {
//			for (int j = 0; j < getSpecies()[i].getGenes().length; j++) {
//				Gene g= getSpecies()[i].getGenes()[j];
//				switch (g.getTranscripts().length) {
//					case 0: System.err.println("Not transcribed gene "+ g.getStableID()
//							+ " ("+ getSpecies()[i].getRealName()+")");
//							break;
//					
//					case 1: for (int k = 0; k < g.getExons().length; k++) 
//								g.getExons()[k].addIdentica(Exon.TYPE_CO);
//							break;
//				
//					default: {
//						int nbTranscripts= g.getTranscripts().length;
//						
//					}
//				}
//			}
//		}
	}
	
	
	/**
	 * @deprecated intransparent loop structure
	 *
	 */
	void initExonHomology_old() {
			
/*		Gene[] as= getSpecies()[0].getGenes();	 
		for (int i = 0; i < as.length; i++) {
			System.out.println(i+": initing exon homology of "+as[i].getStableID());

				// gfx progress
			int sum= 0;
			for (int k = 0; k < as[i].getHomologies().length; k++) {
				for (int j = 0; j < as[i].getHomologies()[k].getExons().length; j++) {				// align pw
					sum+= as[i].getHomologies().length- k;
				}				
			}
			char[] c= new char[sum];
			Arrays.fill(c, '.');
			System.out.println(c);

				// pw alignment of exons
			for (int k = 0; k < as[i].getHomologies().length; k++) {
				
				for (int j = 0; j < as[i].getHomologies()[k].getExons().length; j++) {				// align pw

					for (int ii = 0; ii < as[i].exons.length; ii++) {				// with base gene 					
						
						String[] seqs= new String[] {readSequence(as[i].exons[ii]), readSequence(as[i].getHomologies()[k].getExons()[j])};
						String[] names= new String[] {as[i].exons[ii].getStableID(), as[i].getHomologies()[k].getExons()[j].getStableID()};
						
						PWHit hit= new PWHit(as[i].exons[ii], as[i].getHomologies()[k].getExons()[j]);
						ClustalWrapper cw= ((ClustalWrapper)AlignmentGenerator.alignClustal(AlignmentGenerator.writeOutTemp(names, seqs)));
						hit.setScore(cw.getScore());
						hit.setAlignment(cw.getLayout());
						as[i].exons[ii].addHit(as[i].getHomologies()[k], hit);
						as[i].getHomologies()[k].getExons()[j].addHit(as[i], hit);
					}
					
					System.out.print('*');		// gfx output
					System.out.flush();
					
					for (int kk = (k+1); kk < as[i].getHomologies().length; kk++) {
						
						for (int ii = 0; ii < as[i].getHomologies()[kk].exons.length; ii++) {				// with other homologs 					
							
							String[] seqs= new String[] {readSequence(as[i].getHomologies()[kk].exons[ii]), readSequence(as[i].getHomologies()[k].getExons()[j])};
							String[] names= new String[] {as[i].getHomologies()[kk].exons[ii].getStableID(), as[i].getHomologies()[k].getExons()[j].getStableID()};
							
							PWHit hit= new PWHit(as[i].getHomologies()[kk].exons[ii], as[i].getHomologies()[k].getExons()[j]);
							ClustalWrapper cw= ((ClustalWrapper)AlignmentGenerator.alignClustal(AlignmentGenerator.writeOutTemp(names, seqs)));
							hit.setScore(cw.getScore());
							hit.setAlignment(cw.getLayout());
							as[i].getHomologies()[kk].exons[ii].addHit(as[i].getHomologies()[k], hit);
							as[i].getHomologies()[k].getExons()[j].addHit(as[i].getHomologies()[kk], hit);
						}
						System.out.print('*');
						System.out.flush();
					}
					
				}
			}
			System.out.println();
		}
		

		// output
		for (int i = 0; i < as.length; i++) {
			for (int j = 0; j < as[i].getExons().length; j++) { 
				for (int k = 0; k < as[i].getHomologies().length; k++) {
					// System.out.print(i+ " - "+as[i].getExons()[j].getStableID()+" x "+ as[i].getHomologs()[k].getStableID()+ " : ");
					PWHit[] bestHits= getBRH(as[i].getExons()[j], as[i].getHomologies()[k], true);
					if (bestHits.length== 1) {							// add UBRH
						Exon e1= (Exon) bestHits[0].getObject1();
						Exon e2= (Exon) bestHits[0].getObject2();
						e1.addHomolog(e2.getGene(), e2);
						e2.addHomolog(e1.getGene(), e1);
					
					} else if (bestHits.length> 1)
						System.err.println("MBRHs ("+bestHits.length+") -"+i+ "- "+as[i].getExons()[j].getStableID()+" x "+ as[i].getHomologies()[k].getStableID());
					else 	// < 1, one exon not found
						System.err.println("noBRHs -"+i+ "- "+as[i].getExons()[j].getStableID()+" ("+as[i].getExons().length+") x "+ 
								as[i].getHomologies()[k].getStableID()+ " ("+ as[i].getHomologies()[k].getExons().length+")");
						
//					System.out.print(getBRH(as[i].getExons()[j], as[i].getHomologs()[k], true).length);
//					System.out.println(" BRH / "+ as[i].getExons()[j].getHits(as[i].getHomologs()[k]).length+ " hits.");
				}
			}
		}
*/		
		
	}
	
	public void repairAlignmentErrors() {
		Gene[] ge= getGenes();
		for (int i = 0; i < ge.length; i++) {
			ge[i].repairAlignmentErrors();
		}
	}
	
	/**
	 * 
	 * @param newSpecies
	 * @return <code>true</code> if new <code>Species[]</code> has been added successfully
	 */
	public boolean addSpecies(Species[] newSpecies) {
		
		if (newSpecies== null)
			return false;
		
		boolean result= true;
		for (int i = 0; i < newSpecies.length; i++)  
			result&= addSpecies(newSpecies[i]);
		
		return result;
	}

	private static byte[] readChromosome(CharSequence chromosome) {
		
		String dirPath= getSequenceDirectory(null); 
		String fName= dirPath+File.separator+ chromosome+ Constants2.CHROMOSOME_EXT;
		if (!new File(fName).exists()) {	// try to find file case insensitive
			String[] list= new File(dirPath).list();
			String chrFile= chromosome+ Constants2.CHROMOSOME_EXT;
			for (int i = 0; list!= null&& i < list.length; i++) 
				if (list[i].equalsIgnoreCase(chrFile)) {
					chrFile= list[i];
					break;
				}
			fName= dirPath+ File.separator+ chrFile;
		}
		File f= new File(fName);
		fileSep= "\n"; // TODO: Factory FileHelper.guessFileSep(f);
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(f));
			headerOffset= 0;
			String line= buffy.readLine();
			while (!line.startsWith(">")) {
				headerOffset+= line.length()+ fileSep.length();
				line= buffy.readLine();
			}
			headerOffset+= line.length()+ fileSep.length();
			
			// read first fasta line
			line= buffy.readLine();
			lineLen= line.length();
			buffy.close();
			
			chrLen= f.length()- headerOffset- fileSep.length();	// fsep 1st and last line
			chrLen= chrLen- ((chrLen/ (lineLen+ fileSep.length()))* fileSep.length());
			chrSeq= new byte[(int) chrLen];
			BufferedInputStream in= new BufferedInputStream(new FileInputStream(f));
			in.skip(headerOffset);
			byte[] buf= new byte[lineLen];
			int p= 0, x= 0, s= fileSep.length(); // dbgLctr= 0;
			for (x= in.read(buf); x>= 0; x= in.read(buf)) {
				//++dbgLctr;
				if (x< buf.length|| buf[buf.length- 1]== '\n'|| buf[buf.length- 1]== '\r')	// hg18 chr2 last line 49+ nl
					x-= fileSep.length();
				System.arraycopy(buf, 0, chrSeq, p, x);
				p+= x;
				in.skip(s);
			}
			in.close();
			
			lastChr= chromosome.toString();
			return chrSeq;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	
		return chrSeq;
	}

	private static RandomAccessFile getRAF(Species spe, CharSequence chromosome) {
		
		if (raf == null|| !chromosome.equals(lastChr)) {
			String dirPath= getSequenceDirectory(spe); 
			String fName= dirPath+File.separator+ chromosome+ Constants2.CHROMOSOME_EXT;
			if (!new File(fName).exists()) {	// try to find file case insensitive
				String[] list= new File(dirPath).list();
				String chrFile= chromosome+ Constants2.CHROMOSOME_EXT;
				for (int i = 0; list!= null&& i < list.length; i++) 
					if (list[i].equalsIgnoreCase(chrFile)) {
						chrFile= list[i];
						break;
					}
				fName= dirPath+ File.separator+ chrFile;
			}

            if (!new File(fName).exists()) {	// still no match ... check if there is just a single file
                File[] list= new File(dirPath).listFiles();
                if(list.length == 1){
                    fName = list[0].getAbsolutePath();
                }
            }

            if (!new File(fName).exists()) {	// still nothig ... repoert error...
                throw new RuntimeException("Chromosome file not found! The genome hast to be split into chromosome fasta files, i.e. Chr1.fa. \n" +
                        "I was looking for "+chromosome+ Constants2.CHROMOSOME_EXT);
            }



			File f= new File(fName);
			fileSep= "\n";// TODO: Factory: FileHelper.guessFileSep(f);
			try {
				if (raf!= null)
					raf.close();
				raf= new RandomAccessFile(f, "r");
				
				// init file characteristics
				String read= raf.readLine();
				while (!read.startsWith(">"))
					read= raf.readLine();
				headerOffset= read.length();
				read= "";
				while (read.trim().length()< 1)
					read= raf.readLine();
				lineLen= read.length();
				
				chrLen= raf.length()- headerOffset- 2* fileSep.length();	// fsep 1st and last line
				chrLen= chrLen- ((chrLen/ (lineLen+ fileSep.length()))* fileSep.length());
	
				
				lastChr= chromosome.toString();
				return raf;
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	
		return raf;
	}

	/**
			 * 
			 * @param spe
			 * @param chromosome
			 * @param forwardStrand
			 * @param start 1st position to be read
			 * @param end	1st position not to be read
			 * @return
			 */
			public static void readSequence(Species spe, CharSequence chromosome, boolean forwardStrand,
					long start, long end,
					ByteArrayCharSequence cs, int from, int to) 
						throws RuntimeException {
					
					if (start< 0|| (end- start+ 1> 1000))
						System.currentTimeMillis();
				
					if (end- start+ 1!= to- from)
						throw new RuntimeException("ByteArrayCharSequence too small: "+ (to- from)+ " needs "+ (end- start+ 1));
					
					if (!forwardStrand) {	// WAS: (start< 0), neg strand genes
						start= -start;
						end= -end;
						if (start> end) {	// WAS: for both strands, but shouldnt incur on + strand
							long h= start;
							start= end;
							end= h;
						}
					}
					start--;	// this is ok
					// now cs
					//byte[] seq= new byte[(int) (end- start)];
					//String s= null;
					
					long p= -1;
					try {
		//				System.out.println(getSequenceDirectory(speRealName)+ File.separator+ "chr"+ chromosome+ Constants2.CHROMOSOME_EXT);
		//				System.out.println(start+"-"+end);				
		//				if (spe== null)
		//					return null; // no sequence available
						
						RandomAccessFile raf= getRAF(spe, chromosome);
						
						// this has to be done properly, with the file linebreak-length x
						// offset+x+start+(start/line)*x
						// p= offset+ 1+ start+ (start/line);
						// 100215: should be proper now
						//String pfx= null, sfx= null;
						if (start< 0) {		// circular genomes
							//System.err.println("Neg seek: "+forwardStrand+", "+start+", "+end);
							int diff= (int) (from- start);	// start< 0
							int pfxFrom= forwardStrand? from: to- diff,	// append @start or @end
								pfxTo= forwardStrand? from+ diff: to;
							readSequence(spe, chromosome, forwardStrand, chrLen+ start, chrLen,
									cs, pfxFrom, pfxTo);	// start< 0
							start= 0;
						}
						if (end> chrLen) {
							int diff= (int) (end- chrLen);
							int sfxFrom= forwardStrand? to- diff: from,
								sfxTo= forwardStrand? to: from+ diff;
							readSequence(spe, chromosome, forwardStrand, 1, (end- chrLen)+ 1,
									cs, sfxFrom, sfxTo);
							end= chrLen;
						}
						
						p= headerOffset+ fileSep.length()+ start+ ((start/lineLen)* fileSep.length());
						assert(p>= 0);
						raf.seek(p);	// fpointer is different from reading point!
						int mark= cs.end;
						int nextN= (int) (lineLen- (start% lineLen));				// read (end of) first line
						while (cs.end+ nextN<= to) {		// full lines
							assert(p+ (cs.end- mark)+ nextN< raf.length());
							raf.readFully(cs.chars,cs.end,nextN);
							raf.skipBytes(1);
							cs.end+= nextN;
							nextN= lineLen;
						}
	/*					int a= cs.a.length- pos;
						int b= (int) (raf.length()- p- pos- 1);
						int rest= Math.min(a, b);	// catch EOF (when reading range larger than file)
						if (a < 0)
							rest= b;
						else if (b< 0)
							rest= a;
						if (a< 0&& b< 0)
							rest= 0;
	*/	
						int rest= to- cs.end;	
						if (rest< 0)
							System.currentTimeMillis();
						try {
							raf.readFully(cs.chars,cs.end,rest);	// read start of last line
							cs.end+= rest;
						} catch (Exception e) {	//EOFException, IndexOutOfBoundsException
							System.err.println("Problems reading "+chromosome+": "+ p+ ", "+rest+"> "
									+raf.length()+" into "+cs.length()+": "+e.getMessage());
							System.err.println("check for the right species/genome version!");
							e.printStackTrace();
							return;
						}				
						
	//					s= new String(seq);
	//					if (seq.length- pos> rest)
	//						s= s.substring(0, s.length()- ((seq.length- pos)- rest));
						if (!forwardStrand) 
							ByteArrayCharSequence.reverseComplement(cs, from, to);
						
						
					} catch (Exception e) {
						throw new RuntimeException("Problems reading sequence "+ chromosome+": pos "+p
								+ ", len "+(end- start+ 1)+ " into ["+ cs.start+ ","+ cs.end+ "]\n\t"
								+ e.getMessage());
						//e.printStackTrace();
					}
				}

	/**
		 * 
		 * @param speRealName
		 * @param chromosome
		 * @param forwardStrand
		 * @param start 1st position to be read
		 * @param end	1st position not to be read
		 * @return
		 */
		static byte[] chrSeq= null;
		public static void readSequence(CharSequence chromosome, boolean forwardStrand, 
				long start, long end, ByteArrayCharSequence cs, int from, int to) {
				
				assert(Math.abs(end- start)+ 1== to- from);
				
				if (!forwardStrand) {	// WAS: (start< 0), neg strand genes
					start= -start;
					end= -end;
					if (start> end) {	// WAS: for both strands, but shouldnt incur on + strand
						long h= start;
						start= end;
						end= h;
					}
				}
				start--;	// this is ok
				// now cs
				//byte[] seq= new byte[(int) (end- start)];
				//String s= null;
				
				try {
	//				System.out.println(getSequenceDirectory(speRealName)+ File.separator+ "chr"+ chromosome+ Constants2.CHROMOSOME_EXT);
	//				System.out.println(start+"-"+end);				
	//				if (spe== null)
	//					return null; // no sequence available
					
					
					RandomAccessFile raf= Graph.raf;
					if (!chromosome.equals(lastChr)) {	// chrSeq== null|| , no--tries to read in everytime if chr does not fit in memory
						try {
							chrSeq= readChromosome(chromosome);
							raf= null;
							Graph.raf= null;
						} catch (Throwable t) {
							raf= getRAF(null, chromosome);
						}
					}
					
					// this has to be done properly, with the file linebreak-length x
					// offset+x+start+(start/line)*x
					// p= offset+ 1+ start+ (start/line);
					// 100215: should be proper now
					//String pfx= null, sfx= null;
					if (start< 0) {		// circular genomes
						//System.err.println("Neg seek: "+forwardStrand+", "+start+", "+end);
						start= 0;
/*						NO CIRCULAR GENOMES (Issue 38 and most are non-circular)							
						int diff= (int) (from- start);	// start< 0
						int pfxFrom= forwardStrand? from: to- diff,	// append @start or @end
							pfxTo= forwardStrand? from+ diff: to;
						readSequence(chromosome, forwardStrand, chrLen+ start, chrLen,
								cs, pfxFrom, pfxTo);	// start< 0
						start= 0;
*/						
					}
					if (end> chrLen) {
						end= chrLen;
/*						NO CIRCULAR GENOMES (Issue 38 and most are non-circular)						
						int diff= (int) (end- chrLen);
						int sfxFrom= forwardStrand? to- diff: from,
							sfxTo= forwardStrand? to: from+ diff;
						readSequence(chromosome, forwardStrand, 1, (end- chrLen)+ 1,
								cs, sfxFrom, sfxTo);
						end= chrLen;
*/						
					}
					
					if (raf== null) {
						System.arraycopy(chrSeq, (int) start, cs.chars, from, (int) (end- start)); // 20101210: bugfix
					} else {
						long startP= headerOffset+ fileSep.length()+ start+ ((start/lineLen)* fileSep.length());
						assert(startP>= 0);
						raf.seek(startP);	// fpointer is different from reading point!
						int curr= from;
						int nextN= (int) (lineLen- (start% lineLen));				// read (end of) first line
						int fSepLen= fileSep.length();
						int cnt= 0;
						long p= startP;
						while (curr+ nextN<= to) {		// full lines
							assert(p+ (to- curr)+ nextN< raf.length());
							raf.readFully(cs.chars,curr,nextN);
							p+= nextN;
							raf.skipBytes(fSepLen);
							p+= fSepLen;
							curr+= nextN;
							nextN= lineLen;
							++cnt;
						}
						
						// don't read over end of chromosome
						int rest= (int) Math.min(to- curr, raf.length()- p);	
						try {
							raf.readFully(cs.chars,curr,rest);	// read start of last line
						} catch (Exception e) {	//EOFException, IndexOutOfBoundsException
							
							String msg= "Problems reading chromosome "+chromosome+" from "+ start+ " to "+ end+
									"\ninto array of "+ cs.length()+ " from "+ from+ " to "+ to+
									"\nstarted @ file position "+ startP +", total file length "+ raf.length()+
									"\nread "+cnt+ " batches, first "+ (int) (lineLen- (start% lineLen))+ " chars, then "+ (cnt-1)+ " of "+ lineLen+
									"\nerror occurred @ file position "+ p+
									"\nwhen attempting to copy remaining "+ rest+ " chars @ array position "+ curr;
							throw new RuntimeException(msg, e);
						}
						
					}
					// 20101210: bugfix, changed (to- from) to (end-start): exceeds readlength may exceed end of chromosome
					if ((end- start)< (to- from)) {
						Arrays.fill(cs.chars, (int) (from+ (end- start)), (int) (from+ (to- from)), (byte) 'N');
					}

					if (!forwardStrand) 
						ByteArrayCharSequence.reverseComplement(cs, from, to);
					
				} catch (Exception e) {
					String msg= "Problems reading chromosome "+chromosome+" from "+ start+ " to "+ end+
					"\ninto array of "+ cs.length()+ " from "+ from+ " to "+ to;
					throw new RuntimeException(msg, e);
				}
			}

	int asVariations;
}
