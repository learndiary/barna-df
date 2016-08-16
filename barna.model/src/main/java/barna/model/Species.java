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
 * Created on Mar 31, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package barna.model;

import barna.commons.utils.ArrayUtils;
import barna.model.constants.Constants2;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * 
 * 
 * @author micha
 */
public class Species implements Serializable {
	
	String genomeVersion= null;
	public static final String ENS_ID= "ENS";
	public static final String ENS_IMCB_ID= "SIN";
	public static final String GENOSCOPE_ID= "GS";
	public static final String FLYBASE_ID= "C";
	
	public static final String[] ASSEMBLY_PFX= new String[] {
		"ENSEMBL", "HG", "hg", "NCBI", "golden_path_", "gp", "WASHUC", "JGI", "Zv", "FUGU"
	};
	// org=
	public static final String[][][] SP_GENOME_VERSIONS= new String[][][] {
		new String[][] {new String[] {"GP200405", "NCBI35", "HG17"}, new String[] {"GP200603","NCBI36", "HG18"}},
		new String[][] {new String[] {"GP200603"}},
		new String[][] {new String[] {"GP200603", "MM8", "NCBI36"}},
		new String[][] {new String[] {"GP200411", "RGSC43"}},
		new String[][] {new String[] {"GP200505", "CANFAM2"}},
		new String[][] {new String[] {"GP200503", "BTAU2"}, new String[] {"BTAU31"}},
		new String[][] {new String[] {}},	// opossum
		
		new String[][] {new String[] {"GP200605", "WASHUC2"}},
		new String[][] {new String[] {"GP200508", "JGI41"}},
		new String[][] {new String[] {"GP200603", "ZV6"}},
		new String[][] {new String[] {}},	// fugu
		new String[][] {new String[] {}},	// tetraodon
		
		new String[][] {new String[] {"BDGP43"}},
		new String[][] {new String[] {}},	// mosquito
		new String[][] {new String[] {"GP200501", "APIMEL2"}},
		new String[][] {new String[] {"GP200607", "WS160"}},
		
		new String[][] {new String[] {}},	// yeast
		new String[][] {new String[] {}},	// seasquirt
		new String[][] {new String[] {}},	// seqsquirt2
		new String[][] {new String[] {}},	// platypus
		new String[][] {new String[] {}}	// cress
		
	};

	public static final String[][][] SP_ANNOTATION_VERSIONS= new String[][][] {
		new String[][] {new String[] {"GENCODE200510"}, new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"PANTRO21", "ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42"}, new String[] {"ENSEMBL43"}},	// cow
		new String[][] {new String[] {}},	// opossum

		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {}},	// fugu
		new String[][] {new String[] {}},	// tetraodon
		
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {}},	// mosquito
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		
		new String[][] {new String[] {}},
		new String[][] {new String[] {}},
		new String[][] {new String[] {}},
		new String[][] {new String[] {}},
		new String[][] {new String[] {}}

	};
	public static final String[] SP_UCSC_CGI_STRINGS= new String[] {
		//&clade=vertebrate&org=Mouse&db=mm8
		// org=

			"Homo_sapiens;db=hg17",	// &db=hg18 newest, but for gencode..
			"Chimp",	// db=panTro
			"Mouse",	// db=mm8
			"Rat",		// &db=rn4
			"Dog",		// &db=canFam2
			"Cow",		//&db=bosTau2
			"Opossum",	//&db=monDom4
			
			"Chicken",	//&db=galGal3
			"X.+tropicalis",	//&db=xenTro2
			"Zebrafish",	//&db=danRer4
			"Fugu",			//&db=fr1
			"Tetraodon",	//&db=tetNig1
			
			"D.+melanogaster",	// db=dm2
			"A.+gambiae",	//&db=anoGam1
			"A.+mellifera",  //&db=anoGam1	// not in ensembl
			"C.+elegans",	//&db=ce2
			
			"S.+cerevisiae",	//&db=sacCer1
			"",
			"",
			"",
			""
			
	};
	
	public static final String[][] SP_UCSC_GENOME_VERSIONS= new String[][] {
		//&clade=vertebrate&org=Mouse&db=mm8
		// org=

		new String[] {"hg17", "hg18"},	// &db=hg18 newest, but for gencode..
		new String[] {"panTro"},
		new String[] {"mm8"},
		new String[] {"rn4"},
		new String[] {"canFam2"},
		new String[] {"bosTau2"},
		new String[] {"monDom4"},
		
		new String[] {"galGal3"},
		new String[] {"xenTro2"},
		new String[] {"danRer4"},	// zfish
		new String[] {"fr1"},	// fugu
		new String[] {"tetNig1"},	// tetraodon
		
		new String[] {"dm2"},	// droso
		new String[] {"anoGam1"},	// mosquito
		new String[] {"apiMel2"},  	// honeybee
		new String[] {"ce2"},	// worm
		
		new String[] {"sacCer1"},	// yeast
		new String[] {"ci2"},	// ciona
		new String[] {""},
		new String[] {""},
		new String[] {""}
		
		
};
	
	public static final String[] SP_NAMES_ENS_PFX= new String[] {
			ENS_ID+			"", 
			ENS_ID+			"PTR",
			ENS_ID+			"MUS", 
			ENS_ID+			"RNO",
			ENS_ID+			"CAF",
			ENS_ID+			"BTA",
			ENS_ID+			"MOD", 
			
			ENS_ID+			"GAL",
			ENS_ID+			"XET",
			ENS_ID+			"DAR",
			ENS_IMCB_ID+	"FRU",	
			GENOSCOPE_ID+ 	"TEN", // GSTEN0039005408, HOXA1, HOXA10, HOXBb1, ...
			
			FLYBASE_ID, // CR for genes on forward strand, CG for genes on reverse strand
			ENS_ID+			"ANG",
			ENS_ID+			"APM",
			
			"WO",		// W09B6.4, Y59E9AR.9, ZK994.1, ZK994.3, ZK994.4, ... (no rule)
			"RDN"		// RDN58-1, snRNA, 15S-RNA, YAL016W, ... (no plan)
			
			// species missing
	};
	
	public static final boolean isValidEnsemblPrefix(String pfx) {
		for (int i = 0; i < SP_NAMES_ENS_PFX.length; i++) 
			if (pfx.equalsIgnoreCase(SP_NAMES_ENS_PFX[i]))
				return true;
		return false;
	}
	
	public static final String[] SP_NAMES_BINOMIAL= new String[] {
			"homo_sapiens", 	// mammals
			"pan_troglodytes",
			"mus_musculus", 
			"rattus_norvegicus",
			"canis_familiaris",
			"bos_taurus",
			"monodelphis_domestica",
			
			"gallus_gallus",	// vertebrates
			"xenopus_tropicalis",
			"danio_rerio",
			"takifugu_rubripes",
			"tetraodon_nigroviridis",
			
			"drosophila_melanogaster",	// insects
			"anopheles_gambiae",
			"apis_mellifera",	// not in ensembl
			
			"caenorhabditis_elegans",	// deuterostome?
			"saccharomyces_cerevisiae",
			
			"ciona_intestinalis",	// chordata?
			"ciona_savignyi",
			"ornithorhynchus_anatinus",	// platypus, vertebrate?
			
			"arabidopsis_thaliana"
	};
	
	public static final String[] SP_NAMES_COMMON = new String[] { 
		"human", 
		"chimp", 
		"mouse", 
		"rat", 
		"dog", 
		"cow", 
		"opossum", 
		
		"chicken", "frog", "zebrafish", "fugu", "tetraodon", 
		
		"fruitfly", "mosquito", "honeybee", "worm", 
		
		"yeast", "seasquirt", "seqsquirt2", "platypus", "cress" };
	public static final String[][] SPECIFIC_ANNOTATIONS = new String[][] { 
		new String[] {"Known", "RefSeq", "EnsEmbl", "MGC"}, 
		new String[] {"chimp"}, 
		new String[] {"Known", "RefSeq", "EnsEmbl", "MGC"}, 
		new String[] {"rat"}, 
		new String[] {"dog"}, 
		new String[] {"cow"}, 
		new String[] {"opossum"}, 
		new String[] {"chicken"}, 
		new String[] {"frog"}, 
		new String[] {"zebrafish"}, 
		new String[] {"fugu"}, 
		new String[] {"tetraodon"}, 
		new String[] {"RefSeq", "FlyBase"}, 
		new String[] {"mosquito"}, 
		new String[] {"honeybee"}, 
		new String[] {"RefSeq", "WormBase"}, 
		new String[] {"yeast"}, 
		new String[] {"seasquirt"}, 
		new String[] {"seqsquirt2"}, 
		new String[] {"platypus"}, 
		new String[] {"cress"} 
	};
	public static final String[] SP_NAMES_ANDRE = new String[] { "human", "mouse", "rat", "dog", "cow", "chicken", "tetraodon"};
	// public static final String[] SP_NAMES_ABBREV= new String[] {"hsapiens", "mmusculus", "rnorvegicus"};
	public static final String[] SP_NAMES_METAZOA= new String[] {
			"human",	// mammals 
			"chimp",
			"mouse", 
			"rat",
			"dog",
			"cow",

			"chicken",	// other chordates
			"frog",
			"zebrafish",
			
			"fruitfly",
			"honeybee",		// not in ensembl
			"worm"

	};

	public static Species[] createByBinomialName(String[] binomialNames) {
		
		if (binomialNames== null)
			return null;
		
		Species[] specs= new Species[binomialNames.length];
		for (int i = 0; i < specs.length; i++) 
			specs[i]= new Species(binomialNames[i]);
		
		return specs;
	}
	
	public int spNumber= -1;

	
	public static HashMap basePairMapping= new HashMap(5);
	static {
		basePairMapping.put(new Character('G'), new Character('C'));
		basePairMapping.put(new Character('C'), new Character('G'));
		basePairMapping.put(new Character('A'), new Character('T'));
		basePairMapping.put(new Character('T'), new Character('A'));
		basePairMapping.put(new Character('N'), new Character('N'));
	}

	/**
	 * @return Returns the ensemblPrefix.
	 */
	public String getEnsemblPrefix() {
		return SP_NAMES_ENS_PFX[spNumber];
	}

	public String getBinomialName() {
		return SP_NAMES_BINOMIAL[spNumber];
	}

	public String getAbbreviatedName() {
		String s= SP_NAMES_BINOMIAL[spNumber];
		int p= s.indexOf("_");
		s= Character.toUpperCase(s.charAt(0))+ "."+ s.substring(p+1, s.length());
		return s;
	}
	
	public String getCommonName() {
		return SP_NAMES_COMMON[spNumber];
	}
	
	public String getNameFirst() {
		return SP_NAMES_BINOMIAL[spNumber].substring(0, SP_NAMES_BINOMIAL[spNumber].indexOf('_'));
	}	
	
	public static final String getAbbrevNameForBinomial(String binName) {
		return binName.charAt(0)+ binName.substring(binName.indexOf("_")+ 1);		
	}
	
	public String getNameAbbrev() {
		return getAbbrevNameForBinomial(SP_NAMES_BINOMIAL[spNumber]);
	}
	
	

	HashMap geneHash= null;

	public static final int getSpeNrForPrefix(String ensemblPFX) {
		
		ensemblPFX= ensemblPFX.toUpperCase();
		
		for (int i = SP_NAMES_ENS_PFX.length- 1; i >= 0 ; --i)	// iterate backw, human matches everything 
			if (ensemblPFX.startsWith(SP_NAMES_ENS_PFX[i]))
				return i;
		
		return -1;
	}
	
	public static final String getCommonNameForPrefix(String ensemblPFX) {
		
		int speNr= getSpeNrForPrefix(ensemblPFX);
		if (speNr>= 0)
			return SP_NAMES_COMMON[speNr];
		
		return null;
	}
	
	public static final String getGenomeVersionForAnnotation(String build) {
		for (int i = 0; i < SP_ANNOTATION_VERSIONS.length; i++) 
			for (int j = 0; j < SP_ANNOTATION_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_ANNOTATION_VERSIONS[i][j].length; k++) 
					if (SP_ANNOTATION_VERSIONS[i][j][k].equalsIgnoreCase(build))
						return SP_GENOME_VERSIONS[i][j][0];
		return null;
	}
	
	public static final int getGenomeVerNbForAnnotation(String build) {
		for (int i = 0; i < SP_ANNOTATION_VERSIONS.length; i++) 
			for (int j = 0; j < SP_ANNOTATION_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_ANNOTATION_VERSIONS[i][j].length; k++) 
					if (SP_ANNOTATION_VERSIONS[i][j][k].equalsIgnoreCase(build))
						return j;
		return -1;
	}
	
	public static final int getGenomeVerNb(String genomeVersion) {
		for (int i = 0; i < SP_GENOME_VERSIONS.length; i++) 
			for (int j = 0; j < SP_GENOME_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_GENOME_VERSIONS[i][j].length; k++) 
					if (SP_GENOME_VERSIONS[i][j][k].equalsIgnoreCase(genomeVersion))
						return j;
		return -1;
	}

	public static final int getSpeciesNumberForAnnotation(String build) {
		for (int i = 0; i < SP_ANNOTATION_VERSIONS.length; i++) 
			for (int j = 0; j < SP_ANNOTATION_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_ANNOTATION_VERSIONS[i][j].length; k++) 
					if (SP_ANNOTATION_VERSIONS[i][j][k].equalsIgnoreCase(build))
						return i;
		return -1;
	}
	
	public static final int getSpeciesNumberForGenomeVersion(String genomeVer) {
		for (int i = 0; i < SP_GENOME_VERSIONS.length; i++) 
			for (int j = 0; j < SP_GENOME_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_GENOME_VERSIONS[i][j].length; k++) 
					if (SP_GENOME_VERSIONS[i][j][k].equalsIgnoreCase(genomeVer))
						return i;
		return -1;
	}

	public static final String getCommonNameForGenomeVersion(String build) {
		for (int i = 0; i < SP_GENOME_VERSIONS.length; i++) 
			for (int j = 0; j < SP_GENOME_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_GENOME_VERSIONS[i][j].length; k++) 
					if (SP_GENOME_VERSIONS[i][j][k].equalsIgnoreCase(build))
						return SP_NAMES_COMMON[i];
		return null;
	}
	
	public static final String decodeEnsemblPfx(String ensemblID) {
		Pattern patti= Pattern.compile("^(\\D{4,})\\d{11}.*");
		Matcher m= patti.matcher(ensemblID);
		if (m.matches())
			return m.group(1);
		return null;
	}
	
	public static final String getBinomialForCommonName(String commonName) {
		
		for (int i = 0; i < SP_NAMES_COMMON.length; i++) 
			if (SP_NAMES_COMMON[i].equalsIgnoreCase(commonName))
				return SP_NAMES_BINOMIAL[i];
		
		return null;
	}
	
	public static final String getBinomialForFirstName(String firstName) {
		
		firstName= firstName.toLowerCase();
		for (int i = 0; i < SP_NAMES_BINOMIAL.length; i++) 
			if (SP_NAMES_BINOMIAL[i].substring(0,SP_NAMES_BINOMIAL[i].indexOf('_')).equals(firstName))
				return SP_NAMES_BINOMIAL[i];
		
		return null;
	}	
	
	public static final String getBinomialForEnsemblPfx(String ensemblPfx) {
		
		for (int i = 0; i < SP_NAMES_ENS_PFX.length; i++) 
			if (SP_NAMES_ENS_PFX[i].equalsIgnoreCase(ensemblPfx))
				return SP_NAMES_BINOMIAL[i];
		
		return null;
	}		
	
	public static final String getBinomialForSomeName(String someName) {
		
		String bin= getBinomialForCommonName(someName);
		if (bin!= null)
			return bin;
		
		bin= getBinomialForFirstName(someName);
		if (bin!= null)
			return bin;
		
		bin= getBinomialForEnsemblPfx(someName);
		if (bin!= null)
			return bin;

		return null;
	}	
	
	public static final String getAbbrevNameForPrefix(String ensemblPFX) {
		
		for (int i = 0; i < SP_NAMES_ENS_PFX.length; i++) {
			if (SP_NAMES_ENS_PFX[i].equalsIgnoreCase(ensemblPFX))
				return getAbbrevNameForBinomial(SP_NAMES_BINOMIAL[i]);
		}
		
		return null;
	}
	
	public static int getSpeciesNumber(String someName) {
		
		return getSpeciesNumberForBinomialName(getBinomialForSomeName(someName));
	}
	
	public static int getAnnotationNumber(String someName, String annoName) {
		String[] annotations= SPECIFIC_ANNOTATIONS[getSpeciesNumber(someName)];
		for (int i = 0; i < annotations.length; i++) {
			if (annotations[i].equals(annoName))
				return i;
		}
		return -1;
	}
	
	public static int getSpeciesNumberForBinomialName(String binName) {
		int i;
		for (i = 0; i < SP_NAMES_BINOMIAL.length; i++) 
			if (SP_NAMES_BINOMIAL[i].equals(binName))
				break;
		if (i>= SP_NAMES_BINOMIAL.length) {
			//System.err.println("Unknown name "+ binName);
			return -1;
		}
		
		return i;
	}
	
	public Species(String someName) {

		int i= getSpeciesNumber(someName); 
			// else
		spNumber= i;
	}
	
	public Species(String someName, String genomeVer) {
		this(someName);
		setGenomeVersion(genomeVer);
	}
	

	
	public Gene[] getASGenes() {
		
		if (geneHash== null)
			return null;
	
		Iterator iter= geneHash.values().iterator();
		Gene g;
		Vector o= new Vector();
		while (iter.hasNext()) {
			g= (Gene) iter.next();
			if (g.getTranscripts().length> 1)
				o.add(g);
		}

		Gene[] result= new Gene[o.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Gene) o.elementAt(i);
		return result;	
	}
	
	
	public static char[] invertSequence(char[] seq) {
		
		char[] invertedSeq= new char[seq.length];
		for (int i = 0; i < invertedSeq.length; i++) {
			
			char c= seq[i];
			boolean lowerCase= false;
			if (Character.isLowerCase(c)) {
				c= Character.toUpperCase(c);
				lowerCase= true;
			}
			invertedSeq[i]= ((Character) basePairMapping.get(new Character(c))).charValue();
			if (lowerCase)
				invertedSeq[i]= Character.toLowerCase(invertedSeq[i]);
		}
		
		return invertedSeq;
	}
	
	/**
	 * @deprecated no longer used, see Graph.getSequenceDirectory
	 * @return the species directory
	 */
	public String getSequenceDirectory() {
	
		String seqDirName= Character.toUpperCase(getBinomialName().charAt(0))+ "."
				+ getBinomialName().substring(getBinomialName().indexOf('_')+ 1);
		File speciesGenome= new File(Constants2.DATA_DIR+ File.separator
				+ Constants2.SEQUENCES_SUBDIR+ File.separator
				+ seqDirName);
		Pattern patty= Pattern.compile("(\\D+)(\\d+)(\\D*)");
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
	
	public Species(String newBinomialName, int nbGenes) {
		
		this(newBinomialName);
		setEstimatedGeneNb(nbGenes);
	}
	
	/**
	 * Initilizes the <code>HashMap</code> with the specified capacity and the
	 * specified loading factor for the expected number of genes. 
	 * 
	 * @param nbGenes number of genes
	 * @param loadFactor load factor
	 * @return <code>true</code> when the <code>HashMap</code> was successfully
	 * initialized, <code>false</code> when it had already been filled with
	 * values beforeahead.
	 */
	public boolean setEstimatedGeneNb(int nbGenes, int loadFactor) {
		
		if (geneHash!= null&& geneHash.values().size()> 0)
			return false;
			
			// else
		geneHash= new HashMap(nbGenes, loadFactor);
		return true;
	}
	
	/**
	 * Initilizes the <code>HashMap</code> with the capacity and a default
	 * loading factor of <code>nbGenes/100</code> for the expected number 
	 * of genes.
	 * 
	 * @param nbGenes number of genes
	 * @return <code>true</code> when the <code>HashMap</code> was successfully
	 * initialized, <code>false</code> when it had already been filled with
	 * values beforeahead.
	 */
	public boolean setEstimatedGeneNb(int nbGenes) {
	
		return setEstimatedGeneNb(nbGenes, nbGenes/100);
	}

	
	public Gene getGene(String stableID) {

		return (Gene) geneHash.get(stableID);		// get gene
	}	
	public Gene[] getGenes() {
	
		if (geneHash== null)
			return null;
	
		Object[] o= geneHash.values().toArray();
		return Gene.toGeneArray(o);	
	}
	
	public String getDefaultGenomeVersion() {
		return Species.getDefaultGenomeVersion(getSpNumber());
	}
	
	public static String getDefaultGenomeVersion(int speNb) {
		return SP_GENOME_VERSIONS[speNb][0][0];
	}
	
	
	public Iterator getGeneIterator() {
		return geneHash.values().iterator();
	}

	public String[] getChromosomes() {
		Gene[] ge= getGenes();
		HashMap chrMap= new HashMap();
		for (int i = 0; i < ge.length; i++) {
			Object o= chrMap.get(ge[i].getChromosome());
			if (o== null)
				chrMap.put(ge[i].getChromosome(), ge[i].getChromosome());
		}
		
		return (String[]) ArrayUtils.toField(chrMap.keySet());
	}

	public boolean addGene(Gene newGene) {
		
			// get speciesHash
		if (geneHash== null) 
			geneHash= new HashMap();
		
			// 
		if (geneHash.get(newGene.getStableID())!= null)
			return false;	// gene exists
		geneHash.put(newGene.getStableID(), newGene);	// else
		return true;
	}

	public int getTranscriptCount() {
		if (geneHash== null) 
			return 0;
		Gene[] ge= getGenes();
		int cnt= 0;
		for (int i = 0; i < ge.length; i++) 
			cnt+= ge[i].getTranscriptCount();
		return cnt;
	}
	public Transcript[] getTranscripts(String speciesID) {
		
		Gene[] genes= getGenes();
		Vector transV= new Vector(genes.length);
		int nctr= 0;
		for (int i = 0; i < genes.length; i++) {			// for all genes
			Transcript[] trans= genes[i].getTranscripts();
			if (trans== null) {
				nctr++;
				continue;
			}
			for (int j = 0; j < trans.length; j++) 			// add all transcripts
				transV.add(trans[j]);
		}
		if (nctr> 0)
			System.err.println(nctr+ " genes wo transcipts!");
		
		Transcript[] result= new Transcript[transV.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Transcript) transV.elementAt(i);
		return result;
	}

	static final long serialVersionUID = 5861506554499528191L;

	/**
	 * More effective than <code>getGenes().length</code> since it does not 
	 * convert <code>Object[]</code> to <code>Gene[]</code>.
	 * @return the number of genes
	 */
	public int getGeneNb() {
	
		if (geneHash== null|| geneHash.values()== null)
			return 0;
		
		return geneHash.values().size();
	}
	
	/**
	 * @deprecated deactivated
	 * @param regionType the type of region of interest
	 * @param ssType the type of (splice) site of interest
	 * @return a vector of (splice) sites
	 */
	public SpliceSite[] getSpliceSites(int regionType, int ssType) {
//		int perc= 0;
//		Gene[] ge= getGenes();
//		Vector v= new Vector();
//		for (int i = 0; i < ge.length; i++) {
//			SpliceSite[] ssites= ge[i].getSpliceSites(ssType, regionType);
//			for (int j = 0; ssites!= null&& j < ssites.length; j++) 
//				v.add(ssites[j]);
//		}
//		return (SpliceSite[]) ArrayUtils.toField(v);
		return null;
	}
	/**
	 * @return Returns the buildVersion.
	 */
	public String getAnnotationVersion() {
		return annotationVersion;
	}
	public Exon[] getExons() {

		Vector v= new Vector();
		Iterator iter= geneHash.values().iterator();
		while(iter.hasNext()) {
			Gene ge= (Gene) iter.next();
			Exon[] ex= ge.getExons();
			for (int i = 0; i < ex.length; i++) {
				v.add(ex[i]);
			}
		}
		return (Exon[]) ArrayUtils.toField(v);
	}

	/**
	 * @param buildVersion The buildVersion to set.
	 */
	public void setAnnotationVersion(String buildVersion) {
		this.annotationVersion = buildVersion;
	}

	public String toString() {
		return getCommonName()+" "+getGenomeVersion();
	}
	public int countNonProteinCodingLoci() {
		Gene[] ge= getGenes();
		int cnt= 0;
		for (int i = 0; i < ge.length; i++) 
			if (!ge[i].isProteinCoding())
				++cnt;
		return cnt;
	}
	
	/**
	 * @deprecated deactivated
	 */
	public void initTU() {
//		Gene[] ge= getGenes();
//		for (int i = 0; i < ge.length; i++) {
//			ge[i].initTU();
//		}
	}
	
	public int getSpNumber() {
		return spNumber;
	}

	String annotationVersion = null;

	public String getGenomeVersion() {
		return genomeVersion;
	}

	public void setGenomeVersion(String genomeVersion) {
		this.genomeVersion = genomeVersion;
	}

	public static final String getAnnotation(String speName, String genomeVer, String annoName, String[] keywords) {
		for (int i = 0; keywords!= null&& i < keywords.length; i++) 
			keywords[i]= keywords[i].toUpperCase();
		
		speName= speName.toUpperCase();
		annoName= annoName.toUpperCase();
		String[] list= new File(Constants2.SUBDIR_ANNOTATION).list();		
		Vector v= new Vector();
		for (int i = 0; i < list.length; i++) {
			String s= list[i].toUpperCase();
			if (s.contains(speName)&& s.contains(annoName)) {
				if (genomeVer!= null) {
					genomeVer= genomeVer.toUpperCase();
					if (!s.contains(genomeVer))
						continue;
				}
				if (keywords!= null) {
					int x;
					for (x = 0; x < keywords.length; x++) 
						if (s.contains(keywords[x]))
							break;
					if (x< keywords.length)
						v.add(new File(Constants2.SUBDIR_ANNOTATION+ File.separator+ list[i]).getAbsolutePath());
				} else if (!new barna.model.commons.MyFile(list[i]).getExtension().contains("_"))
					v.add(new File(Constants2.SUBDIR_ANNOTATION+ File.separator+ list[i]).getAbsolutePath());
			}
		}
		
		if (v.size()== 0) {
			System.err.println("No annotation found for "+speName+", "+annoName);
			return null;
		}
		if (v.size()== 1) 
			return (String) v.elementAt(0);
	
		System.out.println("Ambiguous information, press <CR> for fetching newest relase (according to 200x date).");
		try {
			System.in.read();
		} catch (IOException e) {
			e.printStackTrace();
		}
		int max= -1;
		int maxP= -1;
		for (int i = 0; i < v.size(); i++) {
			String s= (String) v.elementAt(i);
			int p= s.indexOf("200");// 200x
			if (p< 0) {
				System.err.println("No date available for "+s);
				continue;
			}
			int x= Integer.parseInt(s.substring(p+3, p+6));
			if (x> max) {
				max= x;
				maxP= i;
			}
		}
		
		return (String) v.elementAt(maxP);
	}
}
