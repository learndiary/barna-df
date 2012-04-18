/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.model.constants;


import barna.model.Species;
import barna.model.commons.MyFile;

import java.awt.*;
import java.util.HashMap;

public class Pedro {
	
	public static Color[] COLS_SAGE= new Color[] {
		new Color(151,27,30),	// mars red 
		new Color(192,39,45),	// ruby
		Color.red, 
		new Color(242, 101, 34)};	// pumpkin
	public static Color[] COLS_MPSS= new Color[] {
		new Color(0,84,166),	// starry night blue 
		new Color(0,113,188), // twilight blue
		new Color(0,143,213), // hawaian blue 
		Color.cyan};
	
	static final String ID_CNT_ALL= "cntAll";
	static final String ID_CNT_3INCOMPLETE= "cnt3Pincl";
	static final String ID_CNT_REAL= "cntReal";
	static final String ID_CNT_MULT= "cntMulti";
	static final String ID_CNT_SINGLE= "cntSingle";
	static final String ID_CNT_MULT_TRPT= "cntMultTrpt";
	static final String ID_CNT_MULT_TAG= "cntMultTag";
	static final String ID_CNT_MULT_TAG_SAGE= "cntMultTagSAGE";
	static final String ID_CNT_REAL_TRPT= "cntRealTrpt";
	static final String ID_CNT_REAL_SAGE= "cntRealSAGE";
	static final String ID_CNT_MULT_SAGE= "cntMultiSAGE";
	static final String ID_CNT_SINGLE_SAGE= "cntSingleSAGE";
	static final String ID_CNT_MULT_TRPT_SAGE= "cntMultTrptSAGE";
	static final String ID_CNT_REAL_TRPT_SAGE= "cntRealTrptSAGE";
	static final String ID_CNT_SGLE_TRUNC= "cntSglTrunc";
	static final String ID_CNT_MULTI_TRUNC= "cntMultiTrunc";
	static final String ID_CNT_SGLE_REPEAT= "cntSglRepeat";
	static final String ID_CNT_MULTI_REPEAT= "cntMultiRepeat";
	static final String ID_CNT_SGLE_TRUNC_SAGE= "cntSglTruncSAGE";
	static final String ID_CNT_MULTI_TRUNC_SAGE= "cntMultiTruncSAGE";
	static final String ID_CNT_SGLE_REPEAT_SAGE= "cntSglRepeatSAGE";
	static final String ID_CNT_MULTI_REPEAT_SAGE= "cntMultiRepeatSAGE";
	static final String ID_DIST_FOUND= "distFound";
	static final String ID_DIST_NOT_FOUND= "distNotFound";
	static final String ID_CNT_SGLE_TRPT= "cntSingleTrpt";
	static final String ID_CNT_SGLE_TRPT_SAGE= "cntSingleTrptSAGE";
	static final String ID_GENE_WO_TAG_SAGE= "cntGeneWOtagSAGE";
	static final String ID_GENE_WO_TAG_MPSS= "cntGeneWOtagMPSS";
	static final String ID_GENE_INCOMPLETE_MPSS= "cntGeneIncompMPSS";
	static final String ID_GENE_INCOMPLETE_SAGE= "cntGeneIncompSAGE";
	static final String ID_GENE_COMPLETE_MPSS= "cntGeneComplMPSS";
	static final String ID_GENE_COMPLETE_SAGE= "cntGeneComplSAGE";
	
	
	public final static String MPSS_DIGEST_SITE_DPNII= "GATC";
	public final static String MPSS_TAG_QUALIFIER= "mpss_tag";
	public final static int MPSS_TAG_LENGTH= 13;
	public final static String SAGE_DIGEST_SITE_NLAIII= "CATG";
	public final static int SAGE_TAG_LENGTH_LONG= 17;
	public final static int SAGE_TAG_LENGTH_SHORT= 10;
	public final static String SAGE_TAG_QUALIFIER= "sage_tag";
	public final static String ID_TAG_TYPE= "tag_type";
	
	
	public final static String ID_DIGEST= "tag_id"; 
	public final static String ID_TAGLEN= "tag_length"; 

	static HashMap realTagHash= null;
	public static boolean isMPSStag(int len) {
		if (len== MPSS_TAG_LENGTH)
			return true;
		return false;
	}

	void checkComplete3PSandro() {
		
		// check polyA:
			// at least 8nt A at the end
		
		// check internal priming: 
		// poly-T for mRNA generation hybridizes with A-rich region
		// in last exon
			// sliding window, size 10, 8nt A
		
	}
	
	
	barna.model.commons.MyFile file;
	Species species;
	boolean writeBED= true, writeGTF= true;
	
	public Pedro(MyFile newFile, String genomeName) {
		this.file= newFile;
		setGenomeVersion(genomeName);
	}
	public void setGenomeVersion(String genomeName) {
		String[] token= genomeName.split("_");
		this.species= new Species(token[0]);
		species.setGenomeVersion(token[1]);
	}

/*
	int checkEJCOverlap(GFFObject[] tags, Exon[] e, int taglen, int[] hits, int[] total) {
		int hitCount= 0;
		for (int i = 0; i < tags.length; i++) {
			if (tags[i]== null)
				continue;

			++total[i];
			
			// too unsafe
			//if (tags[i].getEnd()- tags[i].getStart()+ 1> taglen) {	// < can be truncated tag
			int j;
			for (j = e.length- 1; j > 0; --j) {
				int x= Math.abs(e[j].get5PrimeEdge());
				int y= Math.abs(e[j-1].get3PrimeEdge());
				if (x> tags[i].getStart()&& x< tags[i].getEnd()&& y> tags[i].getStart()&& y< tags[i].getEnd()) {
					++hitCount;
					++hits[i];
					break;
				}
			}
				
		}
		return hitCount;
	}
*/

	
	public static String getDigSiteQualifier(String digSite) {
		if (digSite.equalsIgnoreCase(MPSS_DIGEST_SITE_DPNII))
			return MPSS_TAG_QUALIFIER;
		if (digSite.equalsIgnoreCase(SAGE_DIGEST_SITE_NLAIII))
			return SAGE_TAG_QUALIFIER;
		return null;
	}
}


