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

package fbi.genome.model;

import fbi.commons.StringUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

public class CDSEvent extends ASEvent {

	Transcript[][][] trpts;
	SpliceSite[][][] spliceChains; 
	
	public CDSEvent(Transcript[][][] newTrpts, SpliceSite[][][] ss) {
		super();
		this.trpts= newTrpts;
		this.spliceChains= ss;
	}

	/**
		 * new GTF compliance
		 */
		public String toStringGTF(boolean outSeq) {
		
				// build final string
			StringBuffer result= new StringBuffer(getGene().getChromosome());
			result.append('\t');
	
			result.append(getGTFSource());
			result.append('\t');
	
			byte type= getType();
			if (type== TYPE_AS_EVENT)
				result.append(GTF_FEATURE_ASEVENT);
			else if (type== TYPE_DS_EVENT)
				result.append(GTF_FEATURE_DSEVENT);
			else if (type== TYPE_VS_EVENT)
				result.append(GTF_FEATURE_VSEVENT);
			else 
				result.append(GTF_FEATURE_NOEVENT);
			result.append('\t');
	
	
			int strand= getGene().getStrand();
			DirectedRegion reg= getRegionEvent();
			int srcPos= reg.get5PrimeEdge();
			int snkPos= reg.get3PrimeEdge();
			if (strand> 0) {
				result.append(srcPos);
				result.append('\t');
				result.append(snkPos);
				result.append('\t');
			} else {
				result.append(Math.abs(snkPos));
				result.append('\t');
				result.append(Math.abs(srcPos));
				result.append('\t');
			}
			
				// score
			float score= -1f; 	// getScore(); 091208 SLOW !!! check overlap()
			if (score== -1f)
				result.append(".");
			else
				result.append(StringUtils.fprint(getScore(), 5));
			result.append('\t');
			
			if (strand> 0)
				result.append(ASEvent.STRAND_SYMBOL_POS);
			else
				result.append(ASEvent.STRAND_SYMBOL_NEG);
			result.append('\t');
			
			result.append('.');
			result.append('\t');
				
				// attributes
			result.append(ASEvent.TRANSCRIPT_ID_TAG);
			result.append(' ');
			result.append('\"');
			for (int i = 0; i < trpts.length; i++) {
				for (int j = 0; j < trpts[i].length; j++) {
					for (int k = 0; k < trpts[i][j].length; k++) {
						result.append(trpts[i][j][k].getTranscriptID());
						result.append('/');
					}
					result.replace(result.length()- 1, result.length(), ";");
				}
				result.replace(result.length()- 1, result.length(), ",");
			}
			result.deleteCharAt(result.length()- 1);
			result.append('\"');
			result.append(';');
			result.append(' ');
	
			result.append(ASEvent.GENE_ID_TAG);
			result.append(' ');
			result.append('\"');
			result.append(getGene().getGeneID());	// getGene().getReferenceTranscript()
			result.append('\"');
			result.append(';');
			result.append(' ');
	
			result.append(GTF_ATTRIBUTE_TAG_FLANKS);
			result.append(' ');
			result.append('\"');
			if (getSrc().getPos()== Integer.MIN_VALUE)
				result.append("null");
			else {
				result.append(Math.abs(getSrc().getPos()));
				result.append(getSrc().getSiteSymbol());
			}
			result.append(",");
			if (getSnk().getPos()== Integer.MAX_VALUE)
				result.append("null");
			else {
				result.append(Math.abs(getSnk().getPos()));
				result.append(getSnk().getSiteSymbol());
			}
			result.append('\"');
			result.append(';');
			result.append(' ');
	
			
			result.append(GTF_ATTRIBUTE_TAG_STRUCTURE);
			result.append(' ');
			result.append('\"');
			String evCode= toString();
			result.append(evCode);
			result.append('\"');
			result.append(';');
			result.append(' ');
	
			result.append(GTF_ATTRIBUTE_TAG_SPLICE_CHAIN);	// 5' -> 3'
			result.append(' ');
			result.append('\"');
			for (int i = 0; i < spliceChains.length; i++) {
				for (int j = 0; j < spliceChains[i].length; j++) {
					for (int k = 0; k < spliceChains[i][j].length; k++) {
						result.append(Math.abs(spliceChains[i][j][k].getPos()));
						result.append(spliceChains[i][j][k].getSiteSymbol()); //'/'
					}
					result.append(';');
				}
				result.append(',');
			}
			result.deleteCharAt(result.length()- 1);
			result.append('\"');
			result.append(';');
			result.append(' ');
	
			result.append(GTF_ATTRIBUTE_TAG_SOURCES);
			result.append(' ');
			result.append('\"');
			String[] sources= getGTFSources();
			for (int i = 0; i < sources.length; i++) {
				result.append(sources[i]);
				result.append(',');
			}
			result.deleteCharAt(result.length()-1);
			result.append('\"');
			result.append(';');
			result.append(' ');
	
			addAttribute(GTF_ATTRIBUTE_TAG_DEGREE, new Integer(getDegree()).toString()); 
			if (getAttribute(GTF_ATTRIBUTE_TAG_DIMENSION)== null)
				addAttribute(GTF_ATTRIBUTE_TAG_DIMENSION, new Integer(trpts.length).toString());
			
			if (outSeq) {
				SpliceSite[] su= getSpliceUniverse(true);
				StringBuffer sb= new StringBuffer(su.length* 25);
				for (int i = 0; i < su.length; i++) {
					sb.append(i);
					sb.append(su[i].getSiteSymbol());
					if (su[i].isTSS()|| su[i].isTES()) {
						sb.append(",");
						continue;
					}
					String s= null;
					if (su[i].isDonor())
					 s= Graph.readSequence(su[i], 15, 4); // us: 3exonic ss+ 4 codons
					else
					 s= Graph.readSequence(su[i], 15, 13);	// ds: 1exonic ss+ 4 codons
					sb.append(s);
					sb.append(",");
				}
				sb.deleteCharAt(sb.length()-1);
				addAttribute(GTF_ATTRIBUTE_SITE_SEQ, sb.toString());
			}
			
			if (outputFlankMode) {
				result.append(GTF_ATTRIBUTE_FLANK_MODE);
				result.append(" \"");
				if (getSrc().getPos()== Integer.MIN_VALUE)
					result.append("null");
				else {
					String mode= "constitutive";
					for (int i = 0; i < getGene().getTranscripts().length; i++) {
						Transcript t= getGene().getTranscripts()[i];
						if (t.get5PrimeEdge()<= getSrc().getPos()&&
								t.get3PrimeEdge()>= getSrc().getPos()) {
							int x;
							for (x = 0; x < getSrc().getTranscripts().length; x++) 
								if (getSrc().getTranscripts()[x]== t)
									break;
							if (x== getSrc().getTranscripts().length) {
								mode= "alternative";
								break;
							}
						}
					}
					result.append(mode);
				}
				result.append(",");
				if (getSnk().getPos()== Integer.MAX_VALUE)
					result.append("null");
				else {
					String mode= "constitutive";
					for (int i = 0; i < getGene().getTranscripts().length; i++) {
						Transcript t= getGene().getTranscripts()[i];
						if (t.get5PrimeEdge()<= getSnk().getPos()&&
								t.get3PrimeEdge()>= getSnk().getPos()) {
							int x;
							for (x = 0; x < getSnk().getTranscripts().length; x++) 
								if (getSnk().getTranscripts()[x]== t)
									break;
							if (x== getSnk().getTranscripts().length) {
								mode= "alternative";
								break;
							}
						}
					}
					result.append(mode);
				}
				result.append("\"; ");
			}
			
			if (cdsLoc!= null) {
				result.append(GTF_ATTRIBUTE_CDS_LOC);
				result.append(" \"");
				for (int j = 0; j < cdsLoc.length; j++) {
					if (cdsLoc[j]==null)
						result.append("-,");
					else
						result.append(cdsLoc[j]+",");
				}
				result.deleteCharAt(result.length()- 1);
				result.append("\"; ");
				
				result.append(GTF_ATTRIBUTE_CDS_IMPACT);
				result.append(" \"");
				HashMap<String,String> littleMap= new HashMap<String, String>(8);
				for (int i = 0; i < cdsLoc.length; i++) 
					littleMap.put(cdsLoc[i], null);
				Object[] oo= littleMap.keySet().toArray();
				Arrays.sort(oo);
				for (int i = 0; i < oo.length; i++) 
					result.append(oo[i]+",");
				result.deleteCharAt(result.length()- 1);
				result.append("\"; ");
			}
			
	
			Iterator<String> iter= getAttributeMap().keySet().iterator();
			while (iter.hasNext()) {
				String key= iter.next();
				result.append(key);
				result.append(' ');
				result.append('\"');
				result.append(getAttribute(key));
				result.append('\"');
				result.append(';');
				result.append(' ');
			}
			
			return result.toString();
		}

}
