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

package barna.genome.io.gtf;

import barna.genome.model.*;
import barna.genome.model.constants.Pedro;
import barna.genome.model.gff.GFFObject;

import java.util.Iterator;
import java.util.Vector;

/**
 * GFF utilities to create the appropriate model classes from internal GFF data structures
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class GFF {

    public static ASEvent createASEvent(GFFObject obj){
        ASEvent event = new ASEvent();
        Gene gene= new Gene(obj.getAttribute(GFFObject.GENE_ID_TAG));
            gene.setChromosome(obj.getSeqname());
            gene.setStrand((byte) obj.getStrand());

            if (obj.getFeature().equals(ASEvent.GTF_FEATURE_ASEVENT))
                event.setType(ASEvent.TYPE_AS_EVENT);
            else if(obj.getFeature().equals(ASEvent.GTF_FEATURE_DSEVENT))
                event.setType(ASEvent.TYPE_DS_EVENT);
            else if(obj.getFeature().equals(ASEvent.GTF_FEATURE_VSEVENT))
                event.setType(ASEvent.TYPE_VS_EVENT);

                // flanks, convert null to infinity nodes
            SpliceSite src, snk;
            if (obj.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_FLANKS)!= null) {
                String[] flanks= obj.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_FLANKS).split(",");
                if (flanks[0].equals("null"))
                    src= new SpliceSite(Integer.MIN_VALUE, SpliceSite.TYPE_NOT_INITED, gene);
                else {
                    int pos= Integer.parseInt(flanks[0].substring(0,flanks[0].length()-1));
                    byte type= SpliceSite.getTypeBySymbol(flanks[0].charAt(flanks[0].length()-1));
                    src= new SpliceSite(pos, type, gene);
                }
                if (flanks[1].equals("null"))
                    snk= new SpliceSite(Integer.MAX_VALUE, SpliceSite.TYPE_NOT_INITED, gene);
                else {
                    int pos= Integer.parseInt(flanks[1].substring(0,flanks[1].length()-1));
                    byte type= SpliceSite.getTypeBySymbol(flanks[1].charAt(flanks[1].length()-1));
                    snk= new SpliceSite(pos, type, gene);
                }
            } else {
                int srcPos= obj.getStart();
                int snkPos= obj.getEnd();
                if (obj.getStrand()< 0) {
                    srcPos= -obj.getEnd();
                    snkPos= -obj.getStart();
                }
                src= new SpliceSite(srcPos, SpliceSite.TYPE_NOT_INITED, gene);	// init type after schains
                snk= new SpliceSite(snkPos, SpliceSite.TYPE_NOT_INITED, gene);
            }
            event.setSrc(src);
            event.setSnk(snk);

                // tids, sources..
            String[] tokens= obj.getAttribute(GFFObject.TRANSCRIPT_ID_TAG).split(",");
            String[] tokenSource= obj.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_SOURCES).split(",");
            Transcript[][] trpts= new Transcript [tokens.length][];
            for (int i = 0; i < tokens.length; i++) {
                String[] tokens2= tokens[i].split("/");
                trpts[i]= new Transcript[tokens2.length];
                for (int j = 0; j < tokens2.length; j++) {
                    trpts[i][j]= new Transcript(gene, tokens2[j]);
                    trpts[i][j].setSource(tokenSource[i]);
                }
            }
            event.setTranscripts(trpts);

            SpliceSite[][] spliceChains= new SpliceSite[trpts.length][];
            tokens= obj.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_SPLICE_CHAIN).split(",");
            String[] tokStruct= obj.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_STRUCTURE).split(",");
            for (int i = 0; i < tokens.length; i++) {
                if (tokens[i].length()== 0) {
                    spliceChains[i]= new SpliceSite[0];
                } else {
                    String[] tokens2= tokens[i].split("\\D");	// 080822: changed from "/"
                    String[] tokStruct2= tokStruct[i].split("\\d");
                    spliceChains[i]= new SpliceSite[tokens2.length];
                    for (int j = 0, ts2Pos= 0; j < tokens2.length; ++j, ++ts2Pos) {
                        int pos= Integer.parseInt(tokens2[j]);
                        if (gene.getStrand()< 0)
                            pos= -pos;
                        while (tokStruct2[ts2Pos].length()== 0)
                            ++ts2Pos;
                        byte type= SpliceSite.getTypeBySymbol(tokStruct2[ts2Pos].charAt(0));
                        spliceChains[i][j]= new SpliceSite(pos, type, gene);
                    }
                }
            }
            event.setSpliceChains(spliceChains);	// sort

            if (obj.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_FLANKS)== null) {
                SpliceSite[] su= event.getSpliceUniverse(false);
                if (su[0].isLeftFlank())
                    src.setType(SpliceSite.TYPE_DONOR);
                else
                    src.setType(SpliceSite.TYPE_ACCEPTOR);
                if (su[su.length-1].isLeftFlank())
                    snk.setType(SpliceSite.TYPE_DONOR);
                else
                    snk.setType(SpliceSite.TYPE_ACCEPTOR);
            }

            Iterator<String> iter= obj.getAttributes().keySet().iterator();
            while(iter.hasNext()) {
                String key= iter.next();
                if (key.equals(ASEvent.GTF_ATTRIBUTE_TAG_SPLICE_CHAIN)|| key.equals(GFFObject.TRANSCRIPT_ID_TAG))	// not save redundantly
                    // || key.equals(GTF_ATTRIBUTE_TAG_STRUCTURE)  -->  needed in LaVista for flexible group tag init
                    continue;
                event.addAttribute(key, obj.getAttribute(key));
            }
        return event;
    }



    /**
     * Create a DirectedRegion from a given GFFObject
     *
     * @param obj the source
     * @return region the created region
     */
    public static DirectedRegion createDirectedRegion(GFFObject obj) {
        DirectedRegion region = new DirectedRegion();
        region.setStrand((byte) obj.getStrand());
        region.setStart(obj.getStart());
        region.setEnd(obj.getEnd());
        region.setChromosome(obj.getChromosome());
        region.setID(obj.getFeature());

            // clone?!
        if (obj.getAttributes()!= null)
            region.setAttributes(obj.getAttributes());
        return region;
    }


    public static GFFObject[] fromTranscript(Transcript site) {
        return fromTranscript(site, true, true, false, false, false);
    }

    public static GFFObject[] fromTranscript(Transcript transcript, boolean outputExons, boolean outputCDS, boolean outputStartCodon, boolean outputStopCodon, boolean outputSpliceSites) {

        Vector<GFFObject> v= new Vector<GFFObject>();
        if (outputExons|| outputCDS) {
            Translation ta= null;
            if (transcript.getTranslations()!= null)
                ta= transcript.getTranslations()[0];
            for (int i = 0; i < transcript.getExons().length; i++) {
                GFFObject[] obs= GFFObject.createGTFObjects(transcript.getExons()[i], transcript);

                // (1) exon
                if (outputExons)
                    v.add(obs[0]);

                // (2) start codon
                if (ta!= null&& ta.getCodonStart()!= null
                        && ta.getCodonStart().getPos()>= transcript.getExons()[i].get5PrimeEdge()
                        && ta.getCodonStart().getPos()<= transcript.getExons()[i].get3PrimeEdge())
                    v.add(GFFObject.createGFFObject(ta.getCodonStart(), transcript));

                // (3) cds
                if (outputCDS&& obs.length> 1)
                    v.add(obs[1]);

                // (4) stop codon
                if (ta!= null&& ta.getCodonStop()!= null
                        && ta.getCodonStop().getPos()>= transcript.getExons()[i].get5PrimeEdge()
                        && ta.getCodonStop().getPos()<= transcript.getExons()[i].get3PrimeEdge())
                    v.add(GFFObject.createGFFObject(ta.getCodonStop(), transcript));
            }
        }

        GFFObject[] obs= new GFFObject[v.size()];
        for (int i = 0; i < obs.length; i++)
            obs[i]= v.elementAt(i);
        return obs;
    }


    public static GFFObject[] extract3Tag(Transcript transcript, String digSite, int tagLen, int nb) {
            String seq= transcript.getSplicedSequence();
            GFFObject[] obs= new GFFObject[nb];
            int cnt= 0;
            for (int k = seq.length()-4; k >0; --k) {	// KMP?
                if (seq.substring(k, k+4).equalsIgnoreCase(digSite)) {
                    int end= Math.min(k+4+tagLen, seq.length());
                    String tag= seq.substring(k+4, end);
                    while (tag.length()< tagLen)
                        tag+= "X";	// fill with unknown

                    GFFObject obj= new GFFObject();
                    obj.setStart(transcript.getGenomicPosition(k + 4));
                    obj.setEnd(transcript.getGenomicPosition(end - 1));
                    obj.setStrand(transcript.getStrand());
                    obj.setFeature(Pedro.getDigSiteQualifier(digSite)+"_"+(cnt+1));
                    obj.setSeqname(transcript.getChromosome());
                    obj.setSource(transcript.getSource());
                    //obj.addAttribute(GFFObject.GENE_ID_TAG, getGene().getGeneID());
                    obj.addAttribute(GFFObject.TRANSCRIPT_ID_TAG, transcript.getTranscriptID());
                    obj.addAttribute(GFFObject.ID_ATTRIBUTE_SEQUENCE, tag);

                    boolean ovlEJ= false; SpliceSite ssX= null, ssY= null;
                    DirectedRegion reg= new DirectedRegion(obj.getStart(), obj.getEnd(), obj.getStrand());
                    reg.setChromosome(obj.getChromosome());
                    int j;
                    for (j = transcript.getExons().length- 1; j > 0; --j) {
                        ssX= transcript.getExons()[j].getAcceptor();
                        ssY= transcript.getExons()[j-1].getDonor();
                        int x= ssX.getPos();
                        int y= ssY.getPos();
                        if (reg.contains(x)&& reg.contains(y))
                            ovlEJ= true;
                        if (reg.get5PrimeEdge()> transcript.getExons()[j-1].get3PrimeEdge())
                            break;
                    }
                    obj.addAttribute("endExonNr", new Integer(transcript.getExons().length- j).toString());
                    obj.addAttribute("endDist", new Integer(seq.length()- end).toString());

                    String loc= "UNDEFINED";
                    if (transcript.isCoding()) {
                        Translation tln= transcript.getTranslations()[0];
                        if (tln.overlaps(reg))
                            loc= "CDS";
                        else if (tln.get3PrimeEdge()< reg.get5PrimeEdge())
                            loc= "3UTR";
                        else if (tln.get5PrimeEdge()> reg.get3PrimeEdge())
                            loc= "5UTR";
                    }
                    obj.addAttribute("location", loc);

                    obj.addAttribute("ovlEJ", new Boolean(ovlEJ).toString());
                    if (ovlEJ) {
                        String mod= "constitutive";
                        if (ssX.isAlternative()|| ssY.isAlternative())
                            mod= "alternative";
                        obj.addAttribute("modality", mod);
    //					if (ssX.isAlternative()|| ssY.isAlternative()) {
    //						HashMap map= new HashMap();
    //						for (int i = 0; ssX.isAlternative()&& i < ssX.getAsVars().length; i++)
    //							if (map.get(ssX.getAsVars()[i])== null)
    //								map.put(ssX.getAsVars()[i], ssX.getAsVars()[i]);
    //						for (int i = 0; ssY.isAlternative()&& i < ssY.getAsVars().length; i++)
    //							if (map.get(ssY.getAsVars()[i])== null)
    //								map.put(ssY.getAsVars()[i], ssY.getAsVars()[i]);
    //
    //						StringBuffer ev= new StringBuffer();
    //						Object[] o= map.keySet().toArray();
    //						for (int i = 0; i < o.length; i++)
    //							ev.append(o.toString()+"|");
    //						ev.deleteCharAt(ev.length()- 1);
    //
    //						obj.addAttribute("events", ev.toString());
    //					}
                    }

                    transcript.addAttribute(digSite, tag);
                    obs[cnt++]= obj;
                    if (cnt== obs.length)
                        return obs;
                }
            }
            return obs;
        }



}
