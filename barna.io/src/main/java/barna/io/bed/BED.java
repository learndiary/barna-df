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

package barna.io.bed;

import barna.model.DirectedRegion;
import barna.model.SpliceSite;
import barna.model.Transcript;
import barna.model.bed.BEDobject;
import barna.model.constants.Pedro;

import java.awt.*;
import java.util.HashMap;
import java.util.Vector;

/**
 * BED static utilities
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class BED {

    /**
     * @deprecated change in BEDobject coordinate handling
     * @param digSite
     * @param tagLen
     * @param nb
     * @param colorMap
     * @return
     */
    public static BEDobject[] extract3TagRegions(Transcript transcript, String digSite, int tagLen, int nb, HashMap<String, Color> colorMap) {
        String seq= transcript.getSplicedSequence();
        Vector<BEDobject> v= new Vector<BEDobject>(nb);
        int cnt= 0;

        for (int k = seq.length()-4; k >0; --k) {	// KMP?
            if (seq.substring(k, k+4).equalsIgnoreCase(digSite)) {
                int end= Math.min(k+4+tagLen, seq.length());
                String tag= seq.substring(k+4, end);
                while (tag.length()< tagLen)
                    tag+= "X";	// fill with unknown


                DirectedRegion reg= new DirectedRegion(
                        transcript.getGenomicPosition(k),	// without digestion site: +4
                        transcript.getGenomicPosition(end - 1),
                        transcript.getStrand()
                );
                reg.setChromosome(transcript.getChromosome());
                DirectedRegion[] regs= DirectedRegion.intersect(transcript.getExons(), new DirectedRegion[] {reg});
                if (transcript.getStrand()< 0&& regs.length> 1) {
                    DirectedRegion[] regsRev= new DirectedRegion[regs.length];
                    for (int i = 0; i < regs.length; i++)
                        regsRev[regsRev.length- 1- i]= regs[i];
                    regs= regsRev;
                }

                BEDobject obj= new BEDobject();
                obj.setName(transcript.getTranscriptID());
                obj.setChrom(transcript.getChromosome());
                obj.setStrand(transcript.getStrand());
                obj.setBlockCount(regs.length);
                int st= regs[0].getStart();
                int nd= regs[regs.length-1].getEnd();
                obj.setStart(st);
                obj.setEnd(nd);
                obj.setStart(obj.getStart()-1); // bed shit
                //obj.setEnd(obj.getEnd()+1);
                int[] starts= new int[regs.length], sizes= new int[regs.length];
                for (int i = 0; i < sizes.length; i++) {
                    starts[i]= Math.abs(regs[i].getStart()- st);
                    sizes[i]= regs[i].get3PrimeEdge()- regs[i].get5PrimeEdge()+ 1;	// maybe + 1 !!, bed shit..
                }
                obj.setBlockStarts(BEDobject.toString(starts));
                obj.setBlockSizes(BEDobject.toString(sizes));


                String tagName= Pedro.getDigSiteQualifier(digSite)+"_"+(++cnt);
                if (colorMap!= null&& colorMap.get(tagName)!= null)
                    obj.setCol(colorMap.get(tagName));

                boolean ovlEJ= false; SpliceSite ssX= null, ssY= null;
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
//				obj.setName(getTranscriptID()+"_"+tagName+
//						"("+(exons.length- j)+".last:"+(seq.length()- end)+"nt)");
                obj.setName(tagName);
//				obj.addAttribute(GFFObject.TRANSCRIPT_ID_TAG, getTranscriptID());
//				obj.addAttribute(GTF_ATTRIBUTE_TAG_LAST_EXON, Integer.toString(exons.length- j));
//				obj.addAttribute(GTF_ATTRIBUTE_TAG_ENDDIST, Integer.toString(seq.length()- end));

                v.add(obj);
                if (v.size()== nb)
                    break;
            }
        }

        BEDobject[] beds= new BEDobject[v.size()];
        for (int i = 0; i < beds.length; i++)
            beds[i]= v.elementAt(i);
        return beds;
    }

}
