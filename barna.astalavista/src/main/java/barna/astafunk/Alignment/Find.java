package barna.astafunk.Alignment;

import barna.commons.log.Log;
import barna.astafunk.Astafunk;
import barna.astafunk.DP.DyPMatrix;
import barna.astafunk.DP.Hit;
import barna.astafunk.HMM.ProfileHMM;
import barna.astafunk.utils.Utils;
import barna.model.DirectedRegion;
import barna.model.Transcript;
import barna.model.Translation;

import java.util.*;

/**
 * @version 2
 * @autor Vitor Coelho
 * @since 15/05/14
 */

/**
 * This class searches a profile HMM against a list of transcript sequences. This list
 * is composed by transcripts fused by DELTA extension.
 */
public class Find {

    private int firstSource;
    private int lastSink;

    public Find(int firstSource, int lastSink){
        this.firstSource = firstSource;
        this.lastSink = lastSink;

    }

    /**
     * This methods searches hits of a transcript.
     * @param allHits List of all hits to be updated recursively
     * @param phmm Profile HMM
     * @param t Transcript
     * @return List of updated hits
     */
    public static List<Hit> getHitsDP(List<Hit> allHits, ProfileHMM phmm, Transcript t){

        Translation translation = new Translation(t);
        String aminoSequence = translation.translate();

        if(aminoSequence == null)
            // Do not align if sequence is smaller than model length
            return allHits;

        // Creating a dynamic programming matrix object.
        DyPMatrix dyPMatrix = new DyPMatrix(aminoSequence, phmm);

        //Initializing the cells of the matrix
        dyPMatrix.initializeScoreMatrix();

        //align sequence
        if(!Astafunk.isLocal()){
            dyPMatrix.alignSequenceBnB();
        }else{
            dyPMatrix.alignSequenceLocal();
        }
        allHits = dyPMatrix.getHits(allHits, 0, t, t.getTranscriptID(), null, 0, 0);
        return allHits;

    }

    /**
     * Returns a list of hits founded in the amino acid sequences of transcripts, considering
     * DELTA factor boundaries and source and sink of a event.
     * This method iterates the transcripts to search a profile HMM.
     * @param hitList AS event hits. The first one is the hits
     *                   predicted in the entire sequence. AS event hits are hits predicted only
     *                   on the AS region and DELTA extension.
     * @param transcripts List of transcript (one or more transcripts)
     * @param hmm Profile HMM
     * @return Hashmap of hits/predictions
     */
    public List<Hit> getHitsInPosition(List<Hit> hitList, HashMap<String,Transcript> transcripts,
                                       ProfileHMM hmm, DirectedRegion [] eventRegions){

        //retrieve DELTA of the profile HMM
        int delta = hmm.getDeltaFactor(); //delta

        // Iterate the transcripts with AS events fused by DELTA extension.
        for(Map.Entry<String, Transcript> entry: transcripts.entrySet()){
            Transcript t = entry.getValue();
            String groupID = entry.getKey(); // string with IDs of genes/transcripts in the same variant.
            String aminoSequence = new Translation(t).translate(); // translate genomic sequence to amino acid sequence

            if(aminoSequence==null){ // if there is no valid CDS.
                if(Astafunk.isVerbose()) {
                    Log.message(Utils.printErrorOutput(t, groupID, 2));
                }
                continue;
            }

            int startPosition, endPosition;
            startPosition = t.getPrevCDSPosition(firstSource); // convert genomic to sequence coordinate.
            endPosition = t.getPrevCDSPosition(lastSink); // convert genomic to sequence coordinate.

            if(startPosition < 0)
                startPosition = 0;

            if(endPosition==-1) // before start 5'
                continue;

            // extend by DELTA*3 amino acids the first source and the last sink
            int downSubstringPosition = startPosition/3 - delta >= 0? startPosition/3 - delta : 0;
            int upSubstringPosition = endPosition/3 +delta <=
                    aminoSequence.length()? endPosition/3 +delta: aminoSequence.length();

            String sequenceAux = aminoSequence.substring(downSubstringPosition,upSubstringPosition);

            DyPMatrix dpm = new DyPMatrix(sequenceAux, hmm);
            if(!Astafunk.isLocal()) {
                dpm.alignSequenceBnB(); // glocal search (with branch-and-bound)
            }else{
                dpm.alignSequenceLocal(); // local search (without branch-ad-bound)
            }
            hitList = dpm.getHits(hitList, downSubstringPosition, t, groupID, eventRegions, firstSource, lastSink);

            if(hitList.size() == 0){
                if(Astafunk.isVerbose())
                    Log.message(Utils.printErrorOutput(t, groupID, 3));
            }
        }
        return hitList;
    }
}