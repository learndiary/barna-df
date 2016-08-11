package barna.astafunk;

import barna.astalavista.AStaSettings;
import barna.astalavista.EventExtractor;
import barna.astafunk.Alignment.Find;
import barna.astafunk.DP.DyPMatrix;
import barna.astafunk.DP.Hit;
import barna.astafunk.HMM.ProfileHMM;
import barna.astafunk.utils.OrderingEventComparator;
import barna.astafunk.utils.Utils;
import barna.model.*;
import barna.model.splicegraph.SplicingGraph;

import java.util.*;

/**
 * @author vitorcoelho
 * @since 07/07/15
 * @version 2
 */
public class Tsearch implements Runnable{

    /**
     * Target gene of the search thread.
     */
    private Gene gene;

    /**
     * The splicing graph for this gene.
     */
    private SplicingGraph graph= null;

    /**
     * Maximum group number.
     */
    private int maxGroupNumber = 0;


    /**
     * Constructor of Search.
     * @param gene Gene object of the search thread.
     */
    public Tsearch(Gene gene) {
        this.gene = gene;
        graph = new SplicingGraph(gene);
    }


    /**
     * When an object implementing interface <code>Runnable</code> is used
     * to create a thread, starting the thread causes the object's
     * <code>run</code> method to be called in that separately executing
     * thread.
     * <p/>
     * The general contract of the method <code>run</code> is that it may
     * take any action whatsoever. Search method is selected by command line.
     *
     * @see Thread#run()
     */
    @Override
    public void run() {
        try {
            if (Astafunk.isTranscriptReferencePrint()) {
                printReferenceTranscript();
            } else if (Astafunk.isExhaustive()) {
                // Do not use heuristic table
                exhaustiveSearch();
            }else if(Astafunk.isConstitutive()) {
                searchConstitutiveDomains();
            }else if(Astafunk.isNaive()){
                naiveSearch();
            } else {
                defaultSearch();
            }

        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * Extract list of events from gene.
     * @param gene Gene
     * @return List of alternative splicing events of a gene
     */
    protected List<ASEvent> extractAllEvents(Gene gene) {

        AStaSettings astaSettings = new AStaSettings();

        /**
         * ASI: an internal AS event is flanked to both sides (i.e., at its 5'- and 3'-end) by common sites between the
         * compared transcripts, and does therefore not involve alterntative transcription start
         * ASE: an external AS event involves beside alternative splice site(s) also variability in the 5'- and/or
         * 3'-end of the transcript.
         */
        astaSettings.set(AStaSettings.EVENTS, EnumSet.of(AStaSettings.EventTypes.ASI,
                AStaSettings.EventTypes.ASE));

        astaSettings.set(AStaSettings.EVENTS_DIMENSION, -1);

        EventExtractor extractor = new EventExtractor(gene, astaSettings);

        //Constructing the graph
        extractor.constructGraph();

        //Contract the graph
        extractor.contractGraph(0);

        //Get events by partitions
        extractor.getEventsByPartitions(-1);
        List<ASEvent> l = extractor.getEventV();
        Collections.sort(l, new OrderingEventComparator());

        if(l.size()==0)
            return null;

        return l;
    }

    /**
     * Execute the heuristic search. Main mode. Full heuristic
     *
     * Features:
     * [X] Group+SplitVariant
     * [X] Domain reference file input
     * [ ] Search on entire sequence
     * [ ] Search reference transcript domains (no reference file input)
     *
     */
    private void defaultSearch() throws Exception {

        List<Hit> hitList = new ArrayList<Hit>();
        List<ASEvent> events;

        //Create a list of events of this gene
        events = extractAllEvents(gene);

        int num_events = events.size();

        if (num_events < 1) {
            return;
        }

        List<String> referencePfamList =
                retrieveReferenceList(gene.getTranscripts());
        /*
        tests if have a PFam reference list.
        */
        if (referencePfamList.size() == 0) {
            //Utils.finishedThread++;
            return;
        }

        /* Now, we have a list of Pfam domain hits for this gene.
        So, we must align only this models.
        */

        hitList = groupEvents(events, referencePfamList, hitList);

        Utils.printDefaultHit(hitList);
        //Utils.finishedThread++;
    }

    /**
     * Print reference transcripts from annotation.
     *
     * @throws Exception
     */
    private void printReferenceTranscript() throws Exception {

        List<ASEvent> l = extractAllEvents(gene);

        if (l.size() < 1) {
          //  Utils.finishedThread++; // only genes with AS
            return;
        }
        Transcript referenceTranscript = null;
        String referenceORF = "";

        for (Transcript t : gene.getTranscripts()) {
            String auxSequence = new Translation(t).translate();
            if (auxSequence != null) {
                if (auxSequence.length() > referenceORF.length()) {
                    referenceORF = auxSequence;
                    referenceTranscript = t;
                }
            }
        }

        if (referenceTranscript != null) {
            System.out.println(">" + referenceTranscript.getTranscriptID());
            System.out.println(referenceORF);
        }
      //  Utils.finishedThread++;
    }

    /**
     * Execute partial exhaustive search ( does not use heuristic table).
     * However, it creates a HMM reference list from the reference transcript of each gene.
     * So, the others transcript of this gene will be aligned with this reference list.
     * Features:
     * [X] Group+SplitVariant
     * [ ] Domain reference file input
     * [ ] Search on entire sequence
     * [X] Search reference transcript domains (no reference file input)
     *
     * @throws Exception
     */
    private void exhaustiveSearch() throws Exception {

        List<Hit> hitList = new ArrayList<Hit>();

        /**
         * Parsing profile HMM(s) file.
         */
        List<ASEvent> events;

        List<ProfileHMM> hmmAuxList = new ArrayList<ProfileHMM>();

        //Create a list of events of this gene
        events = extractAllEvents(gene);
        int num_events = events.size();

        if (num_events < 1) {
          //  Utils.finishedThread++;
            return;
        }

        List<String> referencePfamList = createReferenceList(gene);

         /*
        tests if have a PFam reference list.
         */
        if (referencePfamList == null) {
           // Utils.finishedThread++;
            return;
        }

        hitList = groupEvents(events, referencePfamList, hitList);
        Utils.printDefaultHit(hitList);
      //  Utils.finishedThread++;
    }

    /**
     * Retrieve the ACC of hits found by HMMER per transcript.
     * @param listT List of transcripts.
     * @return List of domains for the transcript
     */
    public List<String> retrieveReferenceList(Transcript[] listT) {
    /*
        A PFam reference list reduces our search space.
         */
        List<String> referencePfamList = new ArrayList<String>();

        for(Transcript t : listT){
            if (Astafunk.getHeuristicTable().containsKey(t.getTranscriptID())) {
                for (String acc : Astafunk.getHeuristicTable().get(t.getTranscriptID())) {
                    if (!referencePfamList.contains(acc)) {
                    /*
                    List with Acc's of models that we will search
                     */
                        referencePfamList.add(acc);
                    }
                }
                break;
            }
        }
        return referencePfamList;
    }

    /**
     * Create a domain reference list for a given gene (needs if there is no reference domain file input).
     * (1) Identify the reference transcript rtx of the gene.
     * (2) Search all input HMM against the sequence of this transcript rtx.
     * (3) Create a reference list of these HMMs.
     * @param gene Current gene to search functional impact of AS
     * @return List of profiles for the gene
     */
    public List<String> createReferenceList(Gene gene) {
    /*
        A PFam reference list reduces our search space.
         */
        List<String> referencePfamList = new ArrayList<String>();

    /* If there is not a domains reference list, we need a reference transcript
       A Reference Transcript is a transcript with the largest valid or annotated ORF
       in this gene.
    */
        Transcript referenceTranscript = null;
        String referenceORF = "";

        for (Transcript t : gene.getTranscripts()) {
            String auxSequence = new Translation(t).translate();

            if (auxSequence != null) {
                if (auxSequence.length() > referenceORF.length()) {
                    referenceORF = auxSequence;
                    referenceTranscript = t;
                }
            }
        }
    /*
     Gene does not have Pfam reference list or reference transcript.
    */
        if (referenceTranscript == null) {
            return null;
        }
    /*
    we must align all models from input to the largest ORF (reference ORF)
    The objective is to create a list with the possible domain models that
     can match with the transcripts of this gene and
    align only these models instead of align all models with all transcripts.
    */
        List<Hit> referenceHitList = new ArrayList<Hit>();
        for (Map.Entry<String, ProfileHMM> entry : Astafunk.getHmmHash().entrySet()) {
            referenceHitList = Find.getHitsDP(referenceHitList, entry.getValue(),
                    referenceTranscript);
        }

        // Now we have to create a list of Pfam ACC's to align the events
        for (Hit hit : referenceHitList) {
            String acc = hit.getAcc().split("\\.")[0];
            if (!referencePfamList.contains(acc))
                referencePfamList.add(acc);
        }

        if (referencePfamList.size() == 0) {
            if (Astafunk.isVerbose())
                Utils.printErrorOutput(gene.getTranscripts()[0], null, 4);
            return null;
        }
        return referencePfamList;
    }


    /**
     * Group events by delta overlapping and next calls the method to split the alternative
     * transcripts by variants
     * @param events List of alternative splicing events
     * @return Hashmap of hit lists
     */

    public List<Hit> groupEvents(List<ASEvent> events, List<String> referenceList, List<Hit> hitList){
        // check
        if (events.size()== 0)
            return null;    // nothing to do
        if (events.get(0).getGene()!= graph.getGene())
            throw new IllegalArgumentException("Gene mismatch "+ events.get(0).getGene()+ " vs "+ graph.getGene());

        // sort
        int numEvents = events.size();
        Collections.sort(events, new OrderingEventComparator());

        for(String acc: referenceList){
            ProfileHMM hmm = Astafunk.getHmmHash().get(acc);
            for(int i = 0; i < numEvents ; i++){

                // get transcript with 3'-most end, for heuristic of non-overlapping events
                Transcript[][] tt= events.get(i).getTranscripts();
                int pos3= Integer.MIN_VALUE;
                Transcript refTx= null;
                for(Transcript[] tarray: tt)
                    for(Transcript t : tarray)
                        if (t.get3PrimeEdge()> pos3) {
                            refTx= t;
                            pos3= t.get3PrimeEdge();
                        }
                int sinkE1genomic = events.get(i).getRegionEvent().getEnd();
                int sinkE1exonic= refTx.getExonicPosition(sinkE1genomic);

                // iterate following events, stop at the first that is not to be merged
                int j;
                int startEvent = i; //0
                for (j = i + 1; j < numEvents; j++){

                    Transcript tx= graph.getCommonTranscript(events.get(i), events.get(j));
                    if (tx== null) {
                        // cannot compute distance on transcript level
                        // => estimate how far away from the sink we are
                        // using the reference transcript
                        int sourceE2genomic= events.get(j).getRegionEvent().getStart();
                        int sourceE2exonic= refTx.getPrevExonicPosition(Math.abs(sourceE2genomic));
                        if ((sinkE1exonic+ (hmm.getDeltaFactor()* 3))< sourceE2exonic)
                            break;
                        // else
                        continue;
                    } else {
                        // update reference transcript and sink, for heuristics
                        refTx= tx;
                        sinkE1genomic= events.get(i).getRegionEvent().getEnd();
                        sinkE1exonic= refTx.getExonicPosition(sinkE1genomic);
                    }

                    // else
                    if(!hasDeltaOverlap(
                            events.get(i),             // event i
                            events.get(j),             // event j
                            tx,                        // a common transcript
                            hmm.getDeltaFactor()*3)) {  // amino acids -> nucleotide space
                        break;
                    }
                    i++; //go to next concatenated event
                }
                hitList = splitVariantsToAlign(events, startEvent, j- 1, hmm, hitList);
                i= (j-1);
            }
        }
        return hitList;
    }

    /**
     * Checks the overlap between AS event 1 and 2 with respect to delta.
     * @param asEvent1 an event
     * @param asEvent2 another event with overlapping transcript set
     * @param t a transcript common to <code>asEvent1</code> and <code>asEvent2</code>
     * @param delta the overlap region in <b>nucleotides</b> (!)
     * @return <code>true</code> if both AS events overlap by <code>delta</code>,
     *          <code>false</code> otherwise
     */
    private boolean hasDeltaOverlap(ASEvent asEvent1, ASEvent asEvent2, Transcript t, int delta) {

        if (t== null)
            return false;   // something else to say about that case?

        // coordinates can be negative for genes on reverse strand!
        int sinkE1genomic = asEvent1.getRegionEvent().getEnd();
        int sourceE2genomic = asEvent2.getRegionEvent().getStart();

        // get transcriptomic coordinates
        int sinkE1exonic= t.getExonicPosition(sinkE1genomic);
        int sourceE2exonic= t.getExonicPosition(sourceE2genomic);

        // if the following does not hold, something went bad in coordinate casting
        assert(sinkE1exonic<= sourceE2exonic);

        // decide
        return sinkE1exonic + delta >= sourceE2exonic;
    }


     /**
     * This method receives a list of events and align a transcript of these events.
     * @param events list with <b>all</b> AS events in <code>this</code> gene
     * @param first first index of events to be merged (included)
     * @param last last index of events to be merged (included)
     * @param hmm the HMM model to be aligned
     * @return Hashmap of hits. Key 0 - all hits; Key 1 - AS hits
     */
    public List<Hit> splitVariantsToAlign(List<ASEvent> events, int first, int last,
                                                           ProfileHMM hmm, List<Hit> hitList){

        DirectedRegion [] eventRegionList = getEventRegions(events, first, last);

        // get bounds of merged event
        Find f = createMergedFind(events, first, last);

        // split groups
        HashMap<String, Transcript> groupIDXtranscript; // group Id is transcript id's separeted by commas
        groupIDXtranscript = splitVariants(events, first, last); // return a hash with group id X transcript

        // align
        hitList = f.getHitsInPosition(hitList, groupIDXtranscript, hmm, eventRegionList);
        return hitList;
    }

    /**
     * Instantiates <code>Find</code> with the maximum boundaries of the events to be merged.
     * @param events list of <b>all</b> events in the gene
     * @param first first event to be merged (ordered by genomic positions)
     * @param last last event to be merged
     * @return a new instance of <code>Find</code> with the boundaries of the merged event set
     */
    private Find createMergedFind(List<ASEvent> events, int first, int last) {

        int firstStart = Math.abs(events.get(first).getRegionEvent().getStart());
        int lastEnd = Math.abs(events.get(last).getRegionEvent().getEnd());
        for (int i = (first+ 1); i < last; i++) {
            // this is for the general case of non-complete events
            // where events between first and last can extend beyond
            int start= Math.abs(events.get(first).getRegionEvent().getStart()),
                    end= Math.abs(events.get(first).getRegionEvent().getStart());
            if (start< firstStart)
                firstStart= start;
            if (end> lastEnd)
                lastEnd= end;
        }
        return new Find(firstStart, lastEnd);
    }

        /**
     * Method to return hashmap with [group ID | Transcript] where group ID is
     * a string with transcript id's separated by commas with the same splicing structure.
     * @param eventList List of <b>all</b> events in the gene
     * @param first index of the first event to be merged (included)
     * @param last index of the last event to be merged (included)
     * @return Hashmap with [Transcript ID's | Transcript]
     */
    public HashMap<String, Transcript> splitVariants(List<ASEvent> eventList,
                                                                 int first, int last){
        HashMap<String, Transcript> hashGroup;

        int num_events = (last- first)+ 1;
        if(num_events == 0)
            return null;

        // Hashmap with [Transcript  ID | Set]
        HashMap<String, Integer> transcriptHashMap;

        Transcript[][] transcriptMatrix = mergeTranscriptMatrices(eventList, first, last);

        // Takes all transcripts of gene presents in events
        transcriptHashMap = hashTranscripts(transcriptMatrix);

        // TranscriptHash = [Transcript_ID | Group_ID]
        for(Map.Entry<String, Integer>  element: transcriptHashMap.entrySet()){

            int group = element.getValue();

            for(Transcript[] variant: transcriptMatrix){
                List<String> listTranscript = new ArrayList<String>();
                boolean sameVariant = false; // test if they are in the same variant

                for(Transcript transcript: variant){ // transcripts of a variant
                    if(transcript.getTranscriptID().equals(element.getKey())){
                        sameVariant = true;
                        break;
                    }else{
                        listTranscript.add(transcript.getTranscriptID());
                    }
                }
                if(!sameVariant){
                    transcriptHashMap = updateGroups(transcriptHashMap, listTranscript, group);
                }
            }
        } // iteration of LinkedHashMap

        hashGroup = reduceGroups(transcriptHashMap);
        transcriptHashMap.clear();

        maxGroupNumber = 0;

        // So we take only one per group
        return hashGroup;
    }

    /**
     * Put all variants of concatenated events in a matrix.
     * @param eventList List of events concatenated by delta extension
     * @param first index of the first event sorted by transcriptomic positions
     * @param last index of the last event sorted by transcriptomic positions
     * @return Matrix of variants of all events
     */
    private Transcript[][] mergeTranscriptMatrices(List<ASEvent> eventList, int first, int last) {

        int size = 0;
        for(int i = first; i <= last; i++){
            size+= eventList.get(i).getTranscripts().length;
        }

        Transcript[][] transcriptMatrix = new Transcript[size][];

        int aux = 0;
        for(int i = first; i <= last; i++){
            System.arraycopy(eventList.get(i).getTranscripts(), 0, transcriptMatrix, aux,
                    eventList.get(i).getTranscripts().length);
            aux += eventList.get(i).getTranscripts().length;
        }
        return transcriptMatrix;
    }

    /**
     * Method to return a array with directed regions of a list of AS events.
     * @param events list of concatenated events
     * @param first index of the first event sorted by transcriptomic positions
     * @param last index of the last event sorted by transcriptomic positions
     * @return array with directed regions of events
     */
    private DirectedRegion[] getEventRegions(List<ASEvent> events, int first, int last) {

        DirectedRegion [] eventRegions = new DirectedRegion[(last- first)+1];
        int i =0;
        for(int j = first; j <=last; j++){
            eventRegions[i] = events.get(i).getRegionEvent();
            i++;
        }
        return eventRegions;
    }

    /**
     * Verify with transcripts with different exon/intron structure are in the same variant.
     * If yes, update the variant/group index and create other variant/group in the
     * hash table of transcripts and groups.
     * @param transcriptHashMap Hash Map with Transcript ID x Group ID
     * @param listTranscript List of transcript we must change the group id.
     * @param elementGroup Group id of reference
     * @return Updated hash map
     */
    private HashMap<String, Integer> updateGroups(HashMap<String, Integer> transcriptHashMap,
                                                  List<String> listTranscript, int elementGroup) {

        // Update the group number of a transcript in transcriptHashMap using newGroup if
        // the transcript is in listTranscript

        // maxGroupNumber should be global
        int newGroup = maxGroupNumber + 1;
        boolean hasNewGroup = false;

        for(String transcriptID : listTranscript){
           if(transcriptHashMap.get(transcriptID) == elementGroup) {
               transcriptHashMap.put(transcriptID, newGroup);
                hasNewGroup = true;
            }
        }
        if(hasNewGroup)
            maxGroupNumber++;
        return transcriptHashMap;
    }

    /**
     * Return a string ID of a group/variant. This string is composed by the transcript IDs
     * of the variant/group.
     * @param transcripts List of transcript IDs
     * @return A string of concatenated transcript IDs, representing a ID of a variant/group.
     */
    private static String getGroupID(List<String> transcripts){
        String groupID = transcripts.get(0);
        if(transcripts.size()==1)
            return groupID;
        for(int i = 1; i < transcripts.size(); i++){
            groupID += "," + transcripts.get(i);
        }
        return groupID;
    }

    /**
     * Check if a transcript belongs to a event.
     * @param transcriptID The transcript ID to be checked
     * @param event As event
     * @return true if the AS event contains the transcript
     */
    public static boolean containsTranscript(String transcriptID, ASEvent event){
        Transcript[][] matrixTranscript = event.getTranscripts();
        for(Transcript[] row: matrixTranscript){
            for(Transcript transcript: row){
                if(transcript.getTranscriptID().equals(transcriptID)){
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Method to initialize a hashmap with all transcripts ID of events as keys and
     * group/variant ID as value.
     * @return Hashmap with transcript ID as keys and group Id as value initialized with variant indexes equals
     * to 0.
     */
    private static HashMap<String, Integer> hashTranscripts(Transcript[][] tMatrix){
        HashMap<String, Integer> transcriptHash = new HashMap<String, Integer>();
        for(Transcript[] variant: tMatrix){
            for(Transcript transcript: variant) {
                    if (!transcriptHash.containsKey(transcript.getTranscriptID())){
                        transcriptHash.put(transcript.getTranscriptID(), 0);
                    }
            }
        }
        return transcriptHash;
    }


    /**
     * For each variant, create a groupID string (all transcripts IDs in this variant concatenated)
     * and retrieve a transcript with valid ORF to represent this variant and to be aligned.
     * @param variantHash
     * @return
     */
    //TODO Must be optimized (to remind)
    private HashMap<String, Transcript> reduceGroups(HashMap<String, Integer> variantHash){
        HashMap<Integer, List<String>> variants = new HashMap<Integer, List<String>>();
        // Initialize hash of variant IDs with empty transcript list
        for(Map.Entry<String, Integer> entry : variantHash.entrySet()){
            if(!variants.containsKey(entry.getValue())){
                variants.put(entry.getValue(), new ArrayList<String>());
            }
        }

        // Add a list of transcript IDs for each variant ID
        for(Map.Entry<Integer, List<String>> entry : variants.entrySet()){
            for(Map.Entry<String, Integer> entryAux : variantHash.entrySet()){
                if(entryAux.getValue().equals(entry.getKey())){ // if group ID of variantHash == variants
                    if(!entry.getValue().contains(entryAux.getKey())){ // if list of groups does not contain transcript
                        entry.getValue().add(entryAux.getKey());
                    }
                }
            }
        }

        //For each variant, create a groupID string (all transcripts IDs in this variant concatenated)
        //and retrieve a transcript with valid ORF to represent this variant and to be aligned.
        HashMap<String, Transcript> groupTranscript = new HashMap<String, Transcript>();

        for(Map.Entry<Integer, List<String>> entry : variants.entrySet()){
            String groupID = getGroupID(entry.getValue());
            Transcript t = getValidTranscript(entry.getValue());
            groupTranscript.put(groupID, t);
        }
        return groupTranscript;
    }

    private void showFinalVariants(HashMap<String, Transcript> groupTranscript) {
        System.out.println("#>Show final variants...");
        for(Map.Entry<String, Transcript> entry: groupTranscript.entrySet()){
            System.out.println("#>" + entry.getKey() + "\t" + entry.getValue());
        }
    }

    /**
     * Select a valid transcript, i.e., a transcript with amino acid sequence.
     * @param transcriptIds List of transcript IDs
     * @return Transcript with CDS regions and amino acid sequence.
     */
    private Transcript getValidTranscript(List<String> transcriptIds) {
        Transcript returnTranscript = null;
        boolean hasTranslatedSequence = false;
        for(String transcriptID : transcriptIds){
            for(Transcript t : graph.getGene().getTranscripts()){
                if(t.getTranscriptID().equals(transcriptID)){
                    returnTranscript = t;
                    if(new Translation(t).translate()!= null){
                        hasTranslatedSequence = true;
                        break;
                    }
                }
            }
            if(hasTranslatedSequence) break;
        }
        return returnTranscript;
    }

    /**
     * Method for debugging. Print transcript of an AS event.
     * @param event AS event.
     */
    private static void showTranscripts(ASEvent event){

        Transcript [][] transcMatrix = event.getTranscripts();

        int i = 0;
        for(Transcript [] row: transcMatrix){
            int j = 0;
            for(Transcript transc: row){

                System.out.println("[" + i +"] " +transc.getTranscriptID());

                j++;
            }
            i++;
        }
    }

    /**
     * Method for debugging. Print transcript of an list.
     * @param list List of transcripts.
     */
    private static void showTranscript2align(List<Transcript> list){

        System.out.println("\nBegin showing transcript to align...");

        for(Transcript element: list){

            System.out.println(element.getTranscriptID());

        }
        System.out.println("\nEnd showing transcript to align.\n");
    }

    /**
     * Method for debugging. Print transcript separated by the respective event.
     * @param eventList List of AS events
     */
    private static void showTrancriptsPerEvents(List<ASEvent> eventList){

        System.out.println("\nBegin showing transcripts per event...");
        for(ASEvent event: eventList){
            System.out.println(event.toString());
            showTranscripts(event);
            System.out.println("___________________________________________");

        }

        System.out.println("\nEnd showing transcripts per event.\n");

    }

    /**
     * Method for debugging. Print hashmap of transcript groups.
     * @param list Hashmap with transcript ID as key and group Id as value.
     */
    private static void showHashing(HashMap<String, Integer> list){

        System.out.println("#>Begin showing hash...");

        for(Map.Entry<String, Integer> element: list.entrySet()){
            System.out.println("#>"+element.getKey() + " | " + element.getValue());
        }

        System.out.println("#>End showing hash.\n");
    }

    /**
     * Search domains on constitutive regions of AS reference transcripts and non-AS genes.
     * [Search mode to obtain paper's results]
     */
    private void searchConstitutiveDomains() throws Exception {

        List<Hit> hitList = new ArrayList<Hit>();

        //Create a list of events of the AS genes
        List<ASEvent> events = extractAllEvents(gene);

        // Transcript with the longest ORF. If there is more than one transcript
        // with the longest ORF, select randomly a transcript.
        Transcript referenceTranscript = null;
        String referenceORF = "";

        for (Transcript t : gene.getTranscripts()) {
            String auxSequence = new Translation(t).translate();
            if (auxSequence != null) {
                if (auxSequence.length() > referenceORF.length()) {
                    referenceORF = auxSequence;
                    referenceTranscript = t;
                }
            }
        }

        if(referenceORF.equals("")) {
            // Utils.finishedThread++; //TODO debug
            return;
        }

        List<String> referencePfamList = retrieveReferenceList(gene.getTranscripts());

        /*
        tests if have a PFam reference list.
        */
        if (referencePfamList == null) {
            Utils.finishedThread++;
            return;
        }

        /* Now, we have a list of Pfam domain hits for this gene.
        So, we must align only this models.
        */

        for(String acc: referencePfamList){
            DyPMatrix dp = new DyPMatrix(referenceORF, Astafunk.getHmmHash().get(acc));
            dp.alignSequenceBnB();
            hitList = dp.getConstitutiveHits(hitList, referenceTranscript, events);
        }

        Utils.printDefaultHit(hitList);
        //Utils.finishedThread++; //TODO just for debug
    }

    /**
     * Search all domains in the reference transcript of an alternatively spliced gene.
     * The reference transcript is the one with the longest coding sequence of the gene.
     * Then, it iterates all the transcripts of the same gene searching the domains found
     * in the reference transcript.
     */
    private void naiveSearch(){

        List<ASEvent> events;

        //Create a list of events of this gene
        events = extractAllEvents(gene);
        int num_events = events.size();

        /*
        It does not considered events without AS events.
         */
        if (num_events < 1) {
            //Utils.finishedThread++; //TODO just for debug
            return;
        }

        /*
        Retrieve a domain list of the reference transcript (reference domain list)
        for the reference domain file (Heuristic table) given in the command line.
         */

        List<String> referencePfamList =
                retrieveReferenceList(gene.getTranscripts());
        /*
        Test if the reference transcript has domain.
         */
        if (referencePfamList.size() == 0) {
            //Utils.finishedThread++; //TODO debug
            return;
        }

        /**
         * Iterates transcripts of the AS gene.
         */
        for(Transcript t: gene.getTranscripts()){
            List<Hit> hitList = new ArrayList<Hit>();
            /**
             * Retrieve the HMMs based on reference domain list.
             */
            for(String acc: referencePfamList){
                /**
                 * Translate the coding sequence.
                 */
                String auxSequence = new Translation(t).translate();
                if (auxSequence != null) {
                    // search hmm
                    DyPMatrix dp = new DyPMatrix(auxSequence, Astafunk.getHmmHash().get(acc));
                    dp.alignSequenceBnB();
                    hitList = dp.getHits(hitList, 0, t, null, null, 0, 0);
                }
            }
            Utils.printDefaultHit(hitList);
        }
       // Utils.finishedThread++; //TODO just for debug
    }
}