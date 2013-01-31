package barna.astalavista.gfx;

import barna.commons.IntegerTupleComparator;
import barna.io.FileHelper;
import barna.io.gtf.GTFwrapper;
import barna.model.ASEvent;
import barna.model.DirectedRegion;
import barna.model.Species;
import barna.model.Transcript;
import barna.model.gff.GFFObject;

//import gphase.algo.FilterFactory;
//import gphase.gui.barna.astalavista.gfx.Circle;
//import gphase.gui.barna.astalavista.gfx.ColorMaster;
//
//import gphase.tools.barna.astalavista.gfx.Arrays;
//import gphase.tools.IntegerTupleComparator;
import org.freehep.util.export.ExportFileType;

import java.awt.*;
import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 10/2/12
 * Time: 10:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class LaVista {

    public static class FieldSizeComparator implements Comparator {

        public int compare(Object arg0, Object arg1) {

            int size1= -1, size2= -1;

            if (arg0 instanceof Vector)
                size1= ((Vector) arg0).size();
            else if (arg0 instanceof Object[])
                size1= ((Object[]) arg0).length;

            if (arg1 instanceof Vector)
                size2= ((Vector) arg1).size();
            else if (arg1 instanceof Object[])
                size2= ((Object[]) arg1).length;

            if (size1< size2)
                return (-1);
            if (size2< size1)
                return 1;
            return 0;
        }
    }

    public static class FieldSizeRevComparator extends FieldSizeComparator{
        public int compare(Object arg0, Object arg1) {
            return -super.compare(arg0, arg1);
        }
    }


    public static final int FILTER_NONE= 0;
    public static final int FILTER_HIERARCHICALLY= 1;
    public static final int FILTER_CODING_REDUNDANT= 2;
    public static final int FILTER_STRUCTURALLY= 3;
    public static final int FILTER_CONTAINED_IN_CDS= 4;
    public static final int FILTER_IDENTICALLY= 5;
    public static String[] FILTER_TO_STRING= {"none", "hierarchically", "coding_redundant", "structurally"};


    // /home/msammeth/annotations/hg17_RefSeqGenes_fromUCSC_070611_CDS_mRNAs_fromUCSC_070611_CDS_addDomains.gtf
    public static final String SFX_HTML= ".html", SFX_HEADER="header", SFX_TRAILER="trailer";
    public static final String SUBDIR_PICS= "pics";
    public static int limit= 10000;
    public static final HashMap<String, Color> EVENT_COLOR_MAP= new HashMap<String, Color>(5,1f);
    static {
        EVENT_COLOR_MAP.put("0,1-2^",new Color(58, 83, 164));
        EVENT_COLOR_MAP.put("1-,2-",new Color(237, 34, 36));
        EVENT_COLOR_MAP.put("1^,2^",new Color(16, 129, 64));
        EVENT_COLOR_MAP.put("1-2^,3-4^",new Color(246, 35, 122));
        EVENT_COLOR_MAP.put("0,1^2-",new Color(246, 235, 22));
        EVENT_COLOR_MAP.put("0,1-2^3-4^",new Color(100, 33, 101));
    };
    final static Color OTHERS_COL= new Color(192, 192, 192);
    public final static byte OVL_OPTION_VARIABLE= 0, OVL_OPTION_EXONIC= 1,OVL_OPTION_AREA= 2,OVL_OPTION_VAR_EXONIC= 3, OVL_MODUS_SINGLE= 0, OVL_MODUS_ALL= 1, OVL_MODUS_MANDATORY= 2, OVL_INTERSECT_OVERLAP= 0, OVL_INTERSECT_INCLUDE= 1;


    HashMap<String, Vector<DirectedRegion>>[] trptDomMap= null;

    public class HTMLWriterThread extends Thread {
        GFFObject[] objects= null;
        GFFObject[][] ovlObsEnd= null;
        ASEvent[] events= null;
        HashMap countMap= new HashMap();
        HashMap<String, ASEvent> eventMap= new HashMap<String, ASEvent>();
        HashMap<String, Vector<String>>[] summaryHashes= null;
        HashMap<String, Vector<String>> geneSummaryHash= null;

        public HTMLWriterThread() {
            setName("Writer Thread");
        }
        DirectedRegion[][] regsStart= null, regsEnd= null;

        private HashMap<DirectedRegion, DirectedRegion> getOvlRegionsAround(int pos, DirectedRegion[]  regs, DirectedRegion reg, HashMap<DirectedRegion, DirectedRegion> map) {

            if (pos< 0)
                pos= -(pos+1);

            for (int i = pos; i < regs.length; i++) {
                if (regs[i].overlaps(reg))
                    map.put(regs[i], regs[i]);
            }
            for (int i = pos-1; i>= 0; --i) {
                if (regs[i].overlaps(reg))
                    map.put(regs[i], regs[i]);
            }

            return map;
        }


        /**
         * gets from the sorted genomic region arrays the regions that overlap the event area
         * and belong to one of the transcripts in the event (regions without transcript identifier
         * belong per default to all transcripts).
         * @@param event
         * @@return
         */
        private HashMap<String, Vector<DirectedRegion>>[] getEventRegionMap(ASEvent event) {

            DirectedRegion evRegion= new DirectedRegion(event.getSrc().getPos(), event.getSnk().getPos(), event.getGene().getStrand());
            evRegion.setChromosome(event.getGene().getChromosome());

            // make genomic overlap with event
            DirectedRegion reg= event.getRegionEvent();
//			DirectedRegion[] regs= new DirectedRegion[] {reg};
            HashMap<String, Vector<DirectedRegion>>[] trptRegMaps= new HashMap[regsStart.length];
            String dummyID= "dummyID69";
            for (int i = 0; i < regsStart.length; i++) {

                // get overlap regions
//				if (LaVista.this.getOvlTypes()[i]== LaVista.OVL_OPTION_EXONIC)
//					regs= event.getExonicRegions();
//				else if (LaVista.this.getOvlTypes()[i]== LaVista.OVL_OPTION_VARIABLE)
//					regs= event.getVariableRegions();

                // get overlapping regions
                HashMap<DirectedRegion, DirectedRegion> map= new HashMap<DirectedRegion, DirectedRegion>();
                int p= java.util.Arrays.binarySearch(regsStart[i], reg, DirectedRegion.getDefaultPositionComparator());
                map= getOvlRegionsAround(p, regsStart[i], reg, map);
                p= java.util.Arrays.binarySearch(regsEnd[i], reg, DirectedRegion.getDefaultEndComparator());
                map= getOvlRegionsAround(p, regsEnd[i], reg, map);
//				for (int j = 0; j < regs.length; j++) {
//					int p= java.util.barna.astalavista.gfx.Arrays.binarySearch(regsStart[i], regs[i], DirectedRegion.getDefaultPositionComparator());
//					map= getOvlRegionsAround(p, regsStart[i], regs[i], map);
//					p= java.util.barna.astalavista.gfx.Arrays.binarySearch(regsEnd[i], regs[i], DirectedRegion.getDefaultEndComparator());
//					map= getOvlRegionsAround(p, regsEnd[i], regs[i], map);
//				}

                // all TIDs of event
                HashMap<String, String> tidMap= new HashMap<String, String>();
                for (int j = 0; j < event.getTranscripts().length; j++)
                    for (int k = 0; k < event.getTranscripts()[j].length; k++)
                        tidMap.put(event.getTranscripts()[j][k].getTranscriptID(),
                                event.getTranscripts()[j][k].getTranscriptID());

                // build hash
                trptRegMaps[i]= new HashMap<String, Vector<DirectedRegion>>();
                Iterator<DirectedRegion> iter= map.keySet().iterator();
                while (iter.hasNext()) {
                    DirectedRegion tmpReg= iter.next();
                    String tid= null;
                    if (tmpReg.getAttribute(GFFObject.TRANSCRIPT_ID_TAG)== null)
                        tid= dummyID;
                    else {
                        tid= (String) tmpReg.getAttribute(GFFObject.TRANSCRIPT_ID_TAG);
                        if (tidMap.get(tid)== null)
                            continue;
                    }

                    Vector<DirectedRegion> v= trptRegMaps[i].get(tid);
                    if (v== null) {
                        v= new Vector<DirectedRegion>();
                        trptRegMaps[i].put(tid,v);
                    }
                    v.add(tmpReg);
                }

                // extrapolate genomic features
                if (trptRegMaps[i].get(dummyID)!= null) {
                    Vector<DirectedRegion> v= trptRegMaps[i].remove(dummyID);
                    Iterator<String> it= tidMap.keySet().iterator();
                    while(it.hasNext())
                        trptRegMaps[i].put(it.next(), v);
                }

            }

            return trptRegMaps;
        }

        @Override
        public void run() {

            // finalize
            if (objects== null) {
                long t0= System.currentTimeMillis();
                writeHTMLgroupFinalizeMains();
                writeHTMLgroupFinalizeTrailers(eventMap, countMap);
                LaVista.this.writeHTMLOverview(countMap, eventMap);

                // write gene summary
                LaVista.this.writeHTMLSummary(LaVista.this.inFile.getName(), geneSummaryHash,
                        getLinkLoadMainFrame(
                                FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_"+ SFX_HEADER+ SFX_HTML,
                                FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_1"+ SFX_HTML,
                                FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_1_"+ SFX_TRAILER+ SFX_HTML
                        ));

                // make the summaries
                for (int i = 0; summaryHashes!= null&& i < summaryHashes.length; i++) {
                    LaVista.this.writeHTMLSummary(LaVista.this.ovlFiles[i].getName(), summaryHashes[i],
                            getLinkLoadMainFrame(
                                    FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_"+ SFX_HEADER+ SFX_HTML,
                                    FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_1"+ SFX_HTML,
                                    FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_1_"+ SFX_TRAILER+ SFX_HTML
                            ));
                }

                System.out.println("I painted the overview in "+((System.currentTimeMillis()- t0)/1000)+" sec.");

            } else {
                long t0= System.currentTimeMillis();
                for (int i = 0; i < objects.length; i++) {

                    // kill
                    objects[i].removeAttribute(GFFObject.ATTRIBUTE_TAG_LOCALIZATION);
                    objects[i].removeAttribute(GFFObject.ATTRIBUTE_TAG_LOCAL_REF);

                    ASEvent event= new ASEvent(objects[i]);

                    // check for overlap regions
                    HashMap<String, Vector<DirectedRegion>>[] evRegMaps= null;
                    if (regsStart!= null)
                        evRegMaps= getEventRegionMap(event);
                    if (geneSummaryHash== null)
                        geneSummaryHash= new HashMap<String, Vector<String>>();
                    if (summaryHashes== null)
                        summaryHashes= new HashMap[LaVista.this.ovlFiles.length];

//					if (!event.toStringStructure().equals("0,1-2^3-4^5-6^,3-4^5-6^,5-6^"))
//						continue;
                    String id= event.getAttribute(LaVista.this.getGroupAttributeTag());

                    int cnt= 0;
                    if (countMap.get(id)!= null)
                        cnt= ((Integer) countMap.get(id)).intValue();
                    if (cnt== 0)
                        writeHTMLeventGroupInit(event, evRegMaps, 0,
                                getLinkLoadMainFrame(
                                        FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_"+ SFX_HEADER+ SFX_HTML,
                                        FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_1"+ SFX_HTML,
                                        FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_1_"+ SFX_TRAILER+ SFX_HTML
                                ));
                    ++cnt;
                    int evOnPage= (cnt% eventPerPageLimit);
                    int pgIdx= cnt/ eventPerPageLimit;
                    String res= null;
                    if (evOnPage== 0) {
                        res= writeHTMLeventGroup(event, cnt, pgIdx-1, evRegMaps, geneSummaryHash, summaryHashes);
                        writeHTMLeventGroupInit(event, evRegMaps, pgIdx,
                                getLinkLoadMainFrame(
                                        FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_"+ SFX_HEADER+ SFX_HTML,
                                        FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_1"+ SFX_HTML,
                                        FNAME_OVERVIEW+"_"+getGroupAttributeTag() + "_1_"+ SFX_TRAILER+ SFX_HTML
                                ));
                    } else
                        res= writeHTMLeventGroup(event, cnt, pgIdx, evRegMaps, geneSummaryHash, summaryHashes);

                    if (res!= null) {
                        countMap.put(id, new Integer(cnt));
                        eventMap.put(id, event);
                    }

                    if (i== objects.length-1)
                        System.out.println(event.getGene().getChromosome()+" "+event.getSrc()+"->"+event.getSnk());
                }
                System.out.println("I designed and wrote "+objects.length+" objects in "+((System.currentTimeMillis()- t0)/1000)+" sec.");
            }


            events= null;	// just in case
        }

        public GFFObject[] getObjects() {
            return objects;
        }

        public void setObjects(GFFObject[] objects) {
            this.objects = objects;
        }

        public HashMap getCountMap() {
            return countMap;
        }

        public void setCountMap(HashMap countMap) {
            this.countMap = countMap;
        }

        public HashMap<String, ASEvent> getEventMap() {
            return eventMap;
        }

        public void setEventMap(HashMap<String, ASEvent> eventMap) {
            this.eventMap = eventMap;
        }

        public ASEvent[] getEvents() {
            return events;
        }

        public void setEvents(ASEvent[] events) {
            this.events = events;
        }

        public DirectedRegion[][] getRegsEnd() {
            return regsEnd;
        }

        public void setRegsEnd(DirectedRegion[][] regsEnd) {
            this.regsEnd = regsEnd;
        }

        public DirectedRegion[][] getRegsStart() {
            return regsStart;
        }

        public void setRegsStart(DirectedRegion[][] regsStart) {
            this.regsStart = regsStart;
        }
    }

    public class GroupReaderThread extends Thread {

        DirectedRegion[][] regsStart, regsEnd;
        HashMap<String, Vector<String>>[] summaryHashes = null;
        HashMap<String, Vector<String>> geneSummaryHash= null;

        public GroupReaderThread() {
            setName("Reader Thread");
        }

        private void initRegions(GFFObject[][] ovlObs) {
            if (ovlObs== null)
                return;

            // build up regions
            regsStart= new DirectedRegion[ovlObs.length][];
            regsEnd=  new DirectedRegion[ovlObs.length][];
            for (int i = 0; i < ovlObs.length; i++) {
                if (ovlObs[i]== null) {
                    regsStart[i]= new DirectedRegion[0];
                    regsEnd[i]= new DirectedRegion[0];
                } else {
                    regsStart[i]= new DirectedRegion[ovlObs[i].length];
                    regsEnd[i]= new DirectedRegion[ovlObs[i].length];
                }
                for (int j = 0; ovlObs[i]!= null&& j < ovlObs[i].length; j++) {
                    regsStart[i][j]= new DirectedRegion(ovlObs[i][j]);
                    regsEnd[i][j]= regsStart[i][j];
                }

                // sort
                java.util.Arrays.sort(regsStart[i], DirectedRegion.getDefaultPositionComparator());
                java.util.Arrays.sort(regsEnd[i], DirectedRegion.getDefaultEndComparator());
            }

        }

        public void run() {
            GTFwrapper reader= new GTFwrapper(LaVista.this.getInFile().getAbsolutePath());
            GTFwrapper[] ovlReaders= null;
            GFFObject[][] ovlObs= null;
            if (LaVista.this.ovlFiles!= null) {
                ovlReaders= new GTFwrapper[LaVista.this.ovlFiles.length];
                ovlObs= new GFFObject[LaVista.this.ovlFiles.length][];
            }
            reader.setGeneWise(false);
            reader.setReadFeatures(new String[] {ASEvent.GTF_FEATURE_ASEVENT});
            reader.setReadGene(false);
            reader.setReadGTF(true);
            reader.setLimitGTFObs(LaVista.limit);
            reader.setChromosomeWise(true);
            reader.setStrandWise(true);
            //reader.sweepToChromosome("chr16");
            if (!reader.isApplicable())
                return;
            reader.setSilent(true);


            HTMLWriterThread downstreamThread= new HTMLWriterThread();
            HashMap<String, ASEvent> eventMap= null;
            HashMap countMap= null;
            DirectedRegion[] ovlRegs= null;
            while (true) {
                long t0= System.currentTimeMillis();

                // read event obs
                try {
                    reader.read();
                } catch (Exception e) {
                    e.printStackTrace();
                }
                GFFObject[] obs= reader.getGtfObj();

                // get overlap chr
                if (obs!= null&& LaVista.this.ovlFiles!= null&&
                        (ovlRegs== null|| !ovlRegs[0].getChromosome().equals(obs[0].getSeqname())
                                || ovlRegs[0].getStrand()!= obs[0].getStrand())) {

                    ovlObs= new GFFObject[LaVista.this.ovlFiles.length][];

                    for (int i = 0; i < ovlFiles.length; i++) {
                        if (ovlReaders[i]== null) {
                            ovlReaders[i]= new GTFwrapper(LaVista.this.ovlFiles[i].getAbsolutePath());
                            ovlReaders[i].setGeneWise(false);
                            ovlReaders[i].setChromosomeWise(true);
                            ovlReaders[i].setStrandWise(true);
                            ovlReaders[i].setReadFeatures(null);
                            ovlReaders[i].setReadGene(false);
                            ovlReaders[i].setReadGTF(true);
                            ovlReaders[i].setSilent(true);
                            if (!ovlReaders[i].isApplicable())
                                return;
                        }
                        ovlReaders[i].sweepToChromosomeStrand(obs[0].getSeqname(), (byte) obs[0].getStrand());

                        try {
                            ovlReaders[i].read();
                        } catch (Exception e) {
                            e.printStackTrace();
                        }

                        ovlObs[i]= ovlReaders[i].getGtfObj();
                        ovlReaders[i].setGtfObj(null);
                    }

                    initRegions(ovlObs);
                    ovlObs= null;
                    System.gc();
                }

                // convert
//				ASEvent[] events= null;	// she does not like that..
//				if (obs!= null) {
//					events= new ASEvent[obs.length];
//					for (int i = 0; i < obs.length; i++) {
//						events[i]= new ASEvent(obs[i]);
//						obs[i]= null;
//					}
//					obs= null;
//					reader.setGtfObj(null);
//					System.gc();
//				}

                // stats
                if (obs!= null) {
                    System.out.print("I read "+obs.length+" event objects");
                    if (regsStart!= null) {
                        for (int i = 0; i < regsStart.length; i++)
                            if (regsStart[i]== null)
                                System.out.print(", 0 regions from "+ FileHelper.getFileNameWithoutExtension(LaVista.this.ovlFiles[i]));
                            else
                                System.out.print(", "+regsStart[i].length+" regions from "+ FileHelper.getFileNameWithoutExtension(LaVista.this.ovlFiles[i]));
                    }
                    System.out.println(" in "+((System.currentTimeMillis()- t0)/1000)+" sec.");
                }
//				else
//					System.out.println("finalizing now..");

                if (downstreamThread.isAlive())
                    try {
                        downstreamThread.join();
                    } catch (InterruptedException e1) {
                        ; //:)
                    }

//				if (obs[0].getChromosome().equals("chr10"))
//					obs= null;

                summaryHashes= downstreamThread.summaryHashes;
                geneSummaryHash= downstreamThread.geneSummaryHash;
                countMap= downstreamThread.getCountMap();
                eventMap= downstreamThread.getEventMap();
                downstreamThread= new HTMLWriterThread();
                downstreamThread.setObjects(obs);
                if (regsStart!= null) {
                    downstreamThread.setRegsStart(regsStart);
                    downstreamThread.setRegsEnd(regsEnd);
                }
                //downstreamThread.setEvents(events);
                //downstreamThread.setOvlRegions(ovlRegs);
                downstreamThread.setCountMap(countMap);
                downstreamThread.setEventMap(eventMap);
                downstreamThread.summaryHashes= summaryHashes;
                downstreamThread.geneSummaryHash= geneSummaryHash;
                downstreamThread.start();

                if (obs== null) {
                    try {
                        downstreamThread.join();
                    } catch (InterruptedException e1) {
                        ; //:)
                    }
                    break;
                }
            }
        }
    }

    final static String FNAME_OVERVIEW= "overview", FNAME_GROUP= "group";

    final static boolean LINKOUT_DOMAINS= true;

    final static String HTML_NAME_LANDSCAPE= "landscape.html";

    final static String UCSC_RESET_CART= "http://genome.ucsc.edu/cgi-bin/cartReset?destination=/cgi-bin/hgTables";
    final static String UCSC_LOAD_DEFAULTS= "hgS_doLoadUrl=submit;hgS_loadUrlName=http://genome.imim.es/~msammeth/customTracks/ucsc_def.txt";
    final static String UCSC_LOAD_CT= "hgt.customText=";
    final static String UCSC_CT_DOMAIN_URL= "http://genome.imim.es/~msammeth/customTracks/hg17/domains_Pfam"; // + / + RefTrpt.bed

    final static String FNAME_PIE= "distribution.png";
    final static String FNAME_LANDSCAPE= "landscape.html";
    final static String LOC_FLORET= "../pics/floret_cyan.jpg";
    final static String FNAME_REGIONS= "region.html";

    final static String ANNOTATION_DIR= "annotation";
    final static String UCSC_GB_CGI= "http://genome.ucsc.edu/cgi-bin/hgTracks?";
    public static final String[] SP_UCSC_CGI_STRINGS= new String[] {
            "knownGene=dense;encodeRegions=dense;encodeGencodeGeneOct05=pack;"
                    + "gap=hide;stsMap=hide;mgcGenes=hide;exoniphy=hide;exonWalk=hide;multiz17way=hide;snp125=hide;",	// human
            "",	// chimp
            "knownGene=dense;stsMap=hide;mgcGenes=hide;snp126=hide;multiz17way=hide;uniGene_3=hide;",	// Mouse
            "stsMapRat=hide;mgcGenes=hide;multiz9way=hide;netXenTro2=hide;netHg18=hide;snp125=hide;", // Rat
            "multiz4way=hide;netHg18=hide;blastHg18KG=hide;",		// &db=canFam2
            "all_mrna=pack;",		//Cow
            "multiz7way=hide;",	// Opossum

            "",	//Chicken
            "",	//X.+tropicalis
            "",	//Zebrafish
            "blastHg17KG=hide;cpgIsland=hide;blatHg16=hide;",			//Fugu
            "gaze=pack;netSelf=hide;",	//Tetraodon

            "flyBaseGene=dense;flyBaseNoncoding=pack;"
                    + "multiz17way=hide;blastHg17KG=hide;",	// Drosophila
            "",	//A.+gambiae
            "modelRefGene=dense;brhInparalog=hide;blastDm2FB=hide;netDm2=hide;",  //A.+mellifera

            "sangerGene=pack;sangerGenefinder=hide;c_briggsae_pwMaf=hide;wabaCbr=hide;axtNetCb1=hide;",	//C.+elegans
            "sgdGene=pack;sgdOther=hide;transRegCode=hide;esRegGeneToMotif=hide;multizYeast=hide;"	//S.+cerevisiae
    };
    final static String UCSC_STANDARD_PAR=
            // activate
            "ruler=dense;refGene=pack;mrna=pack;"
                    // deactivate
//		+"gap=hide;nscanGene=hide;intronEst=hide;rmsk=hide;"
                    // all species projection on
                    +"knownGene=pack;flyBaseGene=pack;modelRefGene=pack;sangerGene=pack;sgdGene=pack;" +
                    "gaze=pack;flyBaseNoncoding=pack;all_mrna=pack;" +
                    "encodeRegions=dense;encodeGencodeGeneOct05=pack;"
                    // all species projection off
//		+"gap=hide;cpgIsland=hide;" +
//				"exoniphy=hide;exonWalk=hide;" +
//				"stsMap=hide;stsMapRat=hide;sgdOther=hide;transRegCode=hide;esRegGeneToMotif=hide;" +
//				"mgcGenes=hide;uniGene_3=hide;sangerGenefinder=hide;" +
//				"snp125=hide;snp126=hide;" +
//				"multiz17way=hide;multiz9way=hide;multiz4way=hide;multiz7way=hide;multizYeast=hide;" +
//				"netXenTro2=hide;netHg18=hide;netSelf=hide;brhInparalog=hide;blastDm2FB=hide;netDm2=hide;c_briggsae_pwMaf=hide;wabaCbr=hide;axtNetCb1=hide;" +
//				"blastHg18KG=hide;blastHg17KG=hide;blatHg16=hide;"
                    +"";

    final static String TABLE_EVEN_COLOR= "#FFFFFF";
    final static String TABLE_ODD_COLOR= "#BFDFDF";
    final static String TABLE_HEADER_COLOR= "#00A0A0"; //"#008080";
    final static String HEADER_FILE= "header.ins", HEADER_FRAME_FNAME= "header_frame.ins", HEADER_FRAMES_FNAME= "header_frames.ins", HEADER_FNAME_BODY= "header_body.ins";
    final static String FNAME_MFRAME= "mframe.html";
    final static String TRAILER_FILE= "trailer.ins";
    final static String STYLE_FILE= "style.ins";
    final static int UCSC_FLANK= 50;
    public final static String GTF_ATTRIBUTE_GROUP_ID= "group_id";
    String groupAttributeTag= ASEvent.GTF_ATTRIBUTE_TAG_STRUCTURE, ovlFeature= null;
    boolean writeHTML= false,
            ovlOnly= true, ovlExcluded= false;
    Species species= null;
    File inFile= null, outDir= null;
    File[] ovlFiles= null;
    byte[] overlapTypes= null, overlapModes= null, overlapInter= null;
    String[] refTrptIDs= null;
    HashMap filterMap= new HashMap();

    // rubbish?
    static int filterCode= FILTER_HIERARCHICALLY;
    static int codingCode= ASEvent.TYPE_ALL;
    static String evFilterStr= "";

    static boolean killError= false;
    HashMap mapRegions= null;
    HashMap<String, Color>colTable= null;
    String inputFileName= null, colTableFName= ColorMaster.DEFAULT_COLOR_TABLE_FNAME;
    int eventPerPageLimit= 100;
    String[] customTracks= null;
    String[][] ovlFeatures= null;


    static void include(PrintStream p, String fName) {

        try {
            BufferedReader buffy= new BufferedReader(new FileReader(fName));
            while (buffy.ready())
                p.println(buffy.readLine());
            buffy.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }



    static void include(Writer writer, String fName) {

        try {
            BufferedReader buffy= new BufferedReader(new FileReader(fName));
            while (buffy.ready()) {
                writer.write(buffy.readLine());
                writer.write('\n');
            }
            buffy.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }



    static ASEvent[][] filter(ASEvent[][] vars, int codingC) {
        ASEvent[][] filtClasses= null;
        try {
            String mName= "";
            if (codingC== ASEvent.TYPE_ALL)
                mName= "isTrue";
            else if (codingC== ASEvent.TYPE_CDS)
                mName= "isProteinCoding";
            else if (codingC== ASEvent.TYPE_UTR)
                mName= "isNotAtAllCoding";
            else if (codingC== ASEvent.TYPE_5UTR)
                mName= "isCompletelyIn5UTR";
            else if (codingC== ASEvent.TYPE_3UTR)
                mName= "isCompletelyIn3UTR";

            Method m = vars[0][0].getClass().getMethod(mName, null);
            filtClasses= (ASEvent[][]) barna.astalavista.gfx.Arrays.filter(vars, m);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return filtClasses;
    }

    static String usage= "Hey!\nI'm LaVista, the little sister of Asta. I'm not so good at graph theory like her, I like arts better. \n" +
            "You still do not seem to know me, so I give you a brief summary how to use me:\n\n"+

            "LaVista [options] <input file>\n\n"+

            "-a , --attribute <attribute tag>\nthe attribute tag according to which I group the AS events in the file.\n" +
            "If you say nothing, I will sort them in structural classes as speciefied by the tag \'" +ASEvent.GTF_ATTRIBUTE_TAG_STRUCTURE+"\'\n" +
            "You know that the input file has to be sorted according to the tag provided here!\n\n"+
            "-o , --output <output directory>\nthe place you want me to put my artistic masterpieces. If you say nothing here, " +
            "I will write to the folder of the input file.\n\n"+
            "-x[v,e,w], --overlap[v,e] <gtf file>\noverlap with regions annotated in the gtf file," +
            "\n\t[v,e,w] overlap with Variable, Exonic, or Variable-Exonic regions of AS event (default is variable)\n\n"+
            "-c, --color <file>\nfile with color table\n\n"+
            "-g, --genome <species_genome>\nthe species name (e.g., human) and the genome version (e.g., hg18) in the form \"human_hg18\".\n\n"+
            "-ct, --customtrack <URL(s)>\nURL(s) where custom tracks are available in the form <URL>/chr/geneID\n\n" +

            "You got an basic idea now? If you want to know me better, read the manual on me by Micha.\n" +
            "Have a continued nice day..\n\nand \'asta\'!";

    static protected String missingArgumentError(String par) {
        return ("Dude, if you say \'" +par+ "\', you have to tell me a value for that!"+
                "For now, I'm ignoring \'"+par+"\', ok?!");

    }

    static protected LaVista parseArguments(String[] args) {
        if (args== null|| args.length== 0) {
            System.err.println(usage);
            System.exit(-1);
        }

        File f= new File(args[args.length-1]);
        LaVista myVista= null;
        if (f.exists()&& (FileHelper.getExtension(f).equalsIgnoreCase("gtf")||
                FileHelper.getExtension(f).equalsIgnoreCase("asta"))) {
            myVista= new LaVista(f);
            myVista.setOutDir(f.getParentFile());
        } else {
            System.err.println("Oh noo!! You didn't give me a valid input file as last argument!\n" +
                    f.getAbsolutePath()+" does not exist or does appears not to be a .gtf or .asta file.\n"+
                    "I cannot work like that! Start me without arguments to see what I expect from you.");
            System.exit(-1);
        }

        for (int i = 0; i < args.length-1; i++) {
            if (args[i].equals("-a")|| args[i].equals("--attribute")) {
                if (i+1>= args.length) {
                    System.out.println(missingArgumentError(args[i]));
                    continue;
                }
                myVista.setGroupAttributeTag(args[++i]);
                continue;
            }

            if (args[i].equals("-o")|| args[i].equals("--output")) {
                if (i+1>= args.length) {
                    System.out.println(missingArgumentError(args[i]));
                    continue;
                }
                File fout= new File(args[++i]);
                if (f.exists())
                    myVista.setOutDir(fout);
                else {
                    System.err.println("Ey man, there is no "+f.getAbsolutePath()+" folder!\n" +
                            "What do me expect to do? Byebye..");
                    System.exit(-1);
                }
                continue;
            }

            if (args[i].startsWith("-x")|| args[i].startsWith("--overlap")) {
                f= null;
                if (i+1< args.length)
                    f= new File(args[i+1]);
                if (f== null|| !f.exists()) {
                    System.err.println("Hey man, if you say '"+args[i]+"', you also have to give me a valid file with overlap data!");
                    if (f!= null)
                        System.err.println("I cannot find any "+f.getAbsolutePath()+".");
                    System.exit(-1);
                }

                int parsePos= 2;	// args[i].startsWith("-x")
                if (args[i].startsWith("--overlap"))
                    parsePos= 9;

                byte overlaptype= OVL_OPTION_AREA;
                byte ovlNb= 1;
                if (args[i].length()> parsePos) {
                    if (Character.isLetter(args[i].charAt(parsePos))) {
                        if (args[i].charAt(parsePos)== 'e') {
                            overlaptype= OVL_OPTION_EXONIC;
                            ++parsePos;
                        } else if (args[i].charAt(parsePos)== 'v') {
                            overlaptype= OVL_OPTION_VARIABLE;
                            ++parsePos;
                        } else if (args[i].charAt(parsePos)== 'w') {
                            overlaptype= OVL_OPTION_VAR_EXONIC;
                            ++parsePos;
                        }
                    }
                }
                byte overlapinter= OVL_INTERSECT_OVERLAP;
                if (args[i].length()> parsePos) {
                    if (args[i].charAt(parsePos)== 'i') {
                        overlapinter= OVL_INTERSECT_INCLUDE;
                        ++parsePos;
                    }
                }
                byte overlapmode= OVL_MODUS_SINGLE;
                if (args[i].length()> parsePos) {
                    if (args[i].charAt(parsePos)== 'a') {
                        overlapmode= OVL_MODUS_ALL;
                    } else if (args[i].charAt(parsePos)== 'm') {
                        overlapmode= OVL_MODUS_MANDATORY;
                    } else
                        System.out.println("Eiei, I don't understand overlap option '"+args[i].substring(parsePos)+"'! " +
                                "Allowed are e- exonic, v- variable; a- all. I'm ignoring that.");
                }

                myVista.addOvlFile(f);
                myVista.addOvlType(overlaptype);
                myVista.addOvlMode(overlapmode);
                myVista.addOvlInter(overlapinter);
                ++i;
                continue;
            }

            if (args[i].equalsIgnoreCase("-c")|| args[i].equalsIgnoreCase("--color")) {
                if (i+1== args.length) {
                    System.err.println("Hey man, if you say '"+args[i]+"', you also have to give me an argument!");
                    System.exit(-1);
                }
                myVista.setColTableFName(args[i+1]);
                ++i;
                continue;
            }

            if (args[i].startsWith("-ct")|| args[i].startsWith("--customtrack")) {

                int j = i+1;
                for (; j < args.length; j++)
                    if (!args[j].contains("://"))
                        break;

                if (j== i+1) {
                    System.err.println("Hey, you forgot to specify the customtrack!");
                    System.exit(-1);
                }

                int parsePos= 3;
                if (args[i].equalsIgnoreCase("--customtrack"))
                    parsePos= 13;

                for (++i; i < j; i++)
                    myVista.addCustomTrack(args[i]);
                continue;
            }

            if (args[i].equalsIgnoreCase("-g")|| args[i].equalsIgnoreCase("--genome")) {
                if (i+1== args.length) {
                    System.err.println("Hey man, if you say '"+args[i]+"', you also have to give me an argument!");
                    System.exit(-1);
                }
                String[] tokens= args[i+1].split("_");
                Species spe= new Species(tokens[0]);
                spe.setGenomeVersion(tokens[1]);
                myVista.setSpecies(spe);
                ++i;
                continue;
            }

            System.out.println(":s I do not understand "+args[i]+", sorry!! Maybe you mix me up with someone else.");
        }

        return myVista;
    }


    public static ASEvent[][] sortEvents(ASEvent[] events, Method m) {
        HashMap<Object, Vector<ASEvent>> map= new HashMap<Object, Vector<ASEvent>>();
        for (int i = 0; i < events.length; i++) {
            Object o= null;
            try {
                o= m.invoke(events[i], null);
            } catch (IllegalArgumentException e) {
                e.printStackTrace();
            } catch (IllegalAccessException e) {
                e.printStackTrace();
            } catch (InvocationTargetException e) {
                e.printStackTrace();
            }
            Vector<ASEvent> v= map.get(o);
            if (v== null) {
                v= new Vector<ASEvent>();
                map.put(o, v);
            }
            v.add(events[i]);
        }

        ASEvent[][] ev= new ASEvent[map.size()][];
        Iterator<Vector<ASEvent>> iter= map.values().iterator();
        int ctr= 0;
        while (iter.hasNext()) {
            Vector<ASEvent> v= iter.next();
            ASEvent[] e= new ASEvent[v.size()];
            for (int i = 0; i < e.length; i++)
                e[i]= v.elementAt(i);
            ev[ctr++]= e;
        }
        System.gc();

        java.util.Arrays.sort(ev, new barna.astalavista.gfx.Arrays.FieldSizeRevComparator());
        return ev;
    }

    ASEvent[] events;

    public LaVista(File inFile) {
        this.inFile= inFile;
    }


    public ASEvent[] getEvents() {

        if (events == null) {
            if (inFile== null)
                events= new ASEvent[0];
            else {
                GTFwrapper wrapper = new GTFwrapper(inFile.getAbsolutePath());
                try {
                    wrapper.read();
                } catch (Exception e) {
                    e.printStackTrace();
                }
                GFFObject[] obs= wrapper.getGtfObj();
                events= new ASEvent[obs.length];
                for (int i = 0; i < obs.length; i++) {
                    events[i]= new ASEvent(obs[i]);
                    obs[i]= null;
                }
                System.gc();
            }
        }

        return events;
    }



    public int run() {

        // check clean
        if (!outDir.exists())
            outDir.mkdir();
        String[] f= outDir.list();
        if (f.length> 0) {
            System.out.print("Ooops, there are files in the ouput folder "+outDir.getAbsolutePath()+".\n" +
                    "Do you want me to remove them before writing new files? If so, type \'yes\', if not I will write to that folder anyway, \n" +
                    "it's your fault if I overwrite something by doing that.\n");
            try {
                while (true) {
                    System.out.println("(yes/no/<CR>=stop)");
                    StringBuffer buf= new StringBuffer();
                    int c= System.in.read();
                    while (c!= '\n') {
                        buf.append(new Character((char) c).toString());
                        c= System.in.read();
                    }
                    String s= buf.toString().trim();
                    if (s.length()== 0|| s.equalsIgnoreCase("stop")) {
                        System.err.println("asta");
                        System.exit(-1);
                    }
                    if (s.equalsIgnoreCase("y")|| s.equalsIgnoreCase("yes")) {
                        FileHelper.rmDir(outDir);
                        outDir.mkdir();
                        break;
                    } else  if (s.equalsIgnoreCase("n")|| s.equalsIgnoreCase("no"))
                        break;
                }
            } catch (IOException e1) {
                e1.printStackTrace();
            }
        }

        // create pic subdir
        File d= new File(outDir+File.separator+SUBDIR_PICS);
        if (!d.exists())
            d.mkdir();


        // print stats

        // read
        long t0= System.currentTimeMillis();
        GroupReaderThread readerThread= new GroupReaderThread();
        readerThread.start();
        try {
            readerThread.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        if (colTable!= null)
            writeColTable();

        // write html
//		try {
//			Method m= null;
//			ASEvent.class.getMethod("toString", null);
//			writeHTML(events, m);
//		} catch (SecurityException e2) {
//			e2.printStackTrace();
//		} catch (NoSuchMethodException e2) {
//			e2.printStackTrace();
//		}

        return 0;
    }

    public static void main(String[] args) {

        //
        long t0= System.currentTimeMillis();
        LaVista myVista= parseArguments(args);
        myVista.run();
        long t1= System.currentTimeMillis();

        System.out.println("I finished in : "+ ((t1-t0)/ 1000)+" sec.\nasta!");
        System.exit(0);
    }

    public void writeASTA(ASEvent[][][] vars) {
        PrintStream p= null;
        try {
            p= new PrintStream(outDir+File.separator+ "landscape.asta");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        }
        for (int i = 0; vars!= null&& i < vars.length; i++)
            for (int j = 0; j < vars[i].length; j++)
                //for (int k = 0; k < vars[i][j].length; k++)
                p.println(vars[i][j][0].toStringASTA());
        p.flush(); p.close();
    }

    protected void writeHTMLStats(Writer writer, ASEvent[][] events) {
        // count all, pies
        return;
    }

    protected void writeHTMLOverviewHeader(HashMap countMap, HashMap<String , ASEvent> eventMap) {
        try {
            BufferedWriter writer= new BufferedWriter(new FileWriter(
                    outDir.getAbsolutePath()+ File.separator+ FNAME_OVERVIEW+"_"+getGroupAttributeTag() +"_"+SFX_HEADER+ SFX_HTML));
            include(writer, "header.ins");
            if (countMap.size()== 0)
                headline(writer, "NO AS EVENTS FOUND", null);
            else {
                headline(writer, "Landscape", null);
                writer.write("<div class=\"userspace\" align=\"center\"><img src=\""+SUBDIR_PICS+File.separator+FNAME_PIE+"\"></div>");
            }
            writer.write("</div></body>\n</html>");
            writer.flush();
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    protected void writeHTMLOverviewTrailer() {
        try {
            BufferedWriter writer= new BufferedWriter(new FileWriter(
                    outDir.getAbsolutePath()+ File.separator+ FNAME_OVERVIEW+"_"+getGroupAttributeTag() +"_"+SFX_TRAILER+ SFX_HTML));
            writer.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n" +
                    "<html>\n<head>\n</head>\n<body>\n");
            include(writer, "trailer_area.ins");
            writer.write("</body>\n</html>");
            writer.flush();
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    protected void writeHTMLSummary(String fName, HashMap<String, Vector<String>> hash, String linkback) {

        String relFName= "summary_"+fName;
        String baseFName= outDir+File.separator+relFName;
        String headerFName= baseFName+ "_"+SFX_HEADER;
        int pgCnt= 0;
        try {
            // write header
            BufferedWriter writer= new BufferedWriter(new FileWriter(headerFName+SFX_HTML));
            include(writer, "header.ins");
            headline(writer, "Summary", null);
            if (linkback!= null)
                writer.write("<DIV><a href=\""+linkback+"\"><IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\" alt=\"Landscape\">" +
                        "<FONT face=\"Arial,Lucida Sans\" size=\"1\">Landscape</FONT></a></DIV>\n");
            writer.write("</div><!-- closing main -->\n");
            writer.write("</body>");
            writer.write("</html>");
            writer.flush();
            writer.close();

            writer= null;

            // write mains and redirects
            Object[] keys= hash.keySet().toArray();
            // check for numeric sort
            if (keys.length> 0) {
                int x= 0;
                String s= (String) keys[0];
                for (; x < s.length(); x++)
                    if (Character.isLetter(s.charAt(x)))
                        break;
                if (x< s.length())
                    java.util.Arrays.sort(keys);
                else 	// numeric
                    java.util.Arrays.sort(keys, IntegerTupleComparator.getDefaultIntegerTupleComparator());


            }

            String currFName= null;
            for (int i = 0; i < keys.length; i++) {

                if (i% eventPerPageLimit== 0) {
                    if (writer!= null) {
                        writer.write("</TABLE>\n</div>\n</BODY></HTML>");
                        writer.flush(); writer.close();
                    }
                    ++pgCnt;
                    currFName= baseFName+"_"+pgCnt;
                    writer= new BufferedWriter(new FileWriter(currFName+SFX_HTML));
                    writer.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n"+
                            "<!-- micha: dont know why \"http://www.w3.org/TR/html4/loose.dtd\" makes context.js crash :s \n" +
                            " objects are accessible but no longer modifiable !?! .. java is an island, javascript not .. -->\n"+
                            "<html>\n<head>\n<SCRIPT type=\"text/javascript\" src=\"context.js\"></SCRIPT>\n</head>\n\n<body>\n<DIV align=\"center\">\n");

                    writer.write("<TABLE align=\"center\" cellspacing=\"0\" cellpadding=\"5\" border=\"0\" frame=\"void\">\n\n");
                    writer.write("<TR><TD colspan=\"3\" bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>"+fName+"</b></FONT></TD></TR>\n");
                }

                String colStr= null;
                if (i%2 == 0)
                    colStr= TABLE_EVEN_COLOR;
                else
                    colStr= TABLE_ODD_COLOR;

                // write summary of feature
                writer.write("<TR bgcolor=\""+colStr+"\">\n");
                writer.write("<TD align=\"left\"><a name=\"");
                writer.write(keys[i].toString());
                writer.write("\">");
                writer.write(keys[i].toString());
                writer.write("</a></TD><TD>");
                Vector<String> v= hash.get(keys[i]);
                writer.write(Integer.toString(v.size()));
                writer.write("</TD><TD align=\"left\"><TABLE><TR>");
                for (int j = 0; j < v.size(); j++) {
                    writer.write("<TD><TABLE frame=\"vsides\"><TR><TD>");	// second table for borders only at inner gutter
                    writer.write(v.elementAt(j));
                    writer.write("</TD></TR></TABLE></TD>");
                }
                writer.write("</TABLE></TD></TR>\n");

                // deviate linkback
                String currFooterFName= currFName+ "_"+ SFX_TRAILER;
                BufferedWriter writer2= new BufferedWriter(new FileWriter(outDir+File.separator+keys[i]+SFX_HTML));
                writer2.write("<HTML><HEAD>");
                writer2.write("<meta http-equiv=\"refresh\" content=\"0; URL=");
                writer2.write(getLinkLoadMainFrame(FileHelper.getFileNameOnly(headerFName + SFX_HTML),
                        FileHelper.getFileNameOnly(currFName + SFX_HTML) + "%23" + keys[i].toString(),
                        FileHelper.getFileNameOnly(currFooterFName + SFX_HTML)));
                writer2.write("\">");
                writer2.write("</HEAD></HTML>");
                writer2.flush();
                writer2.close();
            }
            writer.write("</TABLE>\n</div>\n</BODY></HTML>");
            writer.flush(); writer.close();

        } catch (Exception e) {
            e.printStackTrace();
        }


        // write footers
        for (int j = 0; j <= pgCnt; j++) {
            writeHTMLwriteTrailers(relFName, j);
        }

    }

    protected void writeHTMLOverview(HashMap countMap, HashMap<String , ASEvent> eventMap) {

        writeHTMLOverviewHeader(countMap, eventMap);

        if (countMap.size()!= 0)  try {
            //writeHTMLOverviewStats(writer, events);

            Pie pie= new Pie();
            //pie.setSize(new Dimension(300,300));
            //pie.setSize(pie.getPreferredSize());
            pie.setSize(new Dimension(690,300));
            pie.setBackground(OTHERS_COL);
            pie.init();
            int cntRank= 0, rankLimit= 10, sumRanks= 0, sum= 0, sumClasses= 0;

            Object[] oKey= countMap.keySet().toArray();
            for (int i = 0; i < oKey.length; i++) {
                String key= (String) oKey[i];
                Integer val= (Integer) countMap.remove(key);
                sum+= val.intValue();
                Vector v= (Vector) countMap.get(val);
                if (v== null)
                    v= new Vector();
                v.add(eventMap.get(key));
                countMap.put(val, v);	// invert
            }

            // overview
            ArrayList<Integer> freqs= new ArrayList<Integer>(countMap.keySet());
            Collections.sort(freqs);
            int outCnt= 0, innCnt= 0, totCnt= 0, pgCnt= 0;
            int lastEvCnt= Integer.MAX_VALUE;

            String relName= FNAME_OVERVIEW+"_"+ getGroupAttributeTag();
            BufferedWriter writer= null;
            for (int i = freqs.size()- 1; i>= 0;++totCnt) {

                if (totCnt% eventPerPageLimit== 0) {
                    if (writer!= null) {
                        writer.write("</TABLE></DIV>\n");
                        //writer.write("<br><br>");
                        writer.write("</body>\n</html>");
                        writer.flush();
                        writer.close();
                    }

                    ++pgCnt;
                    writer= new BufferedWriter(new FileWriter(outDir.getAbsolutePath()+ File.separator+ relName+ "_"+pgCnt + SFX_HTML));
                    writer.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n"+
                            "<html>\n<head>\n<SCRIPT type=\"text/javascript\" src=\"context.js\"></SCRIPT>\n</head>\n\n<body>\n<DIV align=\"center\">\n");

                    // table header
                    writer.write("<DIV class=\"userspace\" align=\"center\">\n");
                    writer.write("<TABLE bgcolor=\"#FFFFFF\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\" width=\"100%\">\n");
                    writer.write("<TR>\n");
                    writer.write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\" valign=\"middle\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Rank</b></FONT></TH>\n");
                    writer.write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>\n");
                    writer.write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\" valign=\"middle\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Proportion</b></FONT></TH>\n");
                    writer.write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>\n");
                    writer.write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Event<br>Count</b></FONT></TH>\n");
                    writer.write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>\n");
                    writer.write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Event<br>Details</b></FONT></TH>\n");
                    writer.write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>\n");
                    writer.write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"left\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Intron-Exon<br>Structure</b></FONT></TH>\n");
                    writer.write("</TR>\n");
                }

                // line col
                String colStr= null;
                if (totCnt%2 == 0)
                    colStr= TABLE_EVEN_COLOR;
                else
                    colStr= TABLE_ODD_COLOR;
                writer.write("<TR bgcolor=\""+colStr+"\">\n");

                // ranking
                Vector v= (Vector) countMap.get(freqs.get(i));
                if (freqs.get(i).intValue()< lastEvCnt) {
                    lastEvCnt= freqs.get(i).intValue();
                    ++outCnt; innCnt= 1;
                } else
                    ++innCnt;
                writer.write("\t<TD align=\"center\" valign=\"middle\" width=\"60\">"
                        + "<FONT size=\"4\">"+outCnt+"</FONT>");
                if (innCnt!= 1|| v.size()> 1)
                    writer.write("<FONT size=\"2\">."+innCnt+ "</FONT></TD>\n");
                writer.write("\t<TD></TD>\n");

                // init pie
                ASEvent event= (ASEvent) v.remove(0);
                ++sumClasses;
                if (cntRank< rankLimit) {
                    String struct= event.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_STRUCTURE);
                    String desc= ASEvent.EVENT_DESCRIPTION_MAP.get(struct);
                    if (desc== null)
                        desc= struct;
                    Color c= EVENT_COLOR_MAP.get(struct);
                    if (c== null)
                        pie.AddPieSlice(freqs.get(i).intValue(), desc);
                    else
                        pie.AddPieSlice(freqs.get(i).intValue(), desc, c);
                    ++cntRank;
                    sumRanks+= freqs.get(i).intValue();
                }

                // fraction
                float perc= 100f* freqs.get(i).intValue()/ sum;
                String percStr= Float.toString(perc);
                int cutoff= Math.min(percStr.indexOf('.')+ 3, percStr.length());
                percStr= percStr.substring(0, cutoff)+ "%";
                writer.write("\t<TD align=\"right\" valign=\"center\" width=\"60\">"+ percStr+ "</TD>\n");
                writer.write("\t<TD></TD>\n");

                // number
                writer.write("\t<TD align=\"right\" valign=\"middle\" width=\"60\">"+ freqs.get(i)+ "</TD>\n");
                writer.write("\t<TD></TD>\n");

                // structure
                String varFStr= getGroupAbsoluteFileName(event);	//varStr.replace(' ', '_');
                if (varFStr.indexOf(File.separator)>= 0)
                    varFStr= varFStr.substring(varFStr.lastIndexOf(File.separator)+1, varFStr.length());
                writer.write("\t<TD align=\"center\" valign=\"center\">"
                        + "<FONT face=\"Arial,Lucida Sans\" size=\"2\"><a href=\""
                        + getLinkLoadMainFrame(varFStr+"_"+SFX_HEADER+SFX_HTML, varFStr+"_1"+ SFX_HTML, varFStr+"_1_"+SFX_TRAILER+SFX_HTML)
                        +"\"></FONT>Show>></a>\n");
                writer.write("\t<TD></TD>\n");
                writer.write("\t<TD align=\"left\" valign=\"center\">"
                        //+ "<a href=\""+ varFStr+".html\">"
                        + "<img src=\""+ SUBDIR_PICS+File.separator+getGroupFileNameNormal(event)+".gif\"><br>"
                        // </a>
                        + "<FONT size=\"2\">code: </FONT><FONT size=\"2\" face=\"Courier\"><b>"+ event.getAttribute(getGroupAttributeTag())+ "</b></FONT>"
                        +"</TD>\n");

                writer.write("</TR>\n");

                if (v.size()== 0) {
                    countMap.remove(freqs);
                    --i;
                }
            }

            // paint pie
            pie.AddPieSlice((sum- sumRanks), "others ("+(sumClasses- rankLimit)+")", OTHERS_COL);
            pie.SetTitle("Overview");
            pie.setShow_values_on_slices(1);
//			JFrame frame= new JFrame();
//			frame.getContentPane().add(pie);
//			frame.setVisible(true);
            File f= new File(outDir.getAbsolutePath()+File.separator+SUBDIR_PICS+File.separator+FNAME_PIE);
            ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("png").get(0);
            pie.setPreferredSize(pie.getPreferredSize());
            t.exportToFile(f,pie,pie.getParent(),null,null);

            writeHTMLwriteTrailers(relName, pgCnt);

        } catch (IOException e) {
            e.printStackTrace();
        }

    }



    protected void writeHTMLovlLinkbacks(PrintStream p) {
        p.println("<a href=\""+HTML_NAME_LANDSCAPE+"\"><IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\" alt=\"Landscape\">" +
                "<FONT face=\"Arial,Lucida Sans\" size=\"1\">landscape</FONT></a>");
    }

    void writeHTMLOverviewPie(ArrayList<Integer> freqs, int totalCnt, HashMap countMap, int rankLimit) {

        // init pie
        Pie pie= new Pie();
        pie.setSize(new Dimension(300,300));
        pie.init();
        int cntRank= 0;
        for (int i = freqs.size()- 1; cntRank< rankLimit&& i>= 0; --i) {
            double frac= ((double) freqs.get(i).intValue())/ totalCnt;
            String struct= (String) countMap.get(freqs.get(i));
            String desc= ASEvent.EVENT_DESCRIPTION_MAP.get(struct);
            if (desc== null)
                desc= struct;
            pie.AddPieSlice(frac, desc, EVENT_COLOR_MAP.get(struct));
        }
        pie.SetTitle("Overview");
        pie.setShow_values_on_slices(1);
        pie.setSize(pie.getPreferredSize());
//		JFrame frame= new JFrame();
//		frame.getContentPane().add(pie);
//		frame.setVisible(true);
        try {
            File f= new File(outDir.getAbsolutePath()+File.separator+SUBDIR_PICS+File.separator+FNAME_PIE);
            ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("png").get(0);
            pie.setPreferredSize(pie.getPreferredSize());
            t.exportToFile(f,pie,pie.getParent(),null,null);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    static void headline(Writer writer, String headline, String linkBack) {
        //		<TD class="section" border=0 cellpadding=0 cellspacing=0 width=20% ALIGN=RIGHT>
        //		<a href="index.html#TOP"
        //		 onMouseover="window.status='TOP:INDEX';">
        //		<IMG class="pnt" SRC="http://genome.imim.es/g_icons/top.gif" HEIGHT=15 WIDTH=15 BORDER=0></a>
        try {
            writer.write("<TABLE border=0 cellpadding=0 cellspacing=0 width=\"100%\">\n");
            writer.write("<TR border=0 cellpadding=0 cellspacing=0>\n");
            if (linkBack!= null)
                writer.write("<TD class=\"section\" border=0 cellpadding=0 cellspacing=0 width=\"20%\"><a href=\""+linkBack+"\"><IMG class=\"pnt\" SRC=\"http://genome.imim.es/g_icons/top.gif\" HEIGHT=15 WIDTH=15 BORDER=0></a></TD>\n");
            writer.write("<TD class=\"section\" align=\"center\">\n");
            writer.write("<CENTER><FONT size=6 class=\"tgen\">"+headline+"</FONT></CENTER></TD>\n");
            writer.write("</TD></TR></TABLE>\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }



    static void headline(PrintStream p, String headline, String linkBack) {
//		<TD class="section" border=0 cellpadding=0 cellspacing=0 width=20% ALIGN=RIGHT>
//		<a href="index.html#TOP"
//		 onMouseover="window.status='TOP:INDEX';">
//		<IMG class="pnt" SRC="http://genome.imim.es/g_icons/top.gif" HEIGHT=15 WIDTH=15 BORDER=0></a>
        p.println("<TABLE border=0 cellpadding=0 cellspacing=0 width=\"100%\">");
        p.println("<TR border=0 cellpadding=0 cellspacing=0>");
        if (linkBack!= null)
            p.println("<TD class=\"section\" border=0 cellpadding=0 cellspacing=0 width=\"20%\"><a href=\""+linkBack+"\"><IMG class=\"pnt\" SRC=\"http://genome.imim.es/g_icons/top.gif\" HEIGHT=15 WIDTH=15 BORDER=0></a></TD>");
        p.println("<TD class=\"section\" align=\"center\">");
        p.println("<CENTER><FONT size=6 class=\"tgen\">"+headline+"</FONT></CENTER></TD>");
        p.println("</TD></TR></TABLE>");
    }

    void writeHTMLlandscapeLinkbacks(Writer w) {
        if (ovlFeature!= null) {
            try {
                w.write("<a href=\""+FNAME_REGIONS+"\"><IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\">" +
                        "<FONT face=\"Arial,Lucida Sans\" size=\"1\">"+ovlFeature+"</FONT></a>\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    static void _style_writeHTMLStats(ASEvent[][] vars, PrintStream p) {
        int spacer= 80;
        p.println("<HTML>");
        p.println("<HEAD>");
        p.println("<TITLE>AS Landscape</TITLE>");
        include(p, STYLE_FILE);
        p.println("</HEAD>");
        p.println("<BODY>");

        include(p, HEADER_FILE);
        if (vars!= null)
            p.println("<div class=\"title\"><h1>AS Landscape</h1></div><br />");
        else
            p.println("<div class=\"title\"><h1>NO AS EVENTS FOUND</h1></div><br />");

        p.println("<img src=\"distribution.png\">");
        int sum= 0;
        for (int i = 0; vars!= null&& i < vars.length; i++)
            sum+= vars[i].length;
        if (vars!= null) {
            p.println("<TABLE>");
            p.println("<TR>");
            p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"left\"><b>Rank</b></TH>");
            p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\"><b>Count</b></TH>");
            p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\"><b>Fraction</b></TH>");
            p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><b>Details</b></TH>");
            p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><b>Structure</b></TH>");
            p.println("</TR>");
            for (int i = 0; vars!= null&& i < vars.length; i++) {
                String colStr= "";
                if (i%2 == 0)
                    colStr= TABLE_EVEN_COLOR;
                else
                    colStr= TABLE_ODD_COLOR;
                p.println("<TR bgcolor=\""+colStr+"\">");
                p.println("\t<TD align=\"left\" valign=\"middle\" width=\"60\">%23"+ (i+1)+ "</TD>");
                p.println("\t<TD align=\"right\" valign=\"middle\" width=\"60\">"+ vars[i].length+ "</TD>");
                float perc= 100f* vars[i].length/ sum;
                String percStr= Float.toString(perc);
                int cutoff= Math.min(percStr.indexOf('.')+ 3, percStr.length());
                percStr= percStr.substring(0, cutoff)+ "%";
                p.println("\t<TD align=\"right\" valign=\"center\" width=\"60\">"+ percStr+ "</TD>");
                String varStr= vars[i][0].toString();
                String varFStr= convertFName(varStr);	//varStr.replace(' ', '_');
                p.println("\t<TD align=\"center\" valign=\"center\">"
                        + "<a href=\""+ varFStr+".html\">Show>></a>");
                p.println("\t<TD align=\"left\" valign=\"center\">"
                        + "<a href=\""+ varFStr+".html\"><img src=\""+ varFStr+".gif\"></a><br>"
                        + "code "+ varStr+ "</TD>");
                p.println("</TR>");
            }
            p.println("</TABLE>");
        }
        p.println("</div><!-- closing main -->");
        p.println("</BODY>");
        p.println("</HTML>");
    }

    protected static String getFileName(String in) {
        Integer name= (Integer) fileNameMap.get(in);
        if (name== null) {
            name= new Integer(++eventTypeNr);
            fileNameMap.put(in, name);
        }
        return name.toString();
    }



    protected static String convertFName(String in) {
        StringBuffer sb= new StringBuffer(in);
        // kill spaces
        for (int i = 0; i < sb.length(); i++)
            if (sb.charAt(i)== ' ')
                sb.deleteCharAt(i--);
        String out= sb.toString();
        out= out.replace('^', 'd');
        out= out.replace('-', 'a');
        out= out.replace(',', 'I');
        Integer name= (Integer) fileNameMap.get(out);
        if (name== null) {
            name= new Integer(++eventTypeNr);
            fileNameMap.put(out, name);
        }
        return name.toString();
        //		if (out.length()> 245)	{
        //			System.out.println("WARNING: shorted filename "+out);
        //			out= out.substring(0, 245);	// max fname len
        //		}
        //
        //		return out;
    }



    protected static String normalize(String in) {
        StringBuffer sb= new StringBuffer(in);
        // kill spaces
        for (int i = 0; i < sb.length(); i++)
            if (sb.charAt(i)== ' ')
                sb.deleteCharAt(i--);
        String out= sb.toString();
        out= out.replace('^', 'd');
        out= out.replace('-', 'a');
        out= out.replace('*', 's');
        out= out.replace(';', 'e');
        out= out.replace(',', '.');

        String name= null;
        if (in.length()> 100) {
            Integer intg= (Integer) fileNameMap.get(out);
            if (intg== null) {
                intg= new Integer(++eventTypeNr);
                fileNameMap.put(out, intg);
            }
            name= intg.toString();
        } else
            name= out.toString();

        return name;
    }



    protected static String getNormalizedFName(SpliceOSigner2 signer) {

        String name= new Integer(signer.hashCode()).toString();

        return name;
    }

    protected void writeHTMLeventColTPair(ASEvent var, PrintStream p) {
        p.print("\t<TD valign=\"top\"<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b>"+
                var.getTranscripts()[0][0]+" "+var.getTranscripts()[1][0]
                +"</b></TD>\n");
    }

    protected void writeHTMLeventGroupTableHeader(Writer writer, ASEvent event) {
        try{
            writer.write("<TR><TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Nr</b></FONT></TD>" +
                    //"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
                    "\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Position</b></FONT></TD>\n");
            if (event.getAttribute(GFFObject.GENE_ID_TAG)!= null)
                writer.write("\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Gene</b></FONT></TD>\n");
            if (event.getAttribute(GFFObject.ATTRIBUTE_TAG_LOCALIZATION)!= null)
                writer.write("\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Localization</b></FONT></TD>\n");
            writer.write("\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Attr</b></FONT></TD>\n");
            if (ovlFiles!= null) {
                for (int i = 0; i < ovlFiles.length; i++) {
                    writer.write("\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\"><b>");
                    String s= ovlFiles[i].getName();
                    int nbLines= (s.length()/ 15)+1;
                    if (s.length()> 15) {	// split in lines
                        for (int j = 0; j < nbLines; j++) {
                            for (int j2 = 0; j2 < 15; j2++) {
                                int pos= j*15+ j2;
                                if (pos>= s.length())
                                    break;
                                writer.write(s.charAt(pos));
                            }
                            if (j+1< nbLines)
                                writer.write("<br>");
                        }

                    } else
                        writer.write(s);
                    writer.write("</b></FONT></TD>\n");

                    writer.write("\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>IDs</b></FONT></TD>\n");
                }
            }
//			if (event.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_SPLICE_CHAIN)!= null)
//				writer.write("\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Coordinates</b></FONT></TD>\n");
            writer.write("\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Variants</b></FONT></TD>\n");
//			if (event.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_SOURCES)!= null)
//				writer.write("\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Source</b></FONT></TD>\n");
            if (event.getAttribute(Transcript.GTF_ATTRIBUTE_3COMPLETE_TAG)!= null)
                writer.write("\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>[3']</b></FONT></TD>\n");
//			if (event.getAttribute(Transcript.GTF_ATTRIBUTE_NMD_TAG)!= null)
//				writer.write("\t<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>NMD</b></FONT></TD>\n");

            //		if (ovlFeature!= null)
            //			writer.write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Overlap with</b></FONT>" +
            //						"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\">"+ovlFeature+"s</FONT></TD>"+
            //						"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>\n");

            writer.write("</FONT></TD></TR>\n");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    protected int getPos(String s, String q, int occNb) {
        int pos= 0;
        for (int j = 0; pos>= 0&& j < occNb; ++j)
            pos= s.indexOf(q, pos+1);
        if (pos< 0)
            pos= s.length();

        return pos;
    }

    protected String writeHTMLcontextMenu(String bgCol, Transcript[] list) {
        String[] ids= new String[list.length];
        for (int i = 0; i < ids.length; i++)
            ids[i]= list[i].getTranscriptID();
        return writeHTMLcontextMenu(bgCol, ids, true);
    }

    protected String writeHTMLcontextMenu(String bgCol, String[] list, boolean transcripts) {
        StringBuffer sb= new StringBuffer("\"context(event,'"+bgCol+"','");
        for (int j = 0; j < list.length; j++) {
            sb.append("<tr><td>");
            sb.append(list[j]);
            sb.append("</td></tr>");
        }
        sb.append("')\"");

        return sb.toString();
    }

    protected String writeHTMLeventGroup(ASEvent event, int rowNb, int pgIdx, HashMap<String, Vector<DirectedRegion>>[] evDomMap, HashMap<String, Vector<String>> geneSumHash, HashMap<String, Vector<String>>[] summaryMaps) {

//		The only fonts you can use reliably are
//
//		<font face="arial, helv, helvetica, sans-serif">
//		<font face="times, times roman, times new roman, tmsrm, serif">
//		<font face="courier, courier new, fixed">

        String colStr= null;
        if (rowNb%2 == 0)
            colStr= TABLE_EVEN_COLOR;
        else
            colStr= TABLE_ODD_COLOR;

        SpliceOSigner2[] signers= null;
        // evtl write domOverlayPicture
        if (evDomMap!= null) {
            signers= new SpliceOSigner2[evDomMap.length];
            for (int i = 0; i < signers.length; i++)
                signers[i]= writeHTMLEventPicture(event, evDomMap[i],getOvlTypes()[i], getOvlInter()[i], getOvlModes()[i] , true);
            if (signers[0]== null)	// first one is mandatory
                return null;
        }

        String fNameHTML= getGroupAbsoluteFileName(event)+"_"+(pgIdx+1)+SFX_HTML;
        try {
            BufferedWriter[] writers= new BufferedWriter[1];
            writers[0]= new BufferedWriter(new FileWriter(fNameHTML, true));
            //writers[1]= new BufferedWriter(new FileWriter(getGroupFileName(event)+"_fs_"+(pgIdx+1)+SFX_HTML, true));

            int start= (event.getGene().getStrand()>0?event.getSrc().getPos():Math.abs(event.getSnk().getPos()));
            int end= (event.getGene().getStrand()>0?event.getSnk().getPos():Math.abs(event.getSrc().getPos()));
            String posStr= event.getGene().getChromosome()+"&nbsp;"+GFFObject.getStrandSymbol(event.getGene().getStrand())+"&nbsp;"+
                    start+"-"+end;

            for (int x = 0; x < writers.length; x++) {
                for (int i = 0; i < event.getTranscripts().length; i++) {

                    writers[x].write("<TR bgcolor=\""+colStr+"\">");

                    // *** counter, source, pos ***
                    if (i== 0) {
                        String[] tIDs= new String[event.getTranscriptSupport()];
                        int c= 0;
                        for (int j = 0; j < event.getTranscripts().length; j++)
                            for (int k = 0; k < event.getTranscripts()[j].length; k++)
                                tIDs[c++]= event.getTranscripts()[j][k].getTranscriptID();

                        String ct= null;
                        if (customTracks!= null&& customTracks.length> 0)
                            ct= customTracks[0];
                        writers[x].write("<TD rowspan=\""+event.getTranscripts().length+"\" align=\"left\"><a name=\""+rowNb+"\"></a><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\">"+rowNb+"</FONT></TD>" +
                                //"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
                                "\n<TD rowspan=\""+event.getTranscripts().length+"\"><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\"><A href=\""+
                                getUCSCstring(species, ct, event.getGene().getChromosome(), event.getGene().getGeneID(),	// null -> customTrack
                                        start, end, tIDs)+"\" target=\"Genome Browser\">"+ //_blank
                                posStr+ "</A></FONT></TD>\n");

//						if (event.getAttribute(GFFObject.GENE_ID_TAG)!= null)
//							writers[x].write("\t<TD rowspan=\""+event.getTranscripts().length+"\"><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\">"+event.getGene().getGeneID()+"</FONT></TD>\n");

                        if (event.getAttribute(GFFObject.ATTRIBUTE_TAG_LOCALIZATION)!= null) {
                            writers[x].write("\t<TD rowspan=\""+event.getTranscripts().length+"\"><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\">\n\t\t");
                            writers[x].write(event.getAttribute(GFFObject.ATTRIBUTE_TAG_LOCALIZATION));
                            writers[x].write("</FONT></TD>\n");
                        }


                        // gene
                        String id= event.getAttribute(GFFObject.GENE_ID_TAG);
                        writers[x].write("<TD rowspan=\"");
                        writers[x].write(Integer.toString(event.getTranscripts().length));
                        writers[x].write("\" align=\"left\"><a href=\"");
                        writers[x].write(id+ SFX_HTML);
                        writers[x].write("\">");
                        writers[x].write(id);
                        writers[x].write("</a></TD>");
                        Vector<String> vsum= geneSumHash.get(id);
                        if (vsum== null) {
                            vsum= new Vector<String>(2,2);
                            geneSumHash.put(id,vsum);
                        }
                        vsum.add("<a href=\""+
                                getLinkLoadMainFrame(getGroupFileNameNormal(event)+ "_"+ SFX_HEADER+ SFX_HTML,
                                        FileHelper.getFileNameOnly(fNameHTML)+"%23"+rowNb,
                                        getGroupFileNameNormal(event)+"_"+(pgIdx+1)+ "_"+ SFX_TRAILER+ SFX_HTML)+"\">"
                                + "<img src=\""+signers[0].getFName()+"\"><br>"+posStr+"</a>");	// assume there is at least one overlap


                        // attributes
                        writers[x].write("\t<TD rowspan=\""+event.getTranscripts().length+"\"><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\">");
                        if (event.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_DIMENSION)!= null)
                            writers[x].write(event.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_DIMENSION));
                        else
                            writers[x].write(Integer.toString(event.getTranscripts().length)+".&#63;");	// = x.?
                        writers[x].write("&nbsp;");
                        writers[x].write(Integer.toString(event.getTranscriptSupport()));
                        writers[x].write("</FONT></TD>\n");

                    }



                    // *** pictures ***
                    if (i==0) {
                        for (int y = 0; signers!= null&& y < signers.length; y++) {
                            if (signers[y]== null) {
                                writers[x].write("<TD rowspan=\""+event.getTranscripts().length+"\" align=\"left\"></TD>");
                                writers[x].write("<TD rowspan=\""+event.getTranscripts().length+"\" align=\"left\"></TD>");
                            } else {
                                writers[x].write("<TD rowspan=\""+event.getTranscripts().length
                                        +"\" align=\"left\"><img src=\""+signers[y].getFName()+"\"></TD>\n");
                                // IDs
                                writers[x].write("<TD rowspan=\""+event.getTranscripts().length+"\" align=\"left\" nowrap>");
                                if (i== 0&& signers[y]!= null&& signers[y].getTrptDomMap()!= null) {
                                    Iterator<Vector<DirectedRegion>> iter= signers[y].getTrptDomMap().values().iterator();
                                    Vector<DirectedRegion> last= null;
                                    HashMap tmpMap= new HashMap();
                                    while (iter.hasNext()) {
                                        Vector<DirectedRegion> v= iter.next();
                                        for (int j = 0; j < v.size(); j++) {
                                            String id= v.elementAt(j).getID();
                                            if (v.elementAt(j).getAttribute(GTF_ATTRIBUTE_GROUP_ID)!= null)
                                                id= (String) v.elementAt(j).getAttribute(GTF_ATTRIBUTE_GROUP_ID);
                                            if (tmpMap.get(id)!= null)
                                                continue;	// (group) already collected

                                            tmpMap.put(id,id);
                                            Color c= getColTable().get(id);
                                            if (c== null) {
                                                c= ColorMaster.getRandomColor();
                                                getColTable().put(id, c);
                                            }
                                            String imgFName= writeDotPic(c);
                                            if (j!= 0|| last!= null)
                                                writers[x].write("<br>");
                                            writers[x].write("<img src=\""+imgFName+"\">&nbsp;");
                                            writers[x].write("<a href=\""+id+SFX_HTML+"\" target=\"_parent\">"+id+"</a>");
                                            // add to summary
                                            if (summaryMaps[y]== null)
                                                summaryMaps[y]= new HashMap<String, Vector<String>>();
                                            Vector<String> vsum= summaryMaps[y].get(id);
                                            if (vsum== null) {
                                                vsum= new Vector<String>(2,2);
                                                summaryMaps[y].put(id,vsum);
                                            }
                                            // backlink to feature
                                            vsum.add("<a href=\""+
                                                    getLinkLoadMainFrame(getGroupFileNameNormal(event)+ "_"+ SFX_HEADER+ SFX_HTML,
                                                            FileHelper.getFileNameOnly(fNameHTML)+"%23"+rowNb,
                                                            getGroupFileNameNormal(event)+"_"+(pgIdx+1)+ "_"+ SFX_TRAILER+ SFX_HTML)+"\">"
                                                    + "<img src=\""+signers[y].getFName()+"\"><br>"+posStr+"</a>");
                                        }
                                        last= v;
                                    }
                                    writers[x].write("</TD>\n");
                                }
                            }

                        }
                    }

                    // *** splice chains ***
//					if (event.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_SPLICE_CHAIN)!= null) {
//						writers[x].write("\t<TD><TABLE>");
//						writers[x].write("<TR><TD><FONT face=\"courier, courier new, fixed\" size=\"2\">");
//						for (int j = 0; j < event.getSpliceChains()[i].length; j++) {
//							writers[x].write(Integer.toString(Math.abs(event.getSpliceChains()[i][j].getPos())));
//							writers[x].write(event.getSpliceChains()[i][j].getSiteSymbol());
//							if (j< event.getSpliceChains()[i].length- 1)
//								writers[x].write(",");
//						}
//						writers[x].write("</FONT></TD></TR>");
//						writers[x].write("</TABLE></TD>\n");
//					}

                    // *** transcripts ***
                    writers[x].write("<TD><TABLE>");
                    writers[x].write("<TR><TD nowrap><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\">");
                    writers[x].write(Integer.toString(event.getTranscripts()[i].length));
                    writers[x].write(":&nbsp;");
                    StringBuffer buffy= new StringBuffer();
                    for (int j = 0; j < event.getTranscripts()[i].length; j++) {
                        buffy.append(event.getTranscripts()[i][j].getTranscriptID());
                        if (j< event.getTranscripts()[i].length- 1)
                            buffy.append(",");
                    }
                    //writers[x].write("<form><input name=\"\" type=\"text\" size=\"30\" maxlength=\""+buffy.length()+"\" value=\""+buffy.toString()+"\" readonly></form>");
                    if (buffy.length()< 33)
                        writers[x].write(buffy.toString());
                    else {
                        int pos= getPos(buffy.toString(), ",", 3);
                        String s= writeHTMLcontextMenu(colStr, event.getTranscripts()[i]);
                        //s= "\"context(event,'AA0000','<TR><TD>test</TD></TR>')\"";
                        writers[x].write(buffy.substring(0,pos)+"..<A style=\"cursor:default\" onclick="+s+
                                ">[<font color=\"#FF0000\"><b>+</b></FONT><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\">]</A>");
                    }
                    writers[x].write("</FONT></TD></TR>");
                    writers[x].write("</TABLE></TD>\n");

                    // *** sources ***
//					if (event.getAttribute(ASEvent.GTF_ATTRIBUTE_TAG_SOURCES)!= null) {
//						writers[x].write("\t<TD><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\">\n\t\t");
//						String[] sources= event.getGTFSources();
//						writers[x].write(sources[i]);
//						writers[x].write("</FONT></TD>\n");
//					}

                    // *** 3 complete ***
                    if (event.getAttribute(Transcript.GTF_ATTRIBUTE_3COMPLETE_TAG)!= null) {
                        writers[x].write("\t<TD><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\">\n\t\t");
                        String s= event.getAttribute(Transcript.GTF_ATTRIBUTE_3COMPLETE_TAG);
                        if (s.startsWith(","))	// split returning [0] for ","
                            s= " "+s;
                        if (s.endsWith(","))
                            s+= " ";
                        String[] tokens= s.split(",");
                        buffy= new StringBuffer();
                        String[] tokens2= null;
                        if (tokens.length> 0&& tokens[i].length()>0) {
                            tokens2= tokens[i].split("/");
                            for (int k = 0; k < tokens2.length; k++) {
                                buffy.append(tokens2[k]);
                                if (k< tokens2.length- 1)
                                    buffy.append(",");
                            }
                        }
                        if (buffy.length()< 33)
                            writers[x].write(buffy.toString());
                        else {
                            int pos= getPos(buffy.toString(), ",", 3);
                            s= writeHTMLcontextMenu(colStr, tokens2, true);
                            writers[x].write(buffy.substring(0,pos)+"..<A style=\"cursor:default\" onclick="+s+
                                    ">[<font color=\"#FF0000\"><b>+</b></FONT><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\">]</A>");
                        }
                        writers[x].write("</FONT></TD>\n");
                    }

                    // *** nmd ***
//					if (event.getAttribute(Transcript.GTF_ATTRIBUTE_NMD_TAG)!= null) {
//						writers[x].write("\t<TD><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\">\n\t\t");
//						String s= event.getAttribute(Transcript.GTF_ATTRIBUTE_NMD_TAG);
//						if (s.startsWith(","))	// split returning [0] for ","
//							s= " "+s;
//						if (s.endsWith(","))
//							s+= " ";
//						String[] tokens= s.split(",");
//						if (tokens[i].length()>0) {
//							String[] tokens2= tokens[i].split("/");
//							for (int k = 0; k < tokens2.length; k++) {
//								writers[x].write(tokens2[k]);
//								if (k< tokens2.length- 1)
//									writers[x].write(",");
//							}
//						}
//						writers[x].write("</FONT></TD>\n");
//					}


                    //	if (ovlFeature!= null)
                    //		writers[x].write("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"4\" color=\"#FFFFFF\"><b>Overlap with</b></FONT>" +
                    //					"<br><FONT face=\"arial, helv, helvetica, sans-serif\" size=\"2\" color=\"#FFFFFF\">"+ovlFeature+"s</FONT></TD>"+
                    //					"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>\n");

                    writers[x].write("</TR>\n");

                }


                writers[x].flush();
                writers[x].close();
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

        return fNameHTML+"%23"+rowNb;
    }





    protected String getGroupAbsoluteFileName(ASEvent event) {
        String fName= outDir.getAbsolutePath()+ File.separator+
                getGroupFileNameNormal(event);
        return fName;
    }

    protected String getGroupFileNameNormal(String code) {
        String fName=
                FNAME_GROUP+ "_"+ normalize(code);
        return fName;
    }

    protected String getGroupFileNameNormal(ASEvent event) {
        String fName= getGroupFileNameNormal(event.getAttribute(getGroupAttributeTag()));
        return fName;
    }

    public final static String UCSC_URL= "http://genome.ucsc.edu/cgi-bin/hgTracks?", UCSC_PAR_SPECIES= "org=", UCSC_PAR_GENOME_VERSION= "db=",
            UCSC_PAR_CUSTOM_TRACK= "hgt.customText=", UCSC_PAR_MARK_TRANSCRIPTS= "hgFind.matches=", UCSC_PAR_POSITION= "position=";

    public static String getUCSCstring(Species species, String ctURL, String chr, String geneID, int start, int end, String[] markTranscriptIDs) {
        StringBuffer sb= new StringBuffer(UCSC_URL);

        sb.append(UCSC_PAR_POSITION);
        sb.append(chr);
        sb.append(":");
        sb.append(start);
        sb.append("-");
        sb.append(end);
        sb.append(";");

        if (species!= null) {
            sb.append(UCSC_PAR_SPECIES);
            sb.append(species.getCommonName());
            sb.append(";");
            sb.append(UCSC_PAR_GENOME_VERSION);
            sb.append(species.getGenomeVersion());
            sb.append(";");
        }

        if (ctURL!= null) {
            sb.append(UCSC_PAR_CUSTOM_TRACK);
            sb.append(ctURL);
            sb.append("/");
            sb.append(chr);
            sb.append("/");
            sb.append(geneID);
            sb.append(".bed");
            sb.append(";");
        }

        if (markTranscriptIDs!= null) {
            sb.append(UCSC_PAR_MARK_TRANSCRIPTS);
            for (int i = 0; i < markTranscriptIDs.length; i++) {
                sb.append(markTranscriptIDs[i]);
                sb.append(",");
            }
            sb.deleteCharAt(sb.length()-1);
            sb.append(";");
        }

        return sb.toString();
    }


    SpliceOSigner2 writeHTMLEventPicture(ASEvent event, HashMap<String, Vector<DirectedRegion>> ovlMap, byte ovlType, byte ovlInter, byte ovlMode, boolean paintAnyway) {
        try {
            String fName= getGroupFileNameNormal(event);
            SpliceOSigner2 signer= new SpliceOSigner2(event);
            signer.setOverlapType(ovlType);
            signer.setOverlapMode(ovlMode);
            signer.setOverlapInter(ovlInter);
            signer.setColorTable(getColTable());
            signer.setTrptDomMap(ovlMap);
            if (!paintAnyway) {
                if (ovlOnly&& !signer.hasOverlappingRegions())
                    return null;
            }
            if (signer.hasOverlappingRegions())
                fName= getNormalizedFName(signer);
            signer.setSize(signer.getPreferredSize());
            ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
            File f= new File(outDir.getAbsolutePath()+File.separator+SUBDIR_PICS+File.separator+fName+".gif");
//			if (event.getSrc().getPos()== 14922162&& event.getSnk().getPos()== 14925376) {
//				JFrame myFrame= new JFrame();
//				myFrame.getContentPane().add(signer);
//				myFrame.pack();
//				myFrame.setVisible(true);
//			}
            t.exportToFile(f,signer,signer.getParent(),null,null);
            signer.setFName(SUBDIR_PICS+File.separator+fName+".gif");
            return signer;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    protected void writeHTMLeventGroupInitWriteHeader(ASEvent event, HashMap<String, Vector<DirectedRegion>>[] eventRegMaps, String linkback) {
        try {
            for (int i = 0; i < eventRegMaps.length; i++) {
                writeHTMLEventPicture(event, eventRegMaps[i], getOvlTypes()[i] , getOvlInter()[i], getOvlModes()[i], true);
            }

            String fName= getGroupAbsoluteFileName(event)+"_"+SFX_HEADER+SFX_HTML;
            BufferedWriter writer= new BufferedWriter(new FileWriter(fName));
            include(writer, "header.ins");
            headline(writer, "Event Group", null);

            if (linkback!= null)
                writer.write("<DIV><a href=\""+linkback+"\"><IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\" alt=\"Landscape\">" +
                        "<FONT face=\"Arial,Lucida Sans\" size=\"1\">Landscape</FONT></a></DIV>\n");

            // header
            writer.write("<DIV align=\"center\" class=\"userspace\" >");	// class=\"userspace\"  makes it white
            String groupID= event.getAttribute(getGroupAttributeTag());
            if (getGroupAttributeTag().equals(ASEvent.GTF_ATTRIBUTE_TAG_STRUCTURE))
                writer.write("<CENTER><img src=\""+ SUBDIR_PICS+File.separator+getGroupFileNameNormal(event)+".gif\"><br>\n");
            writer.write("<FONT size=\"2\" face=\"Arial,Lucida Sans\">code: </FONT><FONT size=\"2\" face=\"courier, courier new, fixed\"><b>"+ groupID+ "</b></FONT></CENTER><br><br>");
            writer.write("</DIV>\n");

            writer.write("</div><!-- closing main -->\n");
            writer.write("</body>");
            writer.write("</html>");
            writer.flush();
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }


    }

    protected void writeHTMLeventGroupWriteFrameSet(ASEvent event, int pgIdx) {
        try {
            BufferedWriter writer= new BufferedWriter(new FileWriter(getGroupAbsoluteFileName(event)+"_fs_"+(pgIdx+1)+SFX_HTML));
            //include(writer, HEADER_FRAMES_FNAME);
            writer.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Frameset//EN\" \"http://www.w3.org/TR/html4/frameset.dtd\"><html><head></head><frameset>\n");
            writer.write("<frame src=\""+getGroupAbsoluteFileName(event)+"_"+SFX_HEADER+SFX_HTML+"\" MARGINWIDTH=\"2\" MARGINHEIGHT=\"0\" name=\"header\">\n"); //HEADER_FILE
            writer.write("<frame src=\""+getGroupAbsoluteFileName(event)+"_"+(pgIdx+1)+SFX_HTML+"\" MARGINWIDTH=\"2\" MARGINHEIGHT=\"0\" name=\"table\">\n");
            writer.write("<frame src=\""+getGroupAbsoluteFileName(event)+"_"+(pgIdx+1)+"_"+SFX_TRAILER+SFX_HTML+"\" MARGINWIDTH=\"2\" MARGINHEIGHT=\"0\" name=\"trailer\">\n");
            writer.write("<noframes>\n");
            writer.flush();
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
    protected void writeHTMLeventGroupInit(ASEvent event, HashMap<String, Vector<DirectedRegion>>[] evRegMaps,int pgIdx, String linkback) {
        try {
            // header
            writeHTMLeventGroupInitWriteHeader(event, evRegMaps, linkback);

            // start event list
            String fName= getGroupAbsoluteFileName(event);
            BufferedWriter writer= new BufferedWriter(new FileWriter(fName+"_"+(pgIdx+1)+SFX_HTML));
            writer.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n"+
                    "<!-- micha: dont know why \"http://www.w3.org/TR/html4/loose.dtd\" makes context.js crash :s \n" +
                    " objects are accessible but no longer modifiable !?! .. java is an island, javascript not .. -->\n"+
                    "<html>\n<head>\n<SCRIPT type=\"text/javascript\" src=\"context.js\"></SCRIPT>\n</head>\n\n<body>\n<DIV align=\"center\">\n" +
                    "<TABLE align=\"center\" cellspacing=\"0\" cellpadding=\"5\" border=\"0\" frame=\"void\">\n\n");	// border=\"0\"  , width=\"100%\" cellpadding=\"0\"
            writeHTMLeventGroupTableHeader(writer, event);
            writer.flush();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    protected String getLinkLoadMainFrame(String header, String main, String trailer) {
        StringBuffer sb= new StringBuffer(FNAME_MFRAME+"?");
        if (header!= null) {
            sb.append("header=");
            sb.append(header);
            sb.append("&");
        }
        if (main!= null) {
            sb.append("main=");
            sb.append(main);
            sb.append("&");
        }
        if (trailer!= null) {
            sb.append("trailer=");
            sb.append(trailer);
        }
        return sb.toString();
    }

    protected void writeHTMLgroupFinalizeMains() {
        String[] files= outDir.list();
        String mainFile= ".*\\d+[.]html$";
        Matcher matty= null;
        for (int i = 0; i < files.length; i++) {
            if (!Pattern.matches(mainFile, files[i]))
                continue;
            try {
                BufferedWriter writer= new BufferedWriter(new FileWriter(outDir+File.separator+files[i], true));
                writer.write("</TABLE>\n</body>\n</html>");
                writer.flush();
                writer.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    protected void writeHTMLgroupFinalizeTrailers(HashMap<String, ASEvent> eventMap, HashMap countMap)  {

        Iterator<ASEvent> iter= eventMap.values().iterator();
        while (iter.hasNext()) {
            ASEvent event= iter.next();
            int evCnt= ((Integer) countMap.get(event.getAttribute(getGroupAttributeTag()))).intValue();
            int pgCnt= (evCnt/ eventPerPageLimit)+ 1;
            String fName= getGroupAbsoluteFileName(event);
            String relName= fName;
            if (fName.indexOf(File.separator)>= 0)
                relName= fName.substring(fName.lastIndexOf(File.separator)+1, fName.length());

            writeHTMLwriteTrailers(relName, pgCnt);
        }

    }

    void writeHTMLwriteTrailers(String relName,  int pgCnt) {
        String headerFN= relName+"_"+SFX_HEADER+SFX_HTML;
        for (int i = 0; i < pgCnt; i++) {
            try {
                BufferedWriter writer= new BufferedWriter(new FileWriter(outDir+File.separator+relName+"_"+(i+1)+"_"+SFX_TRAILER+SFX_HTML));
                writer.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n" +
                        "<html>\n<head>\n</head>\n\n<body>\n");
                writer.write("<!-- ############### X-INDEX AREA ############### -->\n");
                writer.write("<DIV align=\"center\">");

                if (i> 0)  // +1-1
                    writer.write("<a href=\""+getLinkLoadMainFrame(headerFN, relName+"_"+i+SFX_HTML, relName+"_"+i+"_"+SFX_TRAILER+SFX_HTML)
                            + "\"><img src=\"arrow_left.gif\" alt=\"left\"></a>&nbsp;");

                int min=(i> 10)?i-10:0;
                int max= i+10;
                if (max> pgCnt)
                    max= pgCnt;
                for (int j = min; j < max; j++) {
                    if (j== i)
                        writer.write("<b>");
                    else
                        writer.write("<a href=\""+getLinkLoadMainFrame(headerFN, relName+"_"+(j+1)+SFX_HTML, relName+"_"+(j+1)+"_"+SFX_TRAILER+SFX_HTML)+ "\">");

                    writer.write(Integer.toString(j+1));

                    if (j== i)
                        writer.write("</b>&nbsp;");
                    else
                        writer.write("</a>&nbsp;");
                }
                if (i< (pgCnt-1))
                    writer.write("<a href=\""+getLinkLoadMainFrame(headerFN, relName+"_"+(i+2)+SFX_HTML, relName+"_"+(i+2)+"_"+SFX_TRAILER+SFX_HTML)
                            + "\"><img src=\"arrow_right.gif\" alt=\"right\"></a>&nbsp;");
                writer.write("</DIV>");

                include(writer, "trailer_area.ins");
                writer.write("</body>\n</html>");
                writer.flush();
                writer.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

    }


    static int eventTypeNr= 0;
    static HashMap fileNameMap= new HashMap();

    public File getInFile() {
        return inFile;
    }



    public void setInFile(File inFile) {
        this.inFile = inFile;
    }



    public File getOutDir() {
        return outDir;
    }



    public void setOutDir(File outDir) {
        if (!outDir.exists())
            outDir.mkdir();
        this.outDir = outDir;
    }



    public String getOvlFeature() {
        return ovlFeature;
    }



    public void setOvlFeature(String ovlFeature) {
        this.ovlFeature = ovlFeature;
    }



    public boolean isWriteHTML() {
        return writeHTML;
    }



    public void setWriteHTML(boolean writeHTML) {
        this.writeHTML = writeHTML;
    }



    public String[] getRefTrptIDs() {
        return refTrptIDs;
    }



    public void setRefTrptIDs(String[] refTrptIDs) {
        this.refTrptIDs = refTrptIDs;
    }



    public boolean isOvlExcluded() {
        return ovlExcluded;
    }



    public void setOvlExcluded(boolean ovlExcluded) {
        this.ovlExcluded = ovlExcluded;
    }



    public boolean isOvlOnly() {
        return ovlOnly;
    }



    public void setOvlOnly(boolean ovlOnly) {
        this.ovlOnly = ovlOnly;
    }

    public void addFilter(String className, String filterString) {
        if (filterMap== null)
            filterMap= new HashMap();

        FilterFactory fac= (FilterFactory) filterMap.get(className);
        if (fac== null) {
            fac= new FilterFactory(className, filterString);
            filterMap.put(className, fac);
        } else
            fac.addMethodString(filterString);
    }



    public String getGroupAttributeTag() {
        return groupAttributeTag;
    }



    public void setGroupAttributeTag(String groupAttributeTag) {
        this.groupAttributeTag = groupAttributeTag;
    }



    public File[] getOvlFiles() {
        return ovlFiles;
    }



    public void addOvlFile(File ovlFile) {
        if (ovlFiles== null)
            ovlFiles= new File[] {ovlFile};
        else {
            File[] h= this.getOvlFiles();
            this.ovlFiles= new File[h.length+1];
            for (int i = 0; i < h.length; i++)
                this.ovlFiles[i]= h[i];
            this.ovlFiles[h.length]= ovlFile;
        }
    }



    public HashMap<String,Color> getColTable() {
        if (colTable == null) {
            if (colTableFName== null)
                colTable= ColorMaster.readInColorMap();
            else {
                colTable = new HashMap<String, Color>();
                File f= new File(getColTableFName());
                if (f.exists())
                    try {
                        BufferedReader reader= new BufferedReader(new FileReader(f));
                        while (reader.ready()) {
                            String line= reader.readLine();
                            String[] tokens= line.split("\t");
                            if (tokens.length!= 4)
                                continue;
                            Color c= new Color(Integer.parseInt(tokens[1]),
                                    Integer.parseInt(tokens[2]),
                                    Integer.parseInt(tokens[3]));
                            colTable.put(tokens[0], c);
                        }
                        reader.close();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
            }

        }

        return colTable;
    }

    public void writeColTable() {

        if (colTable== null|| colTableFName== null)
            return;

        try {
            BufferedWriter writer= new BufferedWriter(new FileWriter(getColTableFName()));
            Iterator<String> iter= colTable.keySet().iterator();
            StringBuffer sb;
            while (iter.hasNext()) {
                String key= iter.next();
                sb= new StringBuffer(key);
                sb.append("\t");
                Color c= colTable.get(key);
                sb.append(c.getRed());
                sb.append("\t");
                sb.append(c.getGreen());
                sb.append("\t");
                sb.append(c.getBlue());
                sb.append("\n");
                writer.write(sb.toString());
            }
            writer.flush();
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void setColTable(HashMap colTable) {
        this.colTable = colTable;
    }



    public String getColTableFName() {
        return colTableFName;
    }



    public void setColTableFName(String colTableFName) {
        this.colTableFName = colTableFName;
    }



    public String writeDotPic(Color c) {
        try {
            String imgID= Integer.toHexString(c.getRGB())+".gif";
            String relPath= SUBDIR_PICS+File.separator+imgID;
            File f= new File(outDir.getAbsolutePath()+File.separator+relPath);
            if (f.exists())
                return relPath;
            ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
            Component component= new Circle(c,8);
            component.setSize(component.getPreferredSize());
            t.exportToFile(f,component,component.getParent(),null,null);
            return relPath;
        } catch (Exception e) {
            // TODO: handle exception
        }
        return null;
    }



    public Species getSpecies() {
        return species;
    }



    public void setSpecies(Species species) {
        this.species = species;
    }



    public String[] getCustomTracks() {
        if (customTracks== null) {
            customTracks= new String[0];
        }
        return customTracks;
    }



    public void addCustomTrack(String customTrack) {
        if (customTracks== null)
            customTracks= new String[] {customTrack};
        else {
            String[] h= this.getCustomTracks();
            this.customTracks= new String[h.length+1];
            for (int i = 0; h!= null&& i < h.length; i++)
                this.customTracks[i]= h[i];
            this.customTracks[h.length] = customTrack;
        }
    }



    public byte[] getOvlTypes() {
        return overlapTypes;
    }

    public void addOvlMode(byte overlapMode) {
        if (overlapModes== null)
            overlapModes= new byte[] {overlapMode};
        else {
            byte[] h= this.getOvlModes();
            this.overlapModes= new byte[h.length+1];
            for (int i = 0; h!= null&&  i < h.length; i++)
                this.overlapModes[i]= h[i];
            this.overlapModes[h.length] = overlapMode;
        }
    }

    public void addOvlInter(byte val) {
        if (overlapInter== null)
            overlapInter= new byte[] {val};
        else {
            byte[] h= this.getOvlInter();
            this.overlapInter= new byte[h.length+1];
            for (int i = 0; h!= null&&  i < h.length; i++)
                this.overlapInter[i]= h[i];
            this.overlapInter[h.length] = val;
        }
    }


    public void addOvlType(byte overlapType) {
        if (overlapTypes== null)
            overlapTypes= new byte[] {overlapType};
        else {
            byte[] h= this.getOvlTypes();
            this.overlapTypes= new byte[h.length+1];
            for (int i = 0; h!= null&&  i < h.length; i++)
                this.overlapTypes[i]= h[i];
            this.overlapTypes[h.length] = overlapType;
        }
    }



    public byte[] getOvlModes() {
        return overlapModes;
    }



    public byte[] getOvlInter() {
        return overlapInter;
    }



    public void setOvlInter(byte[] overlapInter) {
        this.overlapInter = overlapInter;
    }
}