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

package fbi.genome.io.gff;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.file.FileHelper;
import fbi.commons.io.IOHandler;
import fbi.commons.io.IOHandlerFactory;
import fbi.commons.tools.Sorter;
import fbi.genome.model.commons.MyArrayHashMap;
import fbi.genome.model.commons.MyFile;
import fbi.genome.model.gff.GFFObject;

import java.io.*;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.Future;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Sort a GTF file. Use the static {@link #sort(java.io.File)}  or {@link #sort(java.io.File, java.io.File)} methods
 * to use this class. The entries are sorted in the following order:
 *
 * <pre>
 *  1. Chromosome (Field 0)
 *  2. Strand (Field 6, + Strand first)
 *  3. Global Transcript Position (the minimal position of a set of transcript with equal chr and transcript ID)
 *  4. Transcript ID
 *  5. Position (Field 3)
 * </pre>
 *
 * @author Thasso Griebel (thasso.griebel@googlemail.com)
 * @author Micha Sammeth (gmicha@googlemail.com)
 */
public class GFFSorter {
    /**
     * Line splitter pattern
     */
    static final Pattern SPLITTER_PATTERN = Pattern.compile("\\s");
    /**
     * Pattern to findTranscriptID the transcript ID
     */
    static final Pattern TRANSCRIPTID_PATTERN = Pattern.compile(GFFObject.TRANSCRIPT_ID_TAG);

    /**
     * The line separator for the current file
     */
    private String eol;

    /**
     * INTERNAL : use the static methods like {@link #sort(java.io.File)} to use this class
     */
    private  GFFSorter() {
    }

    /**
     * Merge the two strings to one byte[]
     *
     * @param tid first string
     * @param chr second string
     * @return bytes merged string as byte[]
     */
	static byte[] encode(ByteArrayCharSequence tid, ByteArrayCharSequence chr) {
		byte[] key = new byte[tid.length() + chr.length()];
		System.arraycopy(tid.chars, tid.start, key, 0, tid.length());
		System.arraycopy(chr.chars, chr.start, key, tid.length(), chr.length());
		return key;
	}

    /**
     * Find the field index of the transcript id or return -1
     *
     * @param input the line
     * @return transcriptIDindex the field index if the transcript ID
     */
	private static int findTranscriptID(CharSequence input) {
    	int index = 0;
		Matcher m= SPLITTER_PATTERN.matcher(input);
	    for(int mCtr= 0;m.find();index= m.end(), mCtr++) {
	    	if (mCtr< 9) 
	    		continue;	    	
	    	CharSequence match = input.subSequence(index, m.start());
	    	Matcher m2= TRANSCRIPTID_PATTERN.matcher(match);
	    	if (m2.find())
	    		return mCtr+1;
	    }
	    return -1;
	}

    /**
     * Split the given line and extract the field specified
     *
     * @param input the source line
     * @param fieldNrs the field to extract
     * @return fields extracted fields
     */
	private static ByteArrayCharSequence[] find(CharSequence input, int[] fieldNrs) {
		int index = 0;
		ByteArrayCharSequence[] result= new ByteArrayCharSequence[fieldNrs.length];
		int mCtr= 0, fCtr= 0;
		Matcher m= SPLITTER_PATTERN.matcher(input);
        while(m.find()) {
        	if (fCtr< fieldNrs.length&& mCtr== fieldNrs[fCtr]) {
        		ByteArrayCharSequence match = (ByteArrayCharSequence) input.subSequence(index, m.start());
        		result[fCtr++]= match;
        	}
            index = m.end();
            ++mCtr;
        }
        // last
        if (fCtr< fieldNrs.length&& fieldNrs[fCtr]== mCtr)
        	result[fCtr]= (ByteArrayCharSequence) input.subSequence(index, input.length());
        return result;
	}

    /**
     * Create a map using the concatenated transcript ID and the chromosome name as key
     * and the minimal transcript position as value. This is later used to sort the
     * entries relative to their global position. We pass the field[] here, where the last position
     * is filled with the identified index of the transcript IDs. An exception is thrown if
     * a problem occurs while parsing the file.
     *
     * @param file the source file
     * @param fieldNrs the fields to check. Position 3 is set the to the transcriptID index
     * @return transcriptMap map from the concatenated transcriptID+chromosome as key and the minimal global position as value
     * @throws Exception in case of any error
     */
    Map<byte[], Integer> createTranscriptMap(File file, int[] fieldNrs) throws Exception{
        IOHandler io = IOHandlerFactory.createDefaultHandler();
        try{
            // guess the file separator from input
			eol = FileHelper.guessFileSep(file);
            // file size
			long size = file.length();
            /*
             estimate line count to estimate map size
             */
			long estLineCount = (size / 100);
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

            FileInputStream reader = new FileInputStream(file);
            io.addStream(reader);

            ByteArrayCharSequence cs= new ByteArrayCharSequence(1000);
            //int[] fieldNrs= new int[]{1,4,7,-1};	// ,-1 for gene
            int bytesRead = 0;
            int lineCounter = 0;
            while(io.readLine(reader, cs) != -1){
                lineCounter++;
				bytesRead += cs.length() + eol.length();

				if (cs.charAt(0)== '#')
					continue;
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

				ByteArrayCharSequence[] fields= find(cs, fieldNrs);	// 1,4,tid,gid
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
				byte[] key= encode(fields[3], fields[0]);
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
        }finally {
            io.close();
        }
    }

    /**
     * Sort the file and write the result to the given file
     *
     * @param f the source file
     * @param outFile the target file
     * @throws Exception in case of any errors
     */
	private void sortFile(File f, File outFile) throws Exception {
        IOHandler io = IOHandlerFactory.createDefaultHandler();
        OutputStream outStr = null;
        PipedOutputStream out = null;
        PipedInputStream in = null;
        BufferedWriter writer = null;
        Future sorterThread = null;


		try {
            // attentionAttention, the nrs have to be sorted - look in find()
            int[] fieldNrs= new int[]{0,3,6,-1};	// ,-1 for transcriptid
            Map<byte[], Integer> transcriptPositions = createTranscriptMap(f, fieldNrs);
			System.gc();
			ByteArrayCharSequence cs= new ByteArrayCharSequence(1000);


            outStr = new BufferedOutputStream(new FileOutputStream(outFile));
            out = new PipedOutputStream();
            in = new PipedInputStream(out);
            writer= new BufferedWriter(new OutputStreamWriter(out));

            sorterThread = Sorter.create(in, outStr, true)
                    .separator("\\s")
                    .field(0, false) // chr
                    .field(6, false) // strand
                    .field(new GlobalTranscriptPositionComparator(transcriptPositions, new int[]{fieldNrs[0], fieldNrs[3]}))
                            //.field(1, true) // tpos
                    .field(fieldNrs[3], false) // tid
                    .field(3, true).sortInBackground();// pos
            fieldNrs= new int[]{fieldNrs[0], fieldNrs[3]};	// 090901: start pos from map
			HashSet<String> setInvalidTx= new HashSet<String>();
			int nrInvalidLines= 0;

            FileInputStream fileInput = new FileInputStream(f);
            io.addStream(fileInput);

            while (io.readLine(fileInput,cs) != -1){
                // issue #56 make sure we skip empty lines and comment lines
                if (cs.length() == 0 || cs.charAt(0)== '#'){
                    continue;
                }

                ByteArrayCharSequence[] fields= find(cs, fieldNrs);	// 1,4,-1
				
                byte[] key= encode(fields[1], fields[0]);
                int val= transcriptPositions.get(key);
                if (val== Integer.MIN_VALUE) {
                    ++nrInvalidLines;
                    setInvalidTx.add(fields[1].toString());
                }else{
				    writer.write(cs.toString());
				    writer.write(eol);
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
		} finally {
            io.close();
            if(outStr != null) try {outStr.close();} catch (IOException e) {}
            if(out != null) try {out.close();} catch (IOException e) {}
            if(in != null) try {in.close();} catch (IOException e) {}
            if(writer != null) try {writer.close();} catch (IOException e) {}
            if(sorterThread != null){
                sorterThread.cancel(true);
            }
		}
	}

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
            ByteArrayCharSequence[] fields= find(cs, fieldNrs);	// 1,4,-1
            byte[] key1= encode(fields[1], fields[0]);
            Integer mapValue1 = map.get(key1);
            int v1= Math.abs(mapValue1);

            cs.clear();
            cs.append(o2);
            fields= find(cs, fieldNrs);	// 1,4,-1
            byte[] key2= encode(fields[1], fields[0]);
            Integer mapValue2 = map.get(key2);
            int v2= Math.abs(mapValue2);
            int result = v1 - v2;
            return result;
        }
    }

    /**
     * Sort the given GTF file, sort it and return the sorted file.
     *
     * @param input the source file
     * @return sorted the sorted file
     */
    public static File sort(File input) {
        File outFile = new File(MyFile.append(input.getAbsolutePath(), "_sorted"));
        sort(input, outFile);
        return outFile;
    }

    /**
     * Sort the given input file and write the sorted result to the given output file
     *
     * @param f the input file
     * @param output the output file
     */
    public static void sort(File f, File output) {
        if(f == null) throw new NullPointerException();
        if(output == null) throw new NullPointerException();
        try {
            new GFFSorter().sortFile(f, output);
        } catch (Exception e) {
            Log.error("Error while sorting " + f.getAbsolutePath() + " : " + e.getMessage(), e);
        }
    }
}
