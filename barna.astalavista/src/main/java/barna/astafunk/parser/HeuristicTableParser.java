package barna.astafunk.parser;

/**
 * @version 1.0
 * @autor vitorlc
 * @since 13/07/14
 */

import barna.astafunk.DP.Hit;
import barna.astafunk.DP.HitComparator;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

// TODO perform tests and check loops/logical iterations (vitor's reminder)


/**
 * This class describes a parser to gene ID X HMMER3.x --dombtblout output file (so-called heuristic table).
 * This table is a optional input for our algorithm, but it helps
 * to reduce the search time.
 *
 * The table format is
 *
 target name \t accession \t tlen \t query name \t accession \t qlen \t E-value \t score \t bias \t #  of  c-Evalue
 \t i-Evalue \t   score \t bias \t from  \t  to \t from \t to \t from \t to \t acc \t description of target
 */
public class HeuristicTableParser {

    /**
     * An BufferRead attribute. Reads text from a character-input stream, buffering characters so as to
     * provide for the efficient reading of characters, arrays, and lines.
     */

    private BufferedReader stdin;

    /**
     * Default constructor initializes BufferedReader object.
     * @param path Directory path to heuristic table
     * @throws FileNotFoundException
     */
    public HeuristicTableParser(String path) throws FileNotFoundException {
        this.stdin = new BufferedReader(new FileReader(path), 128);
    }

    /**
     * Parse HMMER --domtblout output file. Create a hash map
     * with gene/transcript ID and related ACC numbers from best non-overlapped hits.
     * @return Hash map with gene/transcript IDs as keys and Pfam ACC array as values
     * @throws IOException
     */
    public HashMap<String, List<Hit>> parse() throws IOException { //

        if(stdin==null)
            return null;

        String auxLine;
        HashMap<String,List<Hit>> referenceTable = new HashMap<String, List<Hit>>();

        while((auxLine = stdin.readLine()) != null){
            if(auxLine.startsWith("#") || auxLine.contains("n/a"))
                continue;

            String geneID;
            String [] line = auxLine.split("\\s++");
            geneID = line[0];
            String acc = line[4];

            float score = Float.parseFloat(line[13]);
            int startAlign = Integer.parseInt(line[17]);
            int endAlign = Integer.parseInt(line[18]);

            Hit newReferenceDomain = new Hit(null,null, acc, null, score,startAlign,endAlign,0,0,0,0,0,null);
            List <Hit> auxList = new ArrayList<Hit>();

            if(!referenceTable.containsKey(geneID)){
                auxList.add(newReferenceDomain);
                referenceTable.put(geneID, auxList);
            }else{
                auxList = referenceTable.get(geneID);
                auxList = updateDomainReference(auxList, newReferenceDomain);
                referenceTable.put(geneID,auxList);
            }
        }
        return referenceTable;
    }

    public List<Hit> updateDomainReference(List<Hit> domainList, Hit newHit) {
        boolean insert = true;
        List<Hit> hitAuxList = new ArrayList<Hit>();
        HitComparator hc = new HitComparator();

        for (Hit h : domainList) {
            if (newHit.getEndAlignment() < h.getStartAlignment()
                    || newHit.getStartAlignment() > h.getEndAlignment()) {
                continue;
            } else {
                if (h.getScore() >= newHit.getScore()) {
                    insert = false;
                    break;
                } else {
                    insert = true;
                    hitAuxList.add(h);
                }
            }
        }
        if (insert) {
            domainList.removeAll(hitAuxList);
            domainList.add(newHit);
        }
        hitAuxList.clear();

        return domainList;
    }

    public HashMap<String, List<String>> parseAll() throws IOException {

        if(stdin==null)
            return null;

        String auxLine;
        HashMap<String,List<String>> referenceTable = new HashMap<String, List<String>>();

        while((auxLine = stdin.readLine()) != null){
            if(auxLine.startsWith("#") || auxLine.contains("n/a"))
                continue;

            String transcriptID;
            String [] line = auxLine.split("\\s++");
            transcriptID = line[0].trim();
            String acc = line[4].trim();


            List <String> auxList;

            if(!referenceTable.containsKey(transcriptID)){
                auxList = new ArrayList<String>();
                auxList.add(acc);
                referenceTable.put(transcriptID, auxList);
            }else{
                auxList = referenceTable.get(transcriptID);
                // Auxilist contains newReference ACC?
                if(!auxList.contains(acc)){
                    auxList.add(acc);
                    referenceTable.put(transcriptID, auxList);
                }
            }
        }
        return referenceTable;


    }
}