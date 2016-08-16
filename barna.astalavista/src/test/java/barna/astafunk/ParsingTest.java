package barna.astafunk;

import barna.astafunk.DP.Hit;
import barna.astafunk.HMM.ProfileHMM;
import barna.astafunk.parser.HMMParser;
import barna.astafunk.parser.HeuristicTableParser;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author vitorcoelho
 * @version 2
 * @since 16/11/15
 */
public class ParsingTest {


    @Test
    public void hmmerOutputTest() throws IOException {

        // read output hmmer file
        String hmmerOutputPath = new File(getClass().getResource("/example-hmmer-output.txt").getFile()).getAbsolutePath();
        HeuristicTableParser parseH = new HeuristicTableParser(hmmerOutputPath);
        HashMap<String, List<Hit>> heuristicHash = parseH.parse();
        printHeuristicTable(heuristicHash);
    }

    @Test
    public void createReferenceList() throws IOException {

        String hmmerOutputPath = new File(getClass().getResource("/example-hmmer-output.txt").getFile()).getAbsolutePath();

        HeuristicTableParser parseH = new HeuristicTableParser(hmmerOutputPath);
        HashMap<String, List<Hit>> heuristicHash = parseH.parse();
        List<String> referencePfamList = new ArrayList<String>();

        String [] gene_transcript_list = {"NM_001202431", "NM_002541", "NM_001328", "NM_022802", "NM_015328"};
        for(String id: gene_transcript_list) {
            if (heuristicHash.containsKey(id)) {
                for (Hit h : heuristicHash.get(id)) {
                    String acc = h.getAcc().split("\\.")[0];
                    if (!referencePfamList.contains(acc)) {
                        /*
                        List with Acc's of models that we will search
                         */
                        referencePfamList.add(acc);
                    }
                }
            }
        }
        System.out.println(referencePfamList);
    }

    public void printHeuristicTable(HashMap<String, List<Hit>> heuristicHash){
        for(Map.Entry<String, List<Hit>> entry: heuristicHash.entrySet()){
            System.out.println(entry.getKey() + "\t" + entry.getValue());
        }
    }

    @Test
    public void hmmFileTest() throws IOException {

        String hmmFilePath = new File(getClass().getResource("/test.hmm").getFile()).getAbsolutePath();

        HMMParser hmmParser = new HMMParser(hmmFilePath);

        List<ProfileHMM> list = hmmParser.parseAll();




    }
}
