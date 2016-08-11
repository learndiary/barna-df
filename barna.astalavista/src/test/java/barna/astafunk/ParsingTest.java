package barna.astafunk;

import barna.commons.log.Log;
import barna.func.DP.Hit;
import barna.func.HMM.ProfileHMM;
import barna.func.parser.HMMParser;
import barna.func.parser.HeuristicTableParser;
import barna.func.utils.FunkSettings;
import barna.model.Transcript;
import com.sun.corba.se.impl.encoding.OSFCodeSetRegistry;
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
        String hmmerOutputPath = "barna.astalavista/src/test/java/barna/astafunk/example-hmmer-output.txt";
        HeuristicTableParser parseH = new HeuristicTableParser(hmmerOutputPath);
        HashMap<String, List<Hit>> heuristicHash = parseH.parse();
        printHeuristicTable(heuristicHash);
    }

    @Test
    public void createReferenceList() throws IOException {

        String hmmerOutputPath = "barna.astalavista/src/test/java/barna/astafunk/example-hmmer-output.txt";

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

        String hmmFilePath = "barna.astalavista/src/test/java/barna/astafunk/test.hmm";

        HMMParser hmmParser = new HMMParser(hmmFilePath);

        List<ProfileHMM> list = hmmParser.parseAll();




    }
}
