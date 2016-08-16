package barna.astafunk.utils;

/**
 * This documentation describes the FastaReader class.
 *
 * @author Vitor Lima Coelho
 * @since 2013-11-22
 * @version 1
 */


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * The FastaReader class implements methods to read FASTA files.
 */
public final class FASTAReader {

    /**
     * Path to directory to FASTA file.
     */

    private String pathToFasta;

    /**
     * Hashmap of sequences and headers of FASTA file.
     */

    private LinkedHashMap<String, String> hashMapSequence;


    /**
     * Constructor of FastaReader class.
     * @param path Path to FASTA file.
     * @throws Exception
     */
    public FASTAReader(String path) throws Exception{

        this.pathToFasta = path;

        hashMapSequence();


    }

    public LinkedHashMap<String, String> readMultiFasta(File f) {

        LinkedHashMap<String, String> hash= null;
        try {
            hash= new LinkedHashMap<String, String>();
            BufferedReader buffy= new BufferedReader(new FileReader(f));
            String s, seq= null, id= null;
            while ((s= buffy.readLine())!= null) {

                if(s.startsWith(">")) {
                    if(id!= null)
                        hash.put(id, seq);
                    seq= "";
                    id=s.substring(1);
                } else
                    seq+= s;
            }
            if(id!= null)
                hash.put(id, seq);

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            return hash;
        }
    }

    /**
     * Generates hashmap for FASTA file.
     */
    public void hashMapSequence(){

        String path = getPathToFasta();

        System.out.print("\nReading FASTA file, please wait...");

        try{

            hashMapSequence = readMultiFasta(new File(path));
            //FastaReaderHelper.readFastaProteinSequence(new File(path));

            System.out.print("finished.\n");

        } catch (Exception ex) {

            System.err.println("ERROR in creating hashmap of FASTA file.");


            Logger.getLogger(FASTAReader.class.getName()).log(Level.SEVERE, null, ex);

        }

//        catch (CompoundNotFoundError e){
//            System.err.println("ERROR in reading FASTA file. Is FASTA format correct?");
//        }



    }

    /**
     * Show the sequences and headers.
     */

    public void showSequence(){

        for (Map.Entry<String, String> entry : hashMapSequence.entrySet() ) {

            // TODO
            //System.out.println( entry.getValue().getOriginalHeader() + "=" + entry.getValue().getSequenceAsString() );
        }

    }
    /**
     * Get the hashmap
     * @return The hashmap of sequences and headers
     */
    public LinkedHashMap<String, String> getHashMapSequence() {
        return hashMapSequence;
    }

    /**
     * Set the hashmap.
     * @param hashMapSequence Hashmap of sequences and headers.
     */
    public void setHashMapSequence(LinkedHashMap<String, String> hashMapSequence) {
        this.hashMapSequence = hashMapSequence;
    }
    /**
     * Get the path to FASTA file
     * @return The path to the FASTA file.
     */
    public String getPathToFasta(){
        return pathToFasta;
    }

    /**
     * Set the path to FASTA file
     * @param pathToFasta Path to FASTA file.
     */
    public void setPathToFasta(String pathToFasta) {
        this.pathToFasta = pathToFasta;
    }

}

