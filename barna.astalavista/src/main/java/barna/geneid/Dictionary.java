package barna.geneid;

import barna.commons.log.Log;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/8/12
 * Time: 5:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class Dictionary {

    /* Dictionary definitions (hash)            */
    public static final int MAXENTRY= 97;
    public static final int MAXINFO= 100;
    public static final int NOTFOUND= -1;
    public static final char UNKNOWNAA= 'X';

    public static final String sBEGIN=    "Begin";
    public static final String sBEGINFWD= "Begin+";
    public static final String sBEGINRVS= "Begin-";
    public static final String sEND=      "End";
    public static final String sENDFWD=   "End+";
    public static final String sENDRVS=   "End-";

    Node[] T= new Node[GParam.MAXENTRY];
    int nextFree= 0;

    /**
     * Initializing the dictionary: hash table and counter of keys
     */
    public Dictionary() {
        for(int i= 0; i< GParam.MAXENTRY; i++)
            T[i] = null;
        nextFree = 0;
    }

    /**
     * Replicating the gene model rules for every isochore
     */
    static void shareGeneModel(GParam[] isochores, int nIsochores) {

        /* Original gene model is loaded in the first isochore */
        int nTypes = isochores[0].D.nextFree;

        /* Sharing global parameters for any isochore: i */
        for(int i= 1; i< nIsochores; i++) {
            isochores[i].D      = isochores[0].D;
            isochores[i].nc     = isochores[0].nc;
            isochores[i].ne     = isochores[0].ne;
            isochores[i].md     = isochores[0].md;
            isochores[i].Md     = isochores[0].Md;
            isochores[i].nclass = isochores[0].nclass;

            /* Copying info for every feature/exon type */
            for(int j= 0; j< nTypes; j++) {
                /* Frequency as left part on a rule (id. rule) */
                for(int k= 0; k< isochores[i].nc[j]; k++)
                    isochores[i].UC[j][k] = isochores[0].UC[j][k];

                /* Frequency as right part on a rule (id. rule) */
                for(int k= 0; k< isochores[i].ne[j]; k++)
                    isochores[i].DE[j][k] = isochores[0].DE[j][k];
            }

            /* Copy info for every class: rules requiring group checkpoint */
            for(int j= 0; j< isochores[0].nclass; j++)
                isochores[i].block[j]  = isochores[0].block[j];
        }
    }

    /**
     * Loading the gene model rules to build correct genes
     * Every rule is identified by the gm line where it has been found
     * Returns how many rules have been loaded right
     */
    int readGeneModel (BufferedReader buffy,
                        int[] nc, int[] ne,
                        int[][] UC,
                        int[][] DE,
                        long[] md, long[] Md,
                        int[] block) throws Exception {

        /* Identifier for feature (from dictionary) */
        //int a;

        /* Identifier for class (assembling rule) */
        //int nlines;

        /* Format for gene model rules:
 F1:F2:(...):Fn   F1:F2:(...):Fm dmin:dMax  [block] */

        /* Input lines from parameter file */
        int nlines= 0;
        String line= null;
        while((line= buffy.readLine())!= null) {
            /* For every line extracting features (upstream/downstream), */
            /* the minMax distances and the (optional) block */
            /* line number is the class/rule identifier */
            if ((!line.startsWith("#")) && line.trim().length()> 0) {
                /* 0. Backup the line to display errors */
                String lineCopy= line;

                /* 1. Splitting line into 4 parts: UC DE Dist block */
                String[] cols = line.split("(\\s)+");

                /* Three first columns are mandatory, last one is optional */
                if (cols.length< 3)
                    throw new RuntimeException("Wrong format in gene model rule:\n"+ lineCopy);

                /* 2. Processing upstream compatible features list */
                for (String t1 : cols[0].split(":")) {
                    /* Extracting and adding to the dictionary of types */
                    int a = setkeyDict(t1);
                    /* Save: this feature appeared in the UC list of this rule */
                    UC[a][nc[a]++]= nlines;
                }

                /* 3. Processing downstream equivalent features list */
                for (String t1 : cols[1].split(":")) {
                    /* Extracting and adding to the dictionary of types */
                    int a = setkeyDict(t1);
                    /* Save: this feature appeared in the DE list of this rule */
                    DE[a][ne[a]++]= nlines;
                }

                /* 4. Read the distances xx:yy [block] */
                /* a) minimum distance */
                StringTokenizer stoki= new StringTokenizer(cols[2], ":");
                if (!stoki.hasMoreTokens()) {
                    throw new RuntimeException("Wrong distance range (min) in gene model rule:\n"+ lineCopy);
                }
                md[nlines] = Long.parseLong(stoki.nextToken());

                /* b) maximum distance */
                if (!stoki.hasMoreTokens()) {
                    throw new RuntimeException("Wrong distance range (max) in gene model rule:\n"+ lineCopy);
                }

                /* To forget the DMAX requirement use the string SINFI */
                String t1= stoki.nextToken();
                if (t1.equals(GeneIDconstants.SINFI))
                    Md[nlines] = GeneIDconstants.INFI;
                else
                    Md[nlines] = Long.parseLong(t1);

                /* 5. Read the block record (to preserve group)... if exists */
                if (cols.length> 3)
                    block[nlines] = GeneIDconstants.BLOCK;
                else
                    block[nlines] = GeneIDconstants.NONBLOCK;

                nlines++;

            } /* End of if-comment */
        } /* Next rule to read */

        return(nlines);
    }

    /**
     * Fill in the Gene Model with artificial lines to build only one gene
     * Every rule is identified by the gm line where it has been found
     * Returns how many rules have been loaded right
     */
    int forceGeneModel (int nc[], int ne[],
                         int UC[][],
                         int DE[][],
                         long md[], long Md[],
                         int block[]) {


        Log.message("Force one gene prediction");

        /* Identifier for feature (from dictionary) */
        int a;

        /* Identifier for class (assembling rule) */
        int nlines = 0;

        /* 1. First+:Internal+     Internal+:Terminal+       20:40000 block */
        a = setkeyDict("First+");
        UC[a][nc[a]++] = nlines;
        a = setkeyDict("Internal+");
        UC[a][nc[a]++] = nlines;

        a = setkeyDict("Internal+");
        DE[a][ne[a]++]=nlines;
        a = setkeyDict("Terminal+");
        DE[a][ne[a]++]=nlines;

        md[nlines] = 20;
        Md[nlines] = 40000;
        block[nlines] = GeneIDconstants.BLOCK;
        nlines++;

        /* 2. Terminal-:Internal-  First-:Internal-          20:40000 blockr */
        a = setkeyDict("Terminal-");
        UC[a][nc[a]++] = nlines;
        a = setkeyDict("Internal-");
        UC[a][nc[a]++] = nlines;

        a = setkeyDict("First-");
        DE[a][ne[a]++]=nlines;
        a = setkeyDict("Internal-");
        DE[a][ne[a]++]=nlines;

        md[nlines] = 20;
        Md[nlines] = 40000;
        block[nlines] = GeneIDconstants.BLOCK;
        nlines++;

        /* 3. BEGIN+:BEGIN-      First+:Terminal-:Single+:Single-     0:Infinity */
        a = setkeyDict(sBEGINFWD);
        UC[a][nc[a]++] = nlines;
        a = setkeyDict(sBEGINRVS);
        UC[a][nc[a]++] = nlines;

        a = setkeyDict("First+");
        DE[a][ne[a]++]=nlines;
        a = setkeyDict("Terminal-");
        DE[a][ne[a]++]=nlines;
        a = setkeyDict("Single+");
        DE[a][ne[a]++]=nlines;
        a = setkeyDict("Single-");
        DE[a][ne[a]++]=nlines;

        md[nlines] = 0;
        Md[nlines] = GeneIDconstants.INFI;
        nlines++;

        /* 4. Terminal+:First-:Single+:Single-     END+:END-      0:Infinity*/
        a = setkeyDict("Terminal+");
        UC[a][nc[a]++] = nlines;
        a = setkeyDict("First-");
        UC[a][nc[a]++] = nlines;
        a = setkeyDict("Single+");
        UC[a][nc[a]++] = nlines;
        a = setkeyDict("Single-");
        UC[a][nc[a]++] = nlines;

        a = setkeyDict(sENDFWD);
        DE[a][ne[a]++]=nlines;
        a = setkeyDict(sENDRVS);
        DE[a][ne[a]++]=nlines;

        md[nlines] = 0;
        Md[nlines] = GeneIDconstants.INFI;
        nlines++;

        return nlines;
    }


    /**
     * Assign a number-key to the new word and store it
     */
    int setkeyDict(String s) {

        /* If this word exists at the dictionary don't insert */
        int key = getkeyDict(s);
        if (key == NOTFOUND) {
            int i = f(s);

            /* Alloc the new word */
            Node n= null;
            try {
                n= new Node();
            } catch (OutOfMemoryError e) {
                Log.error("Not enough memory: dictionary word");
                throw(e);
            }

            /* Filling the node */
            n.s= s;
            n.key = nextFree++;
            if(T[i] == null) {
                n.next = null;
                T[i] = n;
            } else {
                /* There are more nodes in this position: Colission */
                Node p = T[i];
                /* Insert at the begining of the list */
                T[i] = n;
                n.next = p;
            }
            key = n.key;
        }
        return(key);
    }

    /*
     * Hash Function:: String . Integer between 0..MAXENTRY-1
     */
    public int f(String s) {

        int total= 0;
        for(int i= 0; i< s.length(); i++)
            total = (i+1)* ((int) s.charAt(i)) + total;

        total = total % MAXENTRY;

        return total;
    }

    /*
     * Returns the key for the word request; 
     * NOTFOUND is Not found
     */
    public int getkeyDict(String s) {
        
        /* showDict(d); */
        int key = NOTFOUND;

        /* Computing hash function */
        int i = f(s);

        /* Empty list means not found */
        if(T[i]== null)
            key = NOTFOUND;
        else {
            /* There are more nodes in this position: run the list */
            Node p = T[i];
            /* Searching until the first position not used */
            boolean found= false;
            while(p!= null && !found ) {
                /* Same hash value: compare to see if it is the same string */
                if(s.equals(p.s)) {
                    found= true;
                    key = p.key;
                }
                p= p.next;
            }
            if(!found)
                key = NOTFOUND;
        }
        return(key);
    }

}
