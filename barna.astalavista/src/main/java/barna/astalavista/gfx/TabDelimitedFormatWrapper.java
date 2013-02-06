package barna.astalavista.gfx;

import java.io.*;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 10/6/12
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class TabDelimitedFormatWrapper {
    String[][] table= null;

    /**
     * File name and extension.
     */
    protected String fName= null;

    /**
     * File path without seperator.
     */
    protected String fPath= null;

    /**
     * File last modification
     */
    protected long fLastModified= 0L;

    public TabDelimitedFormatWrapper(String absFilePath) {
        int p= absFilePath.length();
        while ((--p>= 0)&& (absFilePath.charAt(p)!= File.separatorChar))
            ; // count down
        if (p< 0) {
            this.fPath= ".";
            this.fName= absFilePath;
        } else {
            this.fPath= absFilePath.substring(0,p);
            this.fName= absFilePath.substring((p+1), absFilePath.length());
        }
    }

    public void write() throws Exception {
        if (table== null)
            return;
        try {
            BufferedWriter writer= new BufferedWriter(new FileWriter(fPath+
                    File.separator+ fName));
            for (int i = 0; i < table.length; i++) {
                if (table[i]== null)
                    continue;
                for (int j = 0; j < table[i].length; j++) {
                    if (table[i][j]== null)
                        continue;
                    writer.write(table[i][j]);
                    if (j< table[i].length- 1)
                        writer.write("\t");
                }
                writer.write("\n");
            }
            writer.flush(); writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public boolean isApplicable() throws Exception {
        System.err.println("implement..");
        return false;
    }

    public void read() throws Exception {
        Vector v= new Vector();
        try {
            BufferedReader buffy= new BufferedReader(new FileReader(fPath+ File.separator+ fName));
            int fCnt= -1;
            while (buffy.ready()) {
                String line= buffy.readLine();
                if (line.trim().length()< 1)
                    continue;
                StringTokenizer toki= new StringTokenizer(line, "\t");
                if (toki.countTokens()> fCnt)
                    fCnt= toki.countTokens();
//				else
//					if (fCnt!= toki.countTokens()) {
//						System.err.println("Invalid line, expected "+fCnt+" tokens, but found "
//								+toki.countTokens()+"\n"+line);
//						continue;
//					}
                String[] ln= new String[fCnt];
                int len= toki.countTokens();
                for (int i = 0; i < len; i++)
                    ln[i]= toki.nextToken();
                for (int i = len; i < fCnt; i++)
                    ln[i]= "MISS";
                v.add(ln);
            }
            buffy.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        table= (String[][]) Arrays.toField(v);
    }

    public String[][] getTable() {
        return table;
    }

    public String[] getColumn(int colNr) {
        if (table== null)
            return null;
        String[] col= new String[table.length];
        for (int i = 0; i < col.length; i++)
            col[i]= table[i][colNr];
        return col;
    }

    public void setTable(String[][] table) {
        this.table = table;
    }
}
