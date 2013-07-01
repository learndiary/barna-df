package barna.astalavista;

import barna.model.ASEvent;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: Micha
 * Date: 15/06/13
 * Time: 17:53
 * To change this template use File | Settings | File Templates.
 */
public class EventHandler {

    public static void main(String[] args) {
        // "1-2^3-6^,5(7^"
        String pattern= "1-2^3-6^10-,5(7^8-9)";
        reformatEvent(pattern);
    }

    static void reformatEvent(String pattern) {

        System.out.println("pattern "+ pattern);
        String[] variants= pattern.split(",");

        // (1) find trim points
        int trim5= -1, trim3= -1;
        for(int i= 0; i< variants.length; ++i) {

            // 5'
            int p= 0;
            while(Character.isDigit(variants[i].charAt(++p)));
            if (variants[i].charAt(p)== '(') {
                int x= Integer.parseInt(variants[i].substring(0,p));
                trim5= Math.max(x, trim5); // max actually not necessary, if sorted
            }

            // 3'
            if (variants[i].charAt(variants[i].length()- 1)== ')') {
                p= variants[i].length()- 1;
                while(Character.isDigit(variants[i].charAt(--p)));
                int x= Integer.parseInt(variants[i].substring(p+ 1,variants[i].length()- 1));
                trim3= Math.max(x, trim3); // min actually not necessary, if sorted
            }
        }

        // (2) trim
        String[] flank5= new String[variants.length];
        String[] flank3= new String[variants.length];
        String[] var= new String[variants.length];
        for(int i= 0; i< variants.length; ++i) {
            String[] positions= variants[i].split("[\\(\\)\\^\\-]");
            String[] types= variants[i].split("[0-9]");
            types= clearEmptyElements(types);
            int[] pos= new int[positions.length];
            for(int j= 0; j< positions.length; ++j) {
                pos[j]= Integer.parseInt(positions[j]);
            }
            int p= 0, q= positions.length;
            if (trim5> 0) {
                p= Arrays.binarySearch(pos, trim5);
                if (p< 0)
                    p= -(p+ 1)- 1;
                for(int j=0; j<= p; ++j) {
                    flank5[i]= append(flank5[i], positions[j]+ types[j]);
                }
            }
            if (trim3> 0) {
                q= Arrays.binarySearch(pos, trim3);
                if (q< 0)
                    q= -(q+ 1);
                for(int j= q; j< positions.length; ++j) {
                    flank3[i]= append(flank3[i], positions[j]+ types[j]);
                }
            }
            for(int j=p+ 1; j< q; ++j) {
                var[i]= append(var[i], Integer.toString(pos[j]- trim5)+ types[j]);
            }
        }

        // build new Event
        ASEvent event= new ASEvent();
        System.out.println("flank5 "+ flank5[0]+ ","+ flank5[1]);
        System.out.println("flank3 "+ flank3[0]+ ","+ flank3[1]);
        System.out.println("var "+ var[0]+ ","+ var[1]);
    }

    private static String append(String s, String sfx) {
        if (s== null)
            return sfx;
        return s+ sfx;
    }

    private static String[] clearEmptyElements(String[] a) {
        int ee= 0;
        for (int i = 0; i < a.length; i++) {
            if(a[i].trim().length()== 0)
                ++ee;
        }

        if (ee== 0)
            return a;

        String[] b= new String[a.length- ee];
        ee= 0;
        for (int i = 0; i < a.length; i++) {
            if(a[i].trim().length()== 0) {
                ++ee;
            } else {
                b[i- ee]= a[i];
            }
        }

        return b;
    }

}
