package barna.geneid;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/11/12
 * Time: 4:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class ExonGFF {

    static class PExonGFF extends ExonGFF {
             // typedef struct s_exonGFF *pexonGFF;
    }

    Site Acceptor;
    Site Donor;
    char[] Type= new char[GeneIDconstants.MAXTYPE];
    short Frame;
    short Remainder;
    char Strand;
    float PartialScore;
    float HSPScore;
    float R;
    float Score;
    PExonGFF PreviousExon;
    double GeneScore;
    char[] Group= new char[GeneIDconstants.MAXSTRING];
    int offset1;
    int offset2;
    short lValue;
    short rValue;
    short evidence;
    short selected;
    short three_prime_partial;
    short five_prime_partial;
}
